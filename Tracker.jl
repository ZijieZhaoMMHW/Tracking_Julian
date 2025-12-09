"""
    Tracker.jl

This module applies all Ocetrac methods for labeling objects in 3 dimensions.
Translated from Python to Julia.
"""

using Images
using ImageMorphology
using ImageSegmentation
using Statistics
using LinearAlgebra

"""
    apply_mask(binary_images, mask)

Apply a binary mask to the images.

# Arguments
- `binary_images`: Binary image array
- `mask`: Binary mask (1 = valid, 0 = background)

# Returns
- Masked binary images
"""
function apply_mask(binary_images, mask)
    binary_images_with_mask = copy(binary_images)
    binary_images_with_mask[mask .!= 1] .= 0
    return binary_images_with_mask
end


"""
    Tracker

Tracker object for applying Ocetrac tracking algorithm.

# Fields
- `data`: Input data array (3D: time, y, x)
- `mask`: Binary mask (1 = valid, 0 = ignore)
- `radius`: Size of structuring element for morphological operations
- `min_size_quartile`: Quantile threshold for minimum object size (0-1)
- `timedim`: Index of time dimension (default: 1)
- `xdim`: Index of x dimension (default: 3)
- `ydim`: Index of y dimension (default: 2)
- `positive`: Whether to track positive values (true) or negative (false)
"""
mutable struct Tracker
    data::Array{<:Real, 3}
    mask::Array{<:Real, 3}
    radius::Int
    min_size_quartile::Float64
    timedim::Int
    xdim::Int
    ydim::Int
    positive::Bool
    
    function Tracker(data, mask, radius, min_size_quartile; 
                     timedim=1, xdim=3, ydim=2, positive=true)
        
        # Validate inputs
        if ndims(data) != 3
            error("Ocetrac currently only supports 3D arrays with dimensions (time, y, x)")
        end
        
        if size(data) != size(mask)
            error("Data and mask must have the same dimensions")
        end
        
        if radius < 1
            error("radius must be an integer greater than or equal to 1")
        end
        
        if min_size_quartile < 0 || min_size_quartile > 1
            error("min_size_quartile must be between 0 and 1")
        end
        
        new(data, mask, radius, min_size_quartile, timedim, xdim, ydim, positive)
    end
end


"""
    track(tracker::Tracker)

Label and track image features.

# Arguments
- `tracker`: Tracker object containing data and parameters

# Returns
- `labels`: Integer labels of the connected regions (3D array)
- `metadata`: Dictionary containing tracking statistics

# Description
The tracking algorithm:
1. Converts data to binary based on threshold
2. Applies morphological closing and opening
3. Applies mask to filter valid regions
4. Filters objects by area threshold
5. Labels connected components
6. Wraps labels for periodic boundary conditions
"""
function track(tracker::Tracker)
    
    if all(tracker.mask .== 0)
        error("Found only zeros in `mask` input. The mask should indicate valid regions with values of 1")
    end
    
    # Convert data to binary, define structuring element, and perform morphological operations
    println("Performing morphological operations...")
    binary_images = morphological_operations(tracker)
    
    # Apply mask
    println("Applying mask...")
    binary_images_with_mask = apply_mask(binary_images, tracker.mask)
    
    # Filter area
    println("Filtering by area...")
    area, min_area, binary_labels, N_initial = filter_area(tracker, binary_images_with_mask)
    
    # Label objects
    println("Labeling connected components...")
    labels = label_components(binary_labels, trues(3, 3, 3))  # 3D connectivity
    
    # Wrap labels for periodic boundary
    println("Wrapping labels...")
    grid_res = abs(tracker.data[1, 1, 2] - tracker.data[1, 1, 1])
    x_range = size(tracker.data, tracker.xdim)
    
    # Check if data spans full longitude range (360 degrees)
    if x_range * grid_res >= 360 - grid_res
        labels_wrapped, N_final = wrap_labels(labels)
    else
        labels_wrapped = labels
        N_final = maximum(labels)
    end
    
    # Replace 0 with NaN for output
    labels_output = Float64.(labels_wrapped)
    labels_output[labels_output .== 0] .= NaN
    
    # Calculate metadata
    sum_tot_area = sum(area)
    
    reject_area = area[area .<= min_area]
    sum_reject_area = sum(reject_area)
    percent_area_reject = sum_reject_area / sum_tot_area
    
    accept_area = area[area .> min_area]
    sum_accept_area = sum(accept_area)
    percent_area_accept = sum_accept_area / sum_tot_area
    
    metadata = Dict(
        "initial_objects_identified" => N_initial,
        "final_objects_tracked" => N_final,
        "radius" => tracker.radius,
        "size_quantile_threshold" => tracker.min_size_quartile,
        "min_area" => min_area,
        "percent_area_reject" => percent_area_reject,
        "percent_area_accept" => percent_area_accept
    )
    
    println("Initial objects identified: \t", N_initial)
    println("Final objects tracked: \t\t", N_final)
    
    return labels_output, metadata
end


"""
    morphological_operations(tracker::Tracker)

Convert data to binary and perform morphological closing then opening.

# Arguments
- `tracker`: Tracker object

# Returns
- Binary images after morphological operations
"""
function morphological_operations(tracker::Tracker)
    
    # Convert to binary based on positive/negative threshold
    if tracker.positive
        bitmap_binary = tracker.data .> 0
    else
        bitmap_binary = tracker.data .< 0
    end
    
    # Convert boolean to integer (0 or 1)
    bitmap_binary = Int.(bitmap_binary)
    
    # Define structuring element (disk)
    se = strel_disk(tracker.radius)
    
    # Process each time slice
    nt, ny, nx = size(tracker.data)
    mo_binary = similar(bitmap_binary)
    
    for t in 1:nt
        slice = bitmap_binary[t, :, :]
        
        # Pad for wrap-around (periodic boundary)
        diameter = tracker.radius * 2
        padded = wrap_pad(slice, diameter)
        
        # Morphological operations
        if tracker.radius == 1
            # For radius=1, structuring element is single pixel, no-op
            processed = padded
        else
            # Closing then opening
            closed = closing(padded, se)
            processed = opening(closed, se)
        end
        
        # Remove padding
        mo_binary[t, :, :] = processed[diameter+1:end-diameter, diameter+1:end-diameter]
    end
    
    return mo_binary
end


"""
    strel_disk(radius::Int)

Create a disk-shaped structuring element.

# Arguments
- `radius`: Radius of the disk

# Returns
- Boolean matrix representing the structuring element
"""
function strel_disk(radius::Int)
    diameter = radius * 2
    x = -radius:radius
    X = repeat(x', length(x), 1)
    Y = repeat(x, 1, length(x))
    R = X.^2 + Y.^2
    se = R .< radius^2
    return se
end


"""
    wrap_pad(array::Matrix, pad_width::Int)

Pad array with wrap-around (periodic) boundary conditions.

# Arguments
- `array`: 2D array to pad
- `pad_width`: Number of pixels to pad on each side

# Returns
- Padded array
"""
function wrap_pad(array::Matrix, pad_width::Int)
    ny, nx = size(array)
    padded = zeros(eltype(array), ny + 2*pad_width, nx + 2*pad_width)
    
    # Center
    padded[pad_width+1:pad_width+ny, pad_width+1:pad_width+nx] = array
    
    # Edges (wrap around)
    # Left
    padded[pad_width+1:pad_width+ny, 1:pad_width] = array[:, end-pad_width+1:end]
    # Right
    padded[pad_width+1:pad_width+ny, end-pad_width+1:end] = array[:, 1:pad_width]
    # Top
    padded[1:pad_width, pad_width+1:pad_width+nx] = array[end-pad_width+1:end, :]
    # Bottom
    padded[end-pad_width+1:end, pad_width+1:pad_width+nx] = array[1:pad_width, :]
    
    # Corners
    # Top-left
    padded[1:pad_width, 1:pad_width] = array[end-pad_width+1:end, end-pad_width+1:end]
    # Top-right
    padded[1:pad_width, end-pad_width+1:end] = array[end-pad_width+1:end, 1:pad_width]
    # Bottom-left
    padded[end-pad_width+1:end, 1:pad_width] = array[1:pad_width, end-pad_width+1:end]
    # Bottom-right
    padded[end-pad_width+1:end, end-pad_width+1:end] = array[1:pad_width, 1:pad_width]
    
    return padded
end


"""
    filter_area(tracker::Tracker, binary_images)

Filter objects by area using quantile threshold.

# Arguments
- `tracker`: Tracker object
- `binary_images`: Binary image array

# Returns
- `area`: Array of object areas
- `min_area`: Minimum area threshold
- `binary_labels`: Binary labels after filtering
- `N_initial`: Initial number of objects
"""
function filter_area(tracker::Tracker, binary_images)
    
    # Label each time slice independently
    nt, ny, nx = size(binary_images)
    labels = zeros(Int, nt, ny, nx)
    
    max_id = 0
    for t in 1:nt
        slice = binary_images[t, :, :]
        slice_labels = label_components(slice .> 0, trues(3, 3))  # 2D 8-connectivity
        
        # Make labels consecutive across time
        if t > 1
            max_id = maximum([max_id, maximum(labels[t-1, :, :])])
            slice_labels = slice_labels .+ max_id
        end
        
        labels[t, :, :] = slice_labels
    end
    
    # Wrap labels
    labels[labels .== 0] .= 0
    labels_wrapped, N_initial = wrap_labels(labels)
    
    # Calculate area of each object
    props = extract_region_props(labels_wrapped)
    
    if isempty(props)
        error("No objects were detected. Try changing radius or min_size_quartile parameters.")
    end
    
    labelprops = [p.label for p in props]
    area = [p.area for p in props]
    
    # Calculate minimum area threshold
    min_area = quantile(area, tracker.min_size_quartile)
    println("Minimum area: $min_area")
    
    # Keep only objects larger than threshold
    keep_labels = labelprops[area .>= min_area]
    keep_where = [l in keep_labels for l in labels_wrapped]
    out_labels = copy(labels_wrapped)
    out_labels[.!keep_where] .= 0
    
    # Convert to binary
    binary_labels = out_labels .> 0
    binary_labels = Int.(binary_labels)
    
    return area, min_area, binary_labels, N_initial
end


"""
    RegionProp

Structure to hold region properties (similar to skimage regionprops).
"""
struct RegionProp
    label::Int
    area::Int
    centroid::Tuple{Float64, Float64, Float64}
    bbox::Tuple{Int, Int, Int, Int, Int, Int}
end


"""
    extract_region_props(labels::Array{Int, 3})

Extract region properties from labeled array.

# Arguments
- `labels`: Labeled array (3D)

# Returns
- Array of RegionProp structures
"""
function extract_region_props(labels::Array{Int, 3})
    unique_labels = unique(labels)
    filter!(x -> x > 0, unique_labels)  # Remove background (0)
    
    props = RegionProp[]
    
    for label_val in unique_labels
        indices = findall(labels .== label_val)
        area = length(indices)
        
        # Calculate centroid
        coords = [[idx[1], idx[2], idx[3]] for idx in indices]
        centroid_t = mean([c[1] for c in coords])
        centroid_y = mean([c[2] for c in coords])
        centroid_x = mean([c[3] for c in coords])
        centroid = (centroid_t, centroid_y, centroid_x)
        
        # Calculate bounding box
        t_coords = [idx[1] for idx in indices]
        y_coords = [idx[2] for idx in indices]
        x_coords = [idx[3] for idx in indices]
        
        bbox = (minimum(t_coords), minimum(y_coords), minimum(x_coords),
                maximum(t_coords), maximum(y_coords), maximum(x_coords))
        
        push!(props, RegionProp(label_val, area, centroid, bbox))
    end
    
    return props
end


"""
    wrap_labels(labels::Array{Int, 3})

Impose periodic boundary and wrap labels across longitude boundaries.

# Arguments
- `labels`: Labeled array (3D)

# Returns
- `labels_wrapped`: Labels with wrapped boundaries
- `N`: Total number of unique labels
"""
function wrap_labels(labels::Array{Int, 3})
    labels = copy(labels)
    
    # Get first and last columns
    first_column = labels[:, :, 1]
    last_column = labels[:, :, end]
    
    unique_first = unique(first_column[first_column .> 0])
    
    # Iterate over unique values in first column
    for first_val in unique_first
        # Find locations in first column
        first_locs = findall(first_column .== first_val)
        
        # Get corresponding values in last column
        for loc in first_locs
            last_val = last_column[loc[1], loc[2]]
            
            if last_val > 0
                # Find all instances of this label in last column
                last_locs = findall(last_column .== last_val)
                
                for last_loc in last_locs
                    # Get all labels that match
                    bad_labels = unique(labels[:, last_loc[1], last_loc[2]])
                    filter!(x -> x > 0, bad_labels)
                    
                    # Replace all instances with first column value
                    for bad_label in bad_labels
                        labels[labels .== bad_label] .= first_val
                    end
                end
            end
        end
    end
    
    # Relabel consecutively
    unique_vals = sort(unique(labels))
    labels_wrapped = zeros(Int, size(labels))
    
    for (new_label, old_label) in enumerate(unique_vals)
        labels_wrapped[labels .== old_label] .= new_label
    end
    
    # Recalculate total number of labels
    N = maximum(labels_wrapped)
    
    return labels_wrapped, N
end


"""
    save_labels(labels, filename)

Save labeled array to file.

# Arguments
- `labels`: Labeled array
- `filename`: Output filename (supports .jld2, .nc, etc.)
"""
function save_labels(labels, filename)
    # This can be extended with JLD2, NCDatasets, etc.
    println("Saving labels to $filename")
    # Add actual save functionality as needed
end


# ============================================================================
# Example Usage
# ============================================================================

"""
# Basic example
using Tracker

# Load your data (3D array: time, y, x)
# data = ... load your data
# mask = ones(size(data))  # or create appropriate mask

# Create tracker
tracker = Tracker(
    data, 
    mask, 
    radius=2,               # structuring element radius
    min_size_quartile=0.5,  # keep objects larger than median size
    positive=true           # track positive anomalies
)

# Run tracking
labels, metadata = track(tracker)

# Access results
println("Tracked ", metadata["final_objects_tracked"], " objects")

# Save results
# save_labels(labels, "tracked_objects.jld2")
"""

