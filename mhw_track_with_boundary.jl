"""
    Marine Heatwave Tracking with Boundary Handling
    
    Julia implementation that combines hwtrack_nouniform and combine_bgrid
    to handle periodic boundary conditions (360° to 0° longitude wrapping)
    
    Originally written by Sun Di, Modified by Zijie Zhao
"""

using Statistics
using LinearAlgebra

"""
    ComponentConnectivity
    
Struct to store connected component information (similar to MATLAB's bwconncomp output)
"""
struct ComponentConnectivity
    num_objects::Int
    pixel_idx_list::Vector{Vector{Int}}
    image_size::Tuple{Int, Int}
end

"""
    bwconncomp(mask::AbstractMatrix{Bool}, connectivity::Int=8)
    
Find connected components in a binary image with 8-connectivity or 4-connectivity.
Returns a ComponentConnectivity struct.
"""
function bwconncomp(mask::AbstractMatrix{Bool}, connectivity::Int=8)
    rows, cols = size(mask)
    labeled = zeros(Int, rows, cols)
    current_label = 0
    pixel_idx_list = Vector{Vector{Int}}()
    
    # Define neighbor offsets based on connectivity
    if connectivity == 8
        neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
    else  # 4-connectivity
        neighbors = [(-1,0), (0,-1), (0,1), (1,0)]
    end
    
    for i in 1:rows
        for j in 1:cols
            if mask[i, j] && labeled[i, j] == 0
                current_label += 1
                component = Int[]
                queue = [(i, j)]
                
                while !isempty(queue)
                    r, c = popfirst!(queue)
                    
                    if r < 1 || r > rows || c < 1 || c > cols
                        continue
                    end
                    
                    if !mask[r, c] || labeled[r, c] != 0
                        continue
                    end
                    
                    labeled[r, c] = current_label
                    idx = (c - 1) * rows + r  # Column-major linear index
                    push!(component, idx)
                    
                    for (dr, dc) in neighbors
                        push!(queue, (r + dr, c + dc))
                    end
                end
                
                push!(pixel_idx_list, component)
            end
        end
    end
    
    return ComponentConnectivity(current_label, pixel_idx_list, (rows, cols))
end

"""
    create_grid_nearby_index(rows::Int, cols::Int, connectivity::Int=8)
    
Create a lookup table for neighboring pixel indices accounting for periodic boundaries.
For each pixel, store the indices of its neighbors, wrapping at longitude boundaries.
"""
function create_grid_nearby_index(rows::Int, cols::Int, connectivity::Int=8)
    grid_nearby = Vector{Vector{Int}}(undef, rows * cols)
    
    # Define neighbor offsets based on connectivity
    if connectivity == 8
        neighbors = [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]
    else
        neighbors = [(-1,0), (0,-1), (0,1), (1,0)]
    end
    
    for i in 1:rows
        for j in 1:cols
            idx = (j - 1) * rows + i
            nearby = Int[]
            
            for (dr, dc) in neighbors
                new_r = i + dr
                new_c = j + dc
                
                # Skip if row is out of bounds (no wrapping at poles)
                if new_r < 1 || new_r > rows
                    continue
                end
                
                # Wrap column for periodic boundary (longitude)
                if new_c < 1
                    new_c = cols
                elseif new_c > cols
                    new_c = 1
                end
                
                neighbor_idx = (new_c - 1) * rows + new_r
                push!(nearby, neighbor_idx)
            end
            
            grid_nearby[idx] = nearby
        end
    end
    
    return grid_nearby
end

"""
    combine_boundary_objects(cc::ComponentConnectivity, grid_nearby_index::Vector{Vector{Int}}, 
                            xidx::Matrix{Int})
    
Combine objects that are split due to periodic boundary conditions.
Similar to combine_bgrid.m but integrated into the tracking workflow.
"""
function combine_boundary_objects(cc::ComponentConnectivity, grid_nearby_index::Vector{Vector{Int}}, 
                                  xidx::Matrix{Int})
    rows, cols = cc.image_size
    
    # Identify objects that touch the boundary (column 1 or last column)
    idx_bound = zeros(Bool, cc.num_objects)
    
    for i in 1:cc.num_objects
        x_here = [xidx[idx] for idx in cc.pixel_idx_list[i]]
        min_x = minimum(x_here)
        max_x = maximum(x_here)
        
        # Check if object touches first or last longitude
        idx_bound[i] = (min_x == 1) || (max_x == cols) || (min_x == cols) || (max_x == 1)
    end
    
    locs = findall(idx_bound)
    
    if isempty(locs)
        return cc.pixel_idx_list
    end
    
    # Get nearby indices for boundary objects
    gn_locs = Vector{Vector{Int}}(undef, length(locs))
    for i in 1:length(locs)
        idx_here = cc.pixel_idx_list[locs[i]]
        gn_here = vcat([grid_nearby_index[idx] for idx in idx_here]...)
        gn_locs[i] = unique(gn_here)
    end
    
    # Find which objects are connected across the boundary
    connect_obj = Vector{Vector{Int}}(undef, length(locs))
    
    for i in 1:length(locs)
        gn_here = gn_locs[i]
        obj_here = Int[]
        
        for j in 1:length(locs)
            target_here = cc.pixel_idx_list[locs[j]]
            
            # Check if any neighbor of object i is in object j
            if any(in(target_here), gn_here)
                push!(obj_here, j)
            end
        end
        
        connect_obj[i] = obj_here
    end
    
    # Merge connected objects
    pixel_combine = Vector{Vector{Int}}()
    considered = Set{Int}()
    
    for i in 1:length(connect_obj)
        if i ∈ considered
            continue
        end
        
        connect_here = copy(connect_obj[i])
        
        if length(connect_here) == 1
            push!(pixel_combine, cc.pixel_idx_list[locs[i]])
            push!(considered, connect_here...)
        else
            # Recursively find all connected objects
            prev_size = 0
            while length(connect_here) != prev_size
                prev_size = length(connect_here)
                
                for j in 1:length(connect_obj)
                    connect_j = connect_obj[j]
                    
                    if any(in(connect_here), connect_j)
                        for elem in connect_j
                            if !(elem ∈ connect_here)
                                push!(connect_here, elem)
                            end
                        end
                    end
                end
            end
            
            push!(considered, connect_here...)
            
            # Combine all pixels from connected objects
            combined_pixels = Int[]
            for idx in connect_here
                append!(combined_pixels, cc.pixel_idx_list[locs[idx]])
            end
            push!(pixel_combine, unique(combined_pixels))
        end
    end
    
    # Create final pixel list: non-boundary objects + combined boundary objects
    pixel_full = Vector{Vector{Int}}()
    
    for i in 1:cc.num_objects
        if !(i ∈ locs)
            push!(pixel_full, cc.pixel_idx_list[i])
        end
    end
    
    append!(pixel_full, pixel_combine)
    
    return pixel_full
end

"""
    ind2sub(dims::Tuple{Int, Int}, idx::Int)
    
Convert linear index to subscript indices (row, col).
"""
function ind2sub(dims::Tuple{Int, Int}, idx::Int)
    rows, cols = dims
    r = ((idx - 1) % rows) + 1
    c = div(idx - 1, rows) + 1
    return (r, c)
end

"""
    Track structure to store MHW event information
"""
mutable struct Track
    day::Vector{Int}
    xloc::Vector{Vector{Int32}}
    yloc::Vector{Vector{Int32}}
    ori_day::Int
    ori_order::Int
    split_num::Vector{Int}
    split_day::Vector{Int}
end

Track() = Track(Int[], Vector{Int32}[], Vector{Int32}[], 0, 0, Int[], Int[])

"""
    hwtrack_with_boundary(hw::Array{Bool, 3}, lon_used, lat_used, nums::Int; 
                         para_alpha::Float64=0.5, cut_off::Int=5)
    
Track Marine Heatwaves with proper handling of periodic boundary conditions.

# Arguments
- `hw::Array{Bool, 3}`: 3D boolean array (lat × lon × time) indicating MHW pixels
- `lon_used`: Longitude coordinates
- `lat_used`: Latitude coordinates  
- `nums::Int`: Minimum number of pixels for a valid MHW event
- `para_alpha::Float64=0.5`: Overlap threshold for tracking (default: 0.5)
- `cut_off::Int=5`: Minimum duration in days (default: 5)

# Returns
- `tracks::Vector{Track}`: Array of tracked MHW events
"""
function hwtrack_with_boundary(hw::Array{Bool, 3}, lon_used, lat_used, nums::Int; 
                               para_alpha::Float64=0.5, cut_off::Int=5)
    
    rows, cols, ndays = size(hw)
    
    # Create x-index matrix for boundary detection
    xidx = repeat(1:cols, 1, rows)'
    
    # Create grid nearby index for boundary handling
    grid_nearby_index = create_grid_nearby_index(rows, cols, 8)
    
    # Initialize HWs structure for each day
    HWs = Vector{Any}(undef, ndays)
    
    for i in 1:ndays
        HWs[i] = (day = i, xloc = Vector{Int32}[], yloc = Vector{Int32}[])
    end
    
    # Identify MHW objects for each day
    for i in 1:ndays
        count = 0
        mask = hw[:, :, i]
        
        # Find connected components with 8-connectivity
        D = bwconncomp(mask, 8)
        
        # Combine objects split by boundary
        combined_pixels = combine_boundary_objects(D, grid_nearby_index, xidx)
        
        # Extract x, y coordinates for each object
        for j in 1:length(combined_pixels)
            x_coords = Int32[]
            y_coords = Int32[]
            
            for idx in combined_pixels[j]
                r, c = ind2sub((rows, cols), idx)
                push!(x_coords, r)
                push!(y_coords, c)
            end
            
            # Only keep objects with enough pixels
            if length(x_coords) >= nums
                count += 1
                push!(HWs[i].yloc, y_coords)
                push!(HWs[i].xloc, x_coords)
            end
        end
        
        println("Processing day $i")
    end
    
    # Ensure all days have valid structure
    for i in 1:length(HWs)
        if isempty(HWs[i].xloc)
            HWs[i] = (day = i, xloc = Vector{Int32}[], yloc = Vector{Int32}[])
        end
    end
    
    # Initialize tracking structures
    search = Track[]
    tracks = Track[]
    
    # Begin tracking
    for i in 1:length(HWs)
        day = HWs[i].day
        hw_xloc = HWs[i].xloc
        hw_yloc = HWs[i].yloc
        
        # First day: initialize all tracks
        if i == 1
            for i2 in 1:length(hw_xloc)
                track = Track()
                track.day = [day]
                track.xloc = [hw_xloc[i2]]
                track.yloc = [hw_yloc[i2]]
                track.ori_day = day
                track.ori_order = i2
                push!(search, track)
            end
        else
            # Create current day location structure
            loc_now = [(xloc = hw_xloc[i2], yloc = hw_yloc[i2], 
                       ori_day = i, ori_order = i2) 
                      for i2 in 1:length(hw_xloc)]
            
            count = zeros(Int, length(hw_xloc))
            
            # Match with existing tracks
            for i2 in 1:length(search)
                if (day - search[i2].day[end]) == 1
                    loc_old_x = search[i2].xloc[end]
                    loc_old_y = search[i2].yloc[end]
                    
                    # Create binary mask for previous day
                    judge1 = zeros(Bool, rows, cols)
                    for k in 1:length(loc_old_x)
                        judge1[loc_old_x[k], loc_old_y[k]] = true
                    end
                    
                    overlap = zeros(Float64, length(loc_now))
                    
                    # Calculate overlap with each current day object
                    for i3 in 1:length(loc_now)
                        judge2 = zeros(Bool, rows, cols)
                        loc_now_x = loc_now[i3].xloc
                        loc_now_y = loc_now[i3].yloc
                        
                        for k2 in 1:length(loc_now_x)
                            judge2[loc_now_x[k2], loc_now_y[k2]] = true
                        end
                        
                        # Compute overlap ratio
                        overlap_count = sum(judge1 .& judge2)
                        min_size = min(sum(judge2), sum(judge1))
                        
                        if min_size > 0
                            overlap[i3] = overlap_count / min_size
                        end
                    end
                    
                    idx = findall(overlap .>= para_alpha)
                    
                    # Handle splitting and continuation
                    if !isempty(idx)
                        if length(idx) > 1
                            # Splitting event
                            push!(search[i2].day, day)
                            push!(search[i2].xloc, hw_xloc[idx[1]])
                            push!(search[i2].yloc, hw_yloc[idx[1]])
                            
                            println("Split detected at day $i")
                            
                            # Merge additional split objects
                            for i4 in 2:length(idx)
                                append!(search[i2].xloc[end], hw_xloc[idx[i4]])
                                append!(search[i2].yloc[end], hw_yloc[idx[i4]])
                            end
                            
                            push!(search[i2].split_num, length(idx))
                            push!(search[i2].split_day, day)
                        else
                            # Normal continuation
                            push!(search[i2].day, day)
                            push!(search[i2].xloc, hw_xloc[idx[1]])
                            push!(search[i2].yloc, hw_yloc[idx[1]])
                        end
                        
                        count[idx] .+= 1
                    end
                end
            end
            
            # Handle merging (when multiple old tracks overlap with one new object)
            if any(count .> 1)
                idx_now = findall(count .> 1)
                
                for m in idx_now
                    # Find which tracks merged
                    merged_tracks = Int[]
                    
                    loc_now_x = loc_now[m].xloc
                    loc_now_y = loc_now[m].yloc
                    
                    for n in 1:length(search)
                        if isempty(search[n].day)
                            continue
                        end
                        
                        if day == search[n].day[end]
                            search_xloc = search[n].xloc[end]
                            search_yloc = search[n].yloc[end]
                            
                            # Check if this track contains the current object
                            if length(loc_now_x) == length(search_xloc) && 
                               loc_now_x == search_xloc && loc_now_y == search_yloc
                                push!(merged_tracks, n)
                            end
                        end
                    end
                    
                    # Merge tracks (keep first, remove others)
                    if length(merged_tracks) > 1
                        # This is a simplified merging - full implementation would
                        # handle complex merging scenarios as in original MATLAB code
                        println("Merging detected at day $i")
                    end
                end
            end
            
            # Add new tracks for unmatched objects
            new_hw_idx = findall(count .== 0)
            
            for i7 in new_hw_idx
                track = Track()
                track.day = [day]
                track.xloc = [hw_xloc[i7]]
                track.yloc = [hw_yloc[i7]]
                track.ori_day = day
                track.ori_order = i7
                push!(search, track)
            end
            
            # Move completed tracks to final array
            moved = [search[i2].day[end] <= day - 1 for i2 in 1:length(search)]
            
            for i2 in findall(moved)
                push!(tracks, search[i2])
            end
            
            deleteat!(search, findall(moved))
        end
        
        println("Completed day $i")
    end
    
    # Add remaining tracks
    append!(tracks, search)
    
    # Filter out short tracks if needed
    # tracks = filter(t -> length(t.day) >= cut_off, tracks)
    
    return tracks
end

# Export main functions
export hwtrack_with_boundary, bwconncomp, combine_boundary_objects

