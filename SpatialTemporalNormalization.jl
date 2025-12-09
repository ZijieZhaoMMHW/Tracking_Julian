"""
    SpatialTemporalNormalization.jl

Spatial and temporal normalization of tracked heat wave events.
This module normalizes heat wave tracks into a standardized circular grid
and temporal phases for composite analysis.

Based on MATLAB code by Zijie Zhao
Translated to Julia
"""

using Statistics
using Interpolations
using LinearAlgebra

"""
    geodist(lat1, lat2, lon1, lon2)

Calculate geodesic distance in kilometers using the haversine formula.

# Arguments
- `lat1`, `lat2`: Latitude coordinates (degrees)
- `lon1`, `lon2`: Longitude coordinates (degrees)

# Returns
- Distance in kilometers
"""
function geodist(lat1, lat2, lon1, lon2)
    R = 6371.0  # Earth radius in km
    
    # Convert to radians
    lat1_rad = deg2rad.(lat1)
    lat2_rad = deg2rad.(lat2)
    lon1_rad = deg2rad.(lon1)
    lon2_rad = deg2rad.(lon2)
    
    # Haversine formula
    dlat = lat2_rad .- lat1_rad
    dlon = lon2_rad .- lon1_rad
    
    a = sin.(dlat ./ 2).^2 .+ cos.(lat1_rad) .* cos.(lat2_rad) .* sin.(dlon ./ 2).^2
    c = 2 .* asin.(sqrt.(a))
    
    return R .* c
end


"""
    calculate_max_radius(tracks, lon_grid, lat_grid)

Calculate the maximum radius for each track.

# Arguments
- `tracks`: Vector of track structures with fields:
  - `day`: Vector of day indices
  - `xloc`: Vector of x-location vectors
  - `yloc`: Vector of y-location vectors
- `lon_grid`: Longitude grid (1D array)
- `lat_grid`: Latitude grid (1D array)

# Returns
- `radius_max`: Vector of maximum radii for each track (km)
"""
function calculate_max_radius(tracks, lon_grid, lat_grid)
    
    println("Calculating maximum radius for each track...")
    
    # Create 2D grid
    LON, LAT = meshgrid(lon_grid, lat_grid)
    
    n_tracks = length(tracks)
    radius_max = zeros(Float64, n_tracks)
    
    for i in 1:n_tracks
        if mod(i, 100) == 0
            println("Processing track $i / $n_tracks")
        end
        
        t_here = tracks[i].day
        xloc = tracks[i].xloc
        yloc = tracks[i].yloc
        
        radius_track = Float64[]
        
        for j in 1:length(t_here)
            # Get indices for this time step
            x_idx = xloc[j]
            y_idx = yloc[j]
            
            # Calculate centroid
            lon_points = [LON[Int(x), Int(y)] for (x, y) in zip(x_idx, y_idx)]
            lat_points = [LAT[Int(x), Int(y)] for (x, y) in zip(x_idx, y_idx)]
            
            center_lon = mean(lon_points)
            center_lat = mean(lat_points)
            
            # Calculate distances from centroid
            distances = geodist(lat_points, center_lat, lon_points, center_lon)
            
            # Maximum distance is the radius
            radius_j = maximum(distances)
            push!(radius_track, radius_j)
        end
        
        # Maximum radius over all time steps
        radius_max[i] = maximum(radius_track)
    end
    
    return radius_max
end


"""
    meshgrid(x, y)

Create 2D coordinate matrices from coordinate vectors.

# Arguments
- `x`: x-coordinate vector
- `y`: y-coordinate vector

# Returns
- `X, Y`: 2D coordinate matrices
"""
function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end


"""
    spatial_temporal_normalization(tracks, data_anomaly, lon_grid, lat_grid, 
                                   radius_max; res=50, n_phases=5)

Normalize heat wave tracks in space and time.

# Arguments
- `tracks`: Vector of track structures
- `data_anomaly`: 4D array of anomaly data (lon, lat, time, variables)
- `lon_grid`: Longitude grid (1D array)
- `lat_grid`: Latitude grid (1D array)  
- `radius_max`: Vector of maximum radii for each track (km)
- `res`: Resolution of normalized grid (default: 50)
- `n_phases`: Number of temporal phases (default: 5)

# Returns
- `normalized_data`: 5D array (res, res, n_phases, n_tracks, n_variables)
  containing spatially and temporally normalized data
"""
function spatial_temporal_normalization(tracks, data_anomaly, lon_grid, lat_grid, 
                                       radius_max; res=50, n_phases=5)
    
    n_tracks = length(tracks)
    n_vars = size(data_anomaly, 4)
    
    println("Starting spatial-temporal normalization...")
    println("Tracks: $n_tracks, Variables: $n_vars")
    println("Resolution: $res × $res, Phases: $n_phases")
    
    # Initialize output array
    normalized_data = fill(NaN, res, res, n_phases, n_tracks, n_vars)
    
    # Create 2D coordinate grids
    LON, LAT = meshgrid(lon_grid, lat_grid)
    
    # Create normalized circular grid (radius 0 to 1)
    radius_range = range(0, 1, length=res)
    angle_range = range(0, 2π, length=res)
    
    x_norm = zeros(res, res)
    y_norm = zeros(res, res)
    
    for r in 1:res
        for θ in 1:res
            x_norm[r, θ] = radius_range[r] * sin(angle_range[θ])
            y_norm[r, θ] = radius_range[r] * cos(angle_range[θ])
        end
    end
    
    # Temporal phases
    t_phases = range(0, 1, length=n_phases+1)[1:end-1]
    
    # Get time range
    t_min = minimum([minimum(track.day) for track in tracks])
    t_max = maximum([maximum(track.day) for track in tracks])
    
    # Process each variable
    for var_idx in 1:n_vars
        println("\nProcessing variable $var_idx / $n_vars")
        
        anom_data = data_anomaly[:, :, :, var_idx]
        
        # Process each track
        for i in 1:n_tracks
            if mod(i, 100) == 0
                println("  Track $i / $n_tracks")
            end
            
            t_here = tracks[i].day .- t_min .+ 1  # Adjust to 1-indexed
            xloc = tracks[i].xloc
            yloc = tracks[i].yloc
            
            n_timesteps = length(t_here)
            
            # Spatial normalization for each timestep
            spatial_data = fill(NaN, res, res, n_timesteps)
            
            for j in 1:n_timesteps
                # Get anomaly field at this time
                anom_field = anom_data[:, :, t_here[j]]
                
                # Get event pixels
                x_idx = Int.(xloc[j])
                y_idx = Int.(yloc[j])
                
                # Calculate centroid
                lon_points = [LON[x, y] for (x, y) in zip(x_idx, y_idx)]
                lat_points = [LAT[x, y] for (x, y) in zip(x_idx, y_idx)]
                
                center_lon = mean(lon_points)
                center_lat = mean(lat_points)
                
                # Transform to local coordinates (km from center)
                lon_local = zeros(size(LON))
                lat_local = zeros(size(LAT))
                
                for idx in CartesianIndices(LON)
                    # Calculate offset in km
                    lon_local[idx] = geodist(center_lat, center_lat, 
                                           LON[idx], center_lon)
                    lat_local[idx] = geodist(LAT[idx], center_lat, 
                                           center_lon, center_lon)
                    
                    # Apply sign based on direction
                    if LON[idx] < center_lon
                        lon_local[idx] = -lon_local[idx]
                    end
                    if LAT[idx] < center_lat
                        lat_local[idx] = -lat_local[idx]
                    end
                end
                
                # Scale by maximum radius to normalize to [0, 1]
                R_max = radius_max[i]
                lon_local_norm = lon_local ./ R_max
                lat_local_norm = lat_local ./ R_max
                
                # Create mask for event boundary
                event_mask = falses(size(anom_field))
                for (x, y) in zip(x_idx, y_idx)
                    event_mask[x, y] = true
                end
                
                # Find boundary points
                boundary_lon = lon_local_norm[event_mask]
                boundary_lat = lat_local_norm[event_mask]
                
                # Get data points within reasonable range
                valid_mask = (abs.(lon_local_norm) .< 2) .& (abs.(lat_local_norm) .< 2)
                
                lon_valid = lon_local_norm[valid_mask]
                lat_valid = lat_local_norm[valid_mask]
                anom_valid = anom_field[valid_mask]
                
                # Remove NaN values
                not_nan = .!isnan.(anom_valid)
                lon_valid = lon_valid[not_nan]
                lat_valid = lat_valid[not_nan]
                anom_valid = anom_valid[not_nan]
                
                if length(anom_valid) > 3
                    # Interpolate to normalized grid
                    try
                        itp = LinearInterpolation(
                            (lon_valid, lat_valid), 
                            anom_valid,
                            extrapolation_bc=Line()
                        )
                        
                        for r in 1:res
                            for θ in 1:res
                                try
                                    spatial_data[r, θ, j] = itp(x_norm[r, θ], y_norm[r, θ])
                                catch
                                    spatial_data[r, θ, j] = NaN
                                end
                            end
                        end
                    catch
                        # If interpolation fails, leave as NaN
                    end
                end
            end
            
            # Temporal normalization
            tb = range(0, 1, length=n_timesteps+1)[1:end-1]  # Original time
            ta = t_phases  # Target phases
            
            for r in 1:res
                for θ in 1:res
                    ts_original = spatial_data[r, θ, :]
                    
                    # Check if there's valid data
                    if !all(isnan.(ts_original))
                        # Interpolate to normalized phases
                        not_nan_idx = .!isnan.(ts_original)
                        
                        if sum(not_nan_idx) > 1
                            try
                                itp = LinearInterpolation(
                                    tb[not_nan_idx], 
                                    ts_original[not_nan_idx],
                                    extrapolation_bc=Line()
                                )
                                
                                for p in 1:n_phases
                                    normalized_data[r, θ, p, i, var_idx] = itp(ta[p])
                                end
                            catch
                                # Leave as NaN if interpolation fails
                            end
                        end
                    end
                end
            end
        end
    end
    
    println("\nNormalization complete!")
    return normalized_data
end


"""
    save_normalized_data(normalized_data, filename)

Save normalized data to file.

# Arguments
- `normalized_data`: 5D array of normalized data
- `filename`: Output filename (.jld2 or .nc)
"""
function save_normalized_data(normalized_data, filename)
    # Can be extended with JLD2 or NetCDF support
    println("Saving normalized data to $filename")
    println("Data size: $(size(normalized_data))")
    
    # Example with JLD2:
    # using JLD2
    # @save filename normalized_data
end


"""
    compute_composite(normalized_data; method="mean")

Compute composite across all tracks.

# Arguments
- `normalized_data`: 5D array (res, res, phases, tracks, variables)
- `method`: Aggregation method ("mean", "median", "std")

# Returns
- `composite`: 4D array (res, res, phases, variables)
"""
function compute_composite(normalized_data; method="mean")
    
    res, _, n_phases, n_tracks, n_vars = size(normalized_data)
    composite = zeros(res, res, n_phases, n_vars)
    
    for var_idx in 1:n_vars
        for p in 1:n_phases
            for r in 1:res
                for θ in 1:res
                    values = normalized_data[r, θ, p, :, var_idx]
                    values = values[.!isnan.(values)]
                    
                    if !isempty(values)
                        if method == "mean"
                            composite[r, θ, p, var_idx] = mean(values)
                        elseif method == "median"
                            composite[r, θ, p, var_idx] = median(values)
                        elseif method == "std"
                            composite[r, θ, p, var_idx] = std(values)
                        end
                    else
                        composite[r, θ, p, var_idx] = NaN
                    end
                end
            end
        end
    end
    
    return composite
end


# ============================================================================
# Example Usage
# ============================================================================

"""
# Example workflow

# 1. Load tracks from heat wave tracking
# tracks = load_tracks("tracks.jld2")

# 2. Load anomaly data (lon, lat, time, variables)
# data_anomaly = load_anomaly_data("anomalies.nc")
# lon_grid = 0.5:1.0:360
# lat_grid = -89.5:1.0:89.5

# 3. Calculate maximum radius for each track
radius_max = calculate_max_radius(tracks, lon_grid, lat_grid)

# 4. Perform spatial-temporal normalization
normalized_data = spatial_temporal_normalization(
    tracks, 
    data_anomaly,
    lon_grid,
    lat_grid,
    radius_max,
    res=50,        # 50×50 grid
    n_phases=5     # 5 temporal phases (0, 0.25, 0.5, 0.75, 1.0)
)

# 5. Compute composite
composite_mean = compute_composite(normalized_data, method="mean")
composite_std = compute_composite(normalized_data, method="std")

# 6. Save results
save_normalized_data(normalized_data, "normalized_hw_data.jld2")
"""

end  # module
