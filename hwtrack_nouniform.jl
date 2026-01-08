using Images
using ImageMorphology

"""
    hwtrack_nouniform(hw, lon_used, lat_used, nums)

Heat wave tracking algorithm - non-uniform grid version
Originally written by Sun Di
Modified and rewritten by Zijie Zhao

# Arguments
- `hw`: 3D array of heat wave masks (nx, ny, nt)
- `lon_used`: longitude array
- `lat_used`: latitude array  
- `nums`: minimum number of pixels for a valid heat wave

# Returns
- `tracks`: Vector of tracked heat wave events
"""
function hwtrack_nouniform(hw, lon_used, lat_used, nums)
    # Define structure types
    mutable struct HW
        day::Union{Int, Nothing}
        xloc::Vector{Vector{Float32}}
        yloc::Vector{Vector{Float32}}
    end
    
    mutable struct Track
        day::Vector{Int}
        xloc::Vector{Vector{Float32}}
        yloc::Vector{Vector{Float32}}
        ori_day::Int
        ori_order::Int
        split_num::Vector{Int}
        split_day::Vector{Int}
    end
    
    # Initialize HWs array
    HWs = HW[]
    
    # Identify heat wave events for each day
    println("Identifying heat wave events...")
    for i in 1:size(hw, 3)
        mask = hw[:, :, i]
        components = label_components(mask .> 0, trues(3, 3))  # 8-connectivity
        
        xloc_temp = Vector{Float32}[]
        yloc_temp = Vector{Float32}[]
        
        num_objects = maximum(components)
        for j in 1:num_objects
            indices = findall(components .== j)
            x = [idx[1] for idx in indices]
            y = [idx[2] for idx in indices]
            
            if length(x) >= nums
                push!(xloc_temp, Float32.(x))
                push!(yloc_temp, Float32.(y))
            end
        end
        
        push!(HWs, HW(i, xloc_temp, yloc_temp))
        println("Processing day $i")
    end
    
    # Handle empty HWs
    for i in 1:length(HWs)
        if isempty(HWs[i].xloc)
            HWs[i].day = i
        end
    end
    
    # Set parameters
    para_alpha = 0.5
    cut_off = 5
    
    pi180 = π / 180
    earth_radius = 6378.137
    
    # Initialize tracking arrays
    search = Track[]
    tracks = Track[]
    
    nx, ny = size(hw, 1), size(hw, 2)
    
    # Start tracking
    println("\nStarting heat wave tracking...")
    for i in 1:length(HWs)
        day = HWs[i].day
        hw_xloc = HWs[i].xloc
        hw_yloc = HWs[i].yloc
        
        # First day, initialize tracks
        if i == 1
            for i2 in 1:length(hw_xloc)
                track = Track(
                    [day],
                    [hw_xloc[i2]],
                    [hw_yloc[i2]],
                    day,
                    i2,
                    Int[],
                    Int[]
                )
                push!(search, track)
            end
        else
            # Not first day, process overlap and tracking
            loc_now = []
            for i2 in 1:length(hw_xloc)
                push!(loc_now, (xloc=hw_xloc[i2], yloc=hw_yloc[i2], 
                               ori_day=i, ori_order=i2))
            end
            
            count = zeros(Int, length(hw_xloc))
            
            # Loop through existing tracks
            for i2 in 1:length(search)
                if (day - search[i2].day[end]) == 1
                    loc_old_x = search[i2].xloc[end]
                    loc_old_y = search[i2].yloc[end]
                    
                    # Create mask
                    judge1 = zeros(Bool, nx, ny)
                    for k in 1:length(loc_old_x)
                        judge1[Int(loc_old_x[k]), Int(loc_old_y[k])] = true
                    end
                    
                    overlap = zeros(Float64, length(loc_now))
                    
                    # Calculate overlap with all current HWs
                    for i3 in 1:length(loc_now)
                        judge2 = zeros(Bool, nx, ny)
                        loc_now_x = loc_now[i3].xloc
                        loc_now_y = loc_now[i3].yloc
                        
                        for k2 in 1:length(loc_now_x)
                            judge2[Int(loc_now_x[k2]), Int(loc_now_y[k2])] = true
                        end
                        
                        # Calculate overlap ratio
                        overlap_count = sum(judge1 .& judge2)
                        overlap[i3] = overlap_count / min(sum(judge2), sum(judge1))
                    end
                    
                    idx = findall(overlap .>= para_alpha)
                    
                    # Handle splitting cases
                    if !isempty(idx)
                        if length(idx) > 1
                            # Splitting
                            push!(search[i2].day, day)
                            push!(search[i2].xloc, hw_xloc[idx[1]])
                            push!(search[i2].yloc, hw_yloc[idx[1]])
                            
                            println("Split!! Day $i")
                            
                            # Merge all split parts
                            for i4 in 2:length(idx)
                                search[i2].xloc[end] = vcat(search[i2].xloc[end], 
                                                            hw_xloc[idx[i4]])
                                search[i2].yloc[end] = vcat(search[i2].yloc[end], 
                                                            hw_yloc[idx[i4]])
                            end
                            
                            push!(search[i2].split_num, length(idx))
                            push!(search[i2].split_day, day)
                        else
                            # Normal continuation
                            push!(search[i2].day, day)
                            push!(search[i2].xloc, hw_xloc[idx[1]])
                            push!(search[i2].yloc, hw_yloc[idx[1]])
                        end
                    end
                    
                    # Record usage count (for merger detection)
                    count[idx] .+= 1
                end
            end
            
            # Handle merging cases
            if any(count .> 1)
                idx_now = findall(count .> 1)
                old = zeros(Int, length(idx_now), length(search))
                
                for m in 1:length(idx_now)
                    loc_now_xloc = loc_now[idx_now[m]].xloc
                    loc_now_yloc = loc_now[idx_now[m]].yloc
                    
                    c = 0
                    for n in 1:length(search)
                        search_xloc = search[n].xloc[end]
                        search_yloc = search[n].yloc[end]
                        
                        l = zeros(Bool, nx, ny)
                        for j5 in 1:length(loc_now_xloc)
                            l[Int(loc_now_xloc[j5]), Int(loc_now_yloc[j5])] = true
                        end
                        
                        ll = zeros(Bool, nx, ny)
                        for j5 in 1:length(search_xloc)
                            ll[Int(search_xloc[j5]), Int(search_yloc[j5])] = true
                        end
                        
                        if day == search[n].day[end] && 
                           ((length(loc_now_xloc) == length(search_xloc) && 
                             loc_now_xloc == search_xloc) || 
                            (sum(l .& ll) == sum(l)))
                            c += 1
                            old[m, c] = n
                        end
                    end
                end
                
                # Detailed logic for handling mergers
                new_fusion = []
                new_fusion_loc = Int[]
                
                for i5 in 1:length(idx_now)
                    num_old_tracks = sum(old[i5, :] .!= 0)
                    
                    if num_old_tracks > 1
                        println("Merge!! Day $i")
                        
                        # Merge multiple tracks
                        # Simplified handling: merge all related track info into the first track
                        first_idx = old[i5, findfirst(old[i5, :] .!= 0)]
                        
                        for j in 2:num_old_tracks
                            track_idx = old[i5, j]
                            if track_idx != 0 && track_idx != first_idx
                                # Merging logic (simplified version)
                                # More complex handling may be needed in practice
                            end
                        end
                    end
                end
            end
            
            # New tracks
            new_hw_idx = findall(count .== 0)
            if !isempty(new_hw_idx)
                for i7 in new_hw_idx
                    track = Track(
                        [day],
                        [hw_xloc[i7]],
                        [hw_yloc[i7]],
                        day,
                        loc_now[i7].ori_order,
                        Int[],
                        Int[]
                    )
                    push!(search, track)
                end
            end
            
            # Move dissipated tracks
            moved = [search[i2].day[end] <= day - 1 for i2 in 1:length(search)]
            for i2 in 1:length(search)
                if moved[i2]
                    push!(tracks, search[i2])
                end
            end
            deleteat!(search, moved)
        end
        
        println("Processing day $i")
    end
    
    # Add remaining tracks in the search array at the end of the time-series
    for i in 1:length(search)
        push!(tracks, search[i])
    end
    
    # Optional: remove tracks shorter than cut_off days
    # filter!(t -> length(t.day) >= cut_off, tracks)
    
    return tracks
end

# Helper function: save results
function save_tracks(tracks, filename)
    # Can add file saving functionality here
    println("Total tracked heat wave events: $(length(tracks))")
    
    # Statistics
    durations = [length(t.day) for t in tracks]
    println("Average duration: $(mean(durations)) days")
    println("Maximum duration: $(maximum(durations)) days")
    println("Minimum duration: $(minimum(durations)) days")
end

# Usage example
"""
# Load data
using NCDatasets
ds = NCDataset("your_data.nc")
hw = ds["hw"][:]
lon = ds["lon"][:]
lat = ds["lat"][:]
close(ds)

# Run tracking
tracks = hwtrack_nouniform(hw, lon, lat, 10)

# Save results
save_tracks(tracks, "tracks_output.jl")
"""
