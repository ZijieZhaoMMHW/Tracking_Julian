module HWTrackMeshArray

using MeshArrays
using Images: label_components
using Statistics: mean

# =========================
# Data structures
# =========================

mutable struct HWDayObjects
    day::Int
    # Each object: represented by a list of (face, I, J) tuples containing the grid points
    objects::Vector{Vector{NTuple{3,Int}}}
end

mutable struct Track
    day::Vector{Int}
    objects::Vector{Vector{NTuple{3,Int}}}  # Grid point set for each day
    ori_day::Int
    ori_order::Int
    split_num::Vector{Int}
    split_day::Vector{Int}
end

# =========================
# Union-Find (Disjoint Set)
# =========================

mutable struct DSU
    parent::Vector{Int}
    rank::Vector{Int}
end

DSU(n::Int) = DSU(collect(1:n), zeros(Int, n))

function find!(d::DSU, x::Int)
    while d.parent[x] != x
        d.parent[x] = d.parent[d.parent[x]]
        x = d.parent[x]
    end
    return x
end

function union!(d::DSU, a::Int, b::Int)
    ra = find!(d, a)
    rb = find!(d, b)
    ra == rb && return ra
    if d.rank[ra] < d.rank[rb]
        d.parent[ra] = rb
        return rb
    elseif d.rank[ra] > d.rank[rb]
        d.parent[rb] = ra
        return ra
    else
        d.parent[rb] = ra
        d.rank[ra] += 1
        return ra
    end
end

# =========================
# Helpers: work with MeshArray faces
# =========================

"""
Return number of faces and per-face sizes for a 2D MeshArray.
"""
function face_layout(A::AbstractMeshArray)
    nf = length(A.f)                       # A.f is the array-of-faces in MeshArrays.gcmarray
    sizes = [size(A.f[k]) for k in 1:nf]
    return nf, sizes
end

"""
Make an Int MeshArray same grid as template A, initialized to zero.
"""
function zeros_like_int(A::AbstractMeshArray)
    B = similar(A, Int)
    for k in eachindex(B.f)
        fill!(B.f[k], 0)
    end
    return B
end

"""
Iterate all active points (mask==true) and apply a function to (face,i,j).
"""
function each_active(mask::AbstractMeshArray{Bool}, fn::Function)
    for f in eachindex(mask.f)
        M = mask.f[f]
        @inbounds for j in axes(M,2), i in axes(M,1)
            if M[i,j]
                fn(f,i,j)
            end
        end
    end
end

# =========================
# Core: connected components on LLC faces with face-to-face merging
# =========================

"""
    components_mesh(mask::MeshArray{Bool}; minpix=1)

Identify all heat wave objects (connected components) within a day, where cross-face connections 
are merged via MeshArrays.exchange() + Union-Find.

Returns:
- objects :: Vector{Vector{(face,i,j)}}
"""
function components_mesh(mask::AbstractMeshArray{Bool}; minpix::Int=1)
    nf, _ = face_layout(mask)

    # 1) First perform 8-connectivity segmentation within each face
    #    (each face is still a regular array, so label_components is fast)
    labels_face = Vector{Matrix{Int}}(undef, nf)
    nlab_face = zeros(Int, nf)

    for f in 1:nf
        M = mask.f[f]
        lab = label_components(M, trues(3,3))  # 8-connectivity on each face
        labels_face[f] = lab
        nlab_face[f] = maximum(lab)
    end

    # 2) Assign global numbering for each face's label: global_id = offset[f] + local_label
    offsets = cumsum(vcat(0, nlab_face[1:end-1]))
    total_labels = sum(nlab_face)

    # If everything is empty
    total_labels == 0 && return Vector{Vector{NTuple{3,Int}}}()

    # 3) Build an Int MeshArray storing each cell's global label (0 means non-heatwave)
    glab = zeros_like_int(mask)
    for f in 1:nf
        M = mask.f[f]
        L = labels_face[f]
        off = offsets[f]
        G = glab.f[f]
        @inbounds for j in axes(M,2), i in axes(M,1)
            if M[i,j]
                G[i,j] = off + L[i,j]   # L[i,j] >= 1
            end
        end
    end

    # 4) Use exchange() to synchronize boundary information from adjacent faces
    #    Here we exchange Int labels: boundaries of adjacent faces are copied to 
    #    corresponding boundary positions of this face
    #    This is one of the core uses of MeshArrays.
    glab_ex = MeshArrays.exchange(glab)

    # 5) Union-Find: if a grid point is a heatwave and it's adjacent to a label obtained
    #    from exchange at the same position, then union them (this step merges connected components across face boundaries)
    dsu = DSU(total_labels)

    for f in 1:nf
        M = mask.f[f]
        G0 = glab.f[f]
        Gx = glab_ex.f[f]
        @inbounds for j in axes(M,2), i in axes(M,1)
            if M[i,j]
                a = G0[i,j]
                b = Gx[i,j]
                # Only union when exchange provides a different label
                if b != 0 && b != a
                    union!(dsu, a, b)
                end
            end
        end
    end

    # 6) Classify each point to its root label
    buckets = Dict{Int, Vector{NTuple{3,Int}}}()
    for f in 1:nf
        M = mask.f[f]
        G0 = glab.f[f]
        @inbounds for j in axes(M,2), i in axes(M,1)
            if M[i,j]
                a = G0[i,j]
                r = find!(dsu, a)
                if !haskey(buckets, r)
                    buckets[r] = NTuple{3,Int}[]
                end
                push!(buckets[r], (f,i,j))
            end
        end
    end

    # 7) Filter by minpix
    objects = Vector{Vector{NTuple{3,Int}}}()
    for (_k, pts) in buckets
        if length(pts) >= minpix
            push!(objects, pts)
        end
    end

    return objects
end

"""
    overlap_ratio(objA, objB)

Overlap ratio of two objects: |A∩B| / min(|A|,|B|)

Note: overlap here is defined as the intersection of identical (face,i,j) point sets.
"""
function overlap_ratio(objA::Vector{NTuple{3,Int}}, objB::Vector{NTuple{3,Int}})
    setA = Set(objA)
    setB = Set(objB)
    inter = length(intersect(setA, setB))
    return inter / min(length(setA), length(setB))
end

# =========================
# Main tracking
# =========================

"""
    hwtrack_mesh(hw::Vector{MeshArray{Bool}}; nums=10, para_alpha=0.5, cut_off=5)

- hw: A 2D MeshArray(Bool) mask for each day (LLC 13-face)
- nums: Minimum number of pixels
- para_alpha: Overlap threshold for continuation/split detection
- cut_off: Optional: filter events shorter than this many days (disabled by default, consistent with original code)
"""
function hwtrack_mesh(hw::Vector{<:AbstractMeshArray{Bool}};
                      nums::Int=10,
                      para_alpha::Float64=0.5,
                      cut_off::Int=5)

    nt = length(hw)

    # 1) Identify objects day by day
    HWs = Vector{HWDayObjects}(undef, nt)
    println("Identifying heat wave objects (MeshArray/LLC)...")
    for t in 1:nt
        objs = components_mesh(hw[t]; minpix=nums)
        HWs[t] = HWDayObjects(t, objs)
        println("  day $t : $(length(objs)) objects")
    end

    # 2) Tracking (basically keeping your original search/tracks structure)
    search = Track[]
    tracks = Track[]

    for t in 1:nt
        day = HWs[t].day
        objs_now = HWs[t].objects

        if t == 1
            for k in 1:length(objs_now)
                push!(search, Track([day], [objs_now[k]], day, k, Int[], Int[]))
            end
            continue
        end

        # Record how many tracks use each object today (for merge detection)
        used_count = zeros(Int, length(objs_now))

        # 2.1 Try to attach today's objects to yesterday's tracks
        for tr in search
            if day - tr.day[end] != 1
                continue
            end

            obj_old = tr.objects[end]
            overlaps = zeros(Float64, length(objs_now))
            for j in 1:length(objs_now)
                overlaps[j] = overlap_ratio(obj_old, objs_now[j])
            end

            idx = findall(overlaps .>= para_alpha)

            if !isempty(idx)
                if length(idx) > 1
                    # split: merge multiple objects into the same track for today (consistent with original merged storage approach)
                    push!(tr.day, day)
                    merged = Vector{NTuple{3,Int}}()
                    for j in idx
                        append!(merged, objs_now[j])
                        used_count[j] += 1
                    end
                    push!(tr.objects, merged)
                    push!(tr.split_num, length(idx))
                    push!(tr.split_day, day)
                else
                    # normal continuation
                    j = idx[1]
                    push!(tr.day, day)
                    push!(tr.objects, objs_now[j])
                    used_count[j] += 1
                end
            end
        end

        # 2.2 merge: a today's object is pointed to by multiple tracks
        #     Here we provide a usable but not overly complex genealogy implementation:
        #     - Find all tracks pointing to the same object, merge them into the first one, and archive the others
        merged_targets = findall(used_count .> 1)
        if !isempty(merged_targets)
            for j in merged_targets
                # Find all tracks whose last object today == objs_now[j]
                idx_tracks = Int[]
                for (it, tr) in pairs(search)
                    if tr.day[end] == day
                        # Use set equality to judge (by point set)
                        if Set(tr.objects[end]) == Set(objs_now[j])
                            push!(idx_tracks, it)
                        end
                    end
                end
                if length(idx_tracks) > 1
                    # Merge into the first one
                    keeper = idx_tracks[1]
                    for k in idx_tracks[2:end]
                        push!(tracks, search[k])  # Archive merged track (you can also change to record parent/child)
                    end
                    # Delete other tracks from search (note reverse order deletion)
                    for k in sort(idx_tracks[2:end]; rev=true)
                        deleteat!(search, k)
                    end
                end
            end
        end

        # 2.3 New events: unused objects, create new track
        new_idx = findall(used_count .== 0)
        for j in new_idx
            push!(search, Track([day], [objs_now[j]], day, j, Int[], Int[]))
        end

        # 2.4 End: archive tracks that didn't continue from yesterday
        moved = [tr.day[end] <= day - 1 for tr in search]
        for i in reverse(eachindex(search))
            if moved[i]
                push!(tracks, search[i])
                deleteat!(search, i)
            end
        end
    end

    # 3) Wrap up: archive remaining search
    append!(tracks, search)

    # 4) Optional filter for short events (disabled by default; uncomment next line to enable)
    # filter!(tr -> length(tr.day) >= cut_off, tracks)

    return tracks
end

function save_tracks_summary(tracks)
    println("Total tracked heat wave events: $(length(tracks))")
    durations = [length(t.day) for t in tracks]
    if !isempty(durations)
        println("Average duration: $(mean(durations)) days")
        println("Maximum duration: $(maximum(durations)) days")
        println("Minimum duration: $(minimum(durations)) days")
    end
end

end # module
