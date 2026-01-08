module HWTrackMeshArray

using MeshArrays
using Images: label_components
using Statistics: mean

# =========================
# Data structures
# =========================

mutable struct HWDayObjects
    day::Int
    # 每个对象：用 (face, I, J) 三元组列表表示该对象包含的格点集合
    objects::Vector{Vector{NTuple{3,Int}}}
end

mutable struct Track
    day::Vector{Int}
    objects::Vector{Vector{NTuple{3,Int}}}  # 每天的对象格点集合
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

识别一天内所有热浪对象（connected components），其中跨 face 的连接通过 MeshArrays.exchange() + Union-Find 合并。

返回：
- objects :: Vector{Vector{(face,i,j)}}
"""
function components_mesh(mask::AbstractMeshArray{Bool}; minpix::Int=1)
    nf, _ = face_layout(mask)

    # 1) 每个 face 内部先做 8-connectivity 分割
    #    （face 内部还是规则数组，所以 label_components 很快）
    labels_face = Vector{Matrix{Int}}(undef, nf)
    nlab_face = zeros(Int, nf)

    for f in 1:nf
        M = mask.f[f]
        lab = label_components(M, trues(3,3))  # 8-connectivity on each face
        labels_face[f] = lab
        nlab_face[f] = maximum(lab)
    end

    # 2) 给每个 face 的 label 做全局编号：global_id = offset[f] + local_label
    offsets = cumsum(vcat(0, nlab_face[1:end-1]))
    total_labels = sum(nlab_face)

    # 如果全是空
    total_labels == 0 && return Vector{Vector{NTuple{3,Int}}}()

    # 3) 构建一个 Int MeshArray，存放每个 cell 的 global label（0 表示非热浪）
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

    # 4) 通过 exchange() 把边界邻接 face 的信息同步过来
    #    这里用 Int label 做 exchange：相邻 face 的边界会被复制到本 face 的边界对应位置
    #    这是 MeshArrays 的核心用途之一。:contentReference[oaicite:1]{index=1}
    glab_ex = MeshArrays.exchange(glab)

    # 5) Union-Find：如果某格点是热浪且它与“exchange 后获得的邻边 label”相邻（同一格点位置）
    #    就把二者 union（这一步相当于把跨 face 边界的连通块合并起来）
    dsu = DSU(total_labels)

    for f in 1:nf
        M = mask.f[f]
        G0 = glab.f[f]
        Gx = glab_ex.f[f]
        @inbounds for j in axes(M,2), i in axes(M,1)
            if M[i,j]
                a = G0[i,j]
                b = Gx[i,j]
                # 只有当 exchange 给到了“别的 label”时才 union
                if b != 0 && b != a
                    union!(dsu, a, b)
                end
            end
        end
    end

    # 6) 把每个点归类到 root label
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

    # 7) 过滤 minpix
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

两对象重叠比： |A∩B| / min(|A|,|B|)

注意：这里的“重叠”定义为完全相同的 (face,i,j) 点集合交集。
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

- hw: 每天一个 2D MeshArray(Bool) mask（LLC 13-face）
- nums: 最小像素数
- para_alpha: 继续/分裂判定的 overlap 阈值
- cut_off: 可选：短于该天数的事件过滤（默认不启用，和你原代码一致）
"""
function hwtrack_mesh(hw::Vector{<:AbstractMeshArray{Bool}};
                      nums::Int=10,
                      para_alpha::Float64=0.5,
                      cut_off::Int=5)

    nt = length(hw)

    # 1) 逐日识别对象
    HWs = Vector{HWDayObjects}(undef, nt)
    println("Identifying heat wave objects (MeshArray/LLC)...")
    for t in 1:nt
        objs = components_mesh(hw[t]; minpix=nums)
        HWs[t] = HWDayObjects(t, objs)
        println("  day $t : $(length(objs)) objects")
    end

    # 2) 追踪（基本保留你原来的 search/tracks 结构）
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

        # 记录今天每个对象被多少条 track 使用（用于 merge 检测）
        used_count = zeros(Int, length(objs_now))

        # 2.1 尝试把今天对象接到昨天的 tracks 上
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
                    # split：把多个对象合并到同一条 track 的当天（和你原代码一致的“合并存储”方式）
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

        # 2.2 merge：某个今天对象被多条 track 指向
        #     这里给一个“可用、但不做过度复杂谱系”的实现：
        #     - 找到所有指向同一对象的 tracks，把它们合并到第一条，其他条结束入库
        merged_targets = findall(used_count .> 1)
        if !isempty(merged_targets)
            for j in merged_targets
                # 找所有今天末尾对象 == objs_now[j] 的 track
                idx_tracks = Int[]
                for (it, tr) in pairs(search)
                    if tr.day[end] == day
                        # 这里用集合相等判断（按点集）
                        if Set(tr.objects[end]) == Set(objs_now[j])
                            push!(idx_tracks, it)
                        end
                    end
                end
                if length(idx_tracks) > 1
                    # 合并到第一条
                    keeper = idx_tracks[1]
                    for k in idx_tracks[2:end]
                        push!(tracks, search[k])  # 被合并的 track 入库（你也可以改成记录 parent/child）
                    end
                    # 从 search 中删除其他条（注意倒序删）
                    for k in sort(idx_tracks[2:end]; rev=true)
                        deleteat!(search, k)
                    end
                end
            end
        end

        # 2.3 新事件：没被用过的对象，开新 track
        new_idx = findall(used_count .== 0)
        for j in new_idx
            push!(search, Track([day], [objs_now[j]], day, j, Int[], Int[]))
        end

        # 2.4 结束：昨天没续上的 tracks 入库
        moved = [tr.day[end] <= day - 1 for tr in search]
        for i in reverse(eachindex(search))
            if moved[i]
                push!(tracks, search[i])
                deleteat!(search, i)
            end
        end
    end

    # 3) 收尾：剩余 search 入库
    append!(tracks, search)

    # 4) 可选过滤短事件（默认不启用；如要启用就取消下一行注释）
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
