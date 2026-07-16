abstract type Discretisation end

struct NystromDiscretisation{T <: Number} <: Discretisation
    force_pts::Vector{SVector{3,T}}
    velocity::Vector{SVector{3,T}}
end

NystromDiscretisation(N::Int) = NystromDiscretisation(
    Vector{SVector{3,Float64}}(undef, N),
    Vector{SVector{3,Float64}}(undef, N)
)

nf(disc::NystromDiscretisation) = length(disc.force_pts)
nq(disc::NystromDiscretisation) = length(disc.force_pts)

# mutable struct NearestDiscretisation{T <: Number} <: Discretisation
#     N::Int                          # number of force points
#     Q::Int                          # number of quadrature points
#     location::SVector{3,T}          # current location
#     orientation::SMatrix{3,3,T}     # current orientation
#     force_pts::Matrix{T}            # reference configuration force points
#     velocity::Matrix{T}             # reference configuration velocities (at the force points)
#     quad_pts::Matrix{T}             # reference configuration quadrature points
#     nearest::Vector{Int}            # length Q
# end


# n_traction_nodes(points) = size(traction_nodes(points), 2)
# n_unknowns(points) = 3 * n_traction_nodes(points)


# NearestDiscretisation(N::Int, Q::Int; location::SVector{3, T}=SVector(0.,0.,0.), orientation=I3) where {T <: Number} = NearestDiscretisation(
#     N,
#     Q,
#     location,
#     orientation,
#     zeros(T, 3, N),
#     zeros(T, 3, N),
#     zeros(T, 3, Q),
#     zeros(Int, Q)
# )
    
# function NearestDiscretisation(::Type{T}, 
#     N::Int, Q::Int; 
#     location=SVector{3,T}(0, 0, 0), orientation=SMatrix{3,3,T}(I)
#     ) where {T <: Number}
#     NearestDiscretisation(
#         N,
#         Q,
#         location,
#         orientation,
#         zeros(T, 3, N),
#         zeros(T, 3, N),
#         zeros(T, 3, Q),
#         zeros(Int, Q)
#     )
# end
    
# traction_nodes(points::NearestDiscretisation) = points.force_pts  # 3×Nf

# NearestDiscretisation(force_pts::AbstractMatrix{T}, quad_pts::AbstractMatrix{T}; location=SVector(0.,0.,0.), orientation=I3) where {T <: Number} = NearestDiscretisation(
#     size(force_pts, 2),
#     size(quad_pts, 2),
#     location,
#     orientation,
#     force_pts,
#     zeros(eltype(force_pts), size(force_pts)),
#     quad_pts,
#     nearest_neighbour(force_pts, quad_pts)
# )


# NearestDiscretisation(force_pts, quad_pts, nearest; location=SVector(0.,0.,0.), orientation=I3) = NearestDiscretisation(
#     size(force_pts, 2),
#     size(quad_pts, 2),
#     location,
#     orientation,
#     force_pts,
#     zeros(eltype(force_pts), size(force_pts)),
#     quad_pts,
#     nearest
# )
        
# function spacing(points::NearestDiscretisation)
#     @unpack force_pts, quad_pts, N, Q, nearest = points
#     dnn = [minimum(norm.(eachcol(force_pts .- force_pts[:, i]))[setdiff(1:N, i)]) for i in 1:N]
#     hf = (minimum(dnn), median(dnn), maximum(dnn))
#     @info "" findmin(dnn)
#     @info "hf (min, median, max)" hf
#     # hf = median(dnn)

#     dnn = Float64[]
#     for i in 1:N
#         patch_quad_pts = quad_pts[:, nearest .== i]
#         Qp = size(patch_quad_pts, 2)
#         append!(dnn, [minimum(norm.(eachcol(patch_quad_pts .- patch_quad_pts[:, j]))[setdiff(1:Qp, j)]) for j in 1:Qp])
#     end
#     @info "" findmax(dnn)
#     hq = (minimum(dnn), median(dnn), maximum(dnn))
#     @info "hq (min, median, max)" hq
#     hf, hq
# end

# mutable struct NearestDiscretisation{T <: Number} <: Discretisation
#     force_pts::Vector{SVector{3,T}}            # reference configuration force points
#     velocity::Vector{SVector{3,T}}             # reference configuration velocities (at the force points)
#     quad_pts::Vector{SVector{3,T}}             # reference configuration quadrature points
#     nearest::Vector{Int}            # length Q
# end

mutable struct NearestDiscretisation{T <: Number} <: Discretisation
    force_pts::Vector{SVector{3,T}}
    velocity::Vector{SVector{3,T}}
    quad_pts::Vector{SVector{3,T}}
    nearest::Vector{Int}
    force_part_ranges::Vector{UnitRange{Int}}
    quad_part_ranges::Vector{UnitRange{Int}}
end

function ranges_from_sizes(sizes)
    ranges = Vector{UnitRange{Int}}(undef, length(sizes))
    start = 1
    for (i, s) in enumerate(sizes)
        ranges[i] = start:(start + s - 1)
        start += s
    end
    ranges
end
# one-liner equivalent:
# ranges_from_sizes(sizes) = (stop = cumsum(sizes); UnitRange.(stop .- sizes .+ 1, stop))

# size-driven constructor knows the partition, so it fills the ranges
function NearestDiscretisation(nf_sizes::AbstractVector{<:Integer},
                               nq_sizes::AbstractVector{<:Integer})
    N = sum(nf_sizes); Q = sum(nq_sizes)
    NearestDiscretisation(
        Vector{SVector{3,Float64}}(undef, N),
        Vector{SVector{3,Float64}}(undef, N),
        Vector{SVector{3,Float64}}(undef, Q),
        zeros(Int, Q),
        ranges_from_sizes(nf_sizes),
        ranges_from_sizes(nq_sizes),
    )
end

# keep (N, Q) as a single-part fallback
NearestDiscretisation(N::Int, Q::Int) = NearestDiscretisation([N], [Q])

NearestDiscretisation() = NearestDiscretisation(
    SVector{3,Float64}[],
    SVector{3,Float64}[],
    SVector{3,Float64}[],
    Int[],
    UnitRange{Int}[],
    UnitRange{Int}[],
)

# subview(d::NearestDiscretisation, f_rng, q_rng) = NearestDiscretisation(
#     view(d.force_pts, f_rng), 
#     view(d.velocity, f_rng), 
#     view(d.quad_pts, q_rng),
#     view(d.nearest, q_rng)
# )

# NearestDiscretisation() = NearestDiscretisation(
#     SVector{3,Float64}[],
#     SVector{3,Float64}[],
#     SVector{3,Float64}[],
#     Int[]
# )

# NearestDiscretisation(N::Int, Q::Int) = NearestDiscretisation(
#     Vector{SVector{3,Float64}}(undef, N),
#     Vector{SVector{3,Float64}}(undef, N),
#     Vector{SVector{3,Float64}}(undef, Q),
#     zeros(Int, Q)
# )

nf(disc::NearestDiscretisation) = length(disc.force_pts)
nq(disc::NearestDiscretisation) = length(disc.quad_pts)
    

# spacing between force points
function hf(disc::NearestDiscretisation)
    fps = disc.force_pts
    N   = nf(disc)
    dnn = [minimum(norm(fps[i] - fps[j]) for j in 1:N if j != i) for i in 1:N]
    @info "N=$N" findmin(dnn) median(dnn) maximum(dnn)
    (minimum(dnn), median(dnn), maximum(dnn))
end

# spacing between quadrature points

function hq(disc::NearestDiscretisation)
    @unpack quad_pts, nearest = disc
    N   = nf(disc)
    dnn = Float64[]
    for i in 1:N
        patch = [quad_pts[j] for j in eachindex(quad_pts) if nearest[j] == i]
        Qp = length(patch)
        Qp > 1 && append!(dnn, [minimum(norm(patch[j] - patch[k]) for k in 1:Qp if k != j) for j in 1:Qp])
    end
    @info "Q=$(nq(disc))" minimum(dnn) median(dnn) findmax(dnn)
    (minimum(dnn), median(dnn), maximum(dnn))
end

# # Forward unknown properties to `points`
# @inline function Base.getproperty(f::VanedFlagellumNearestDiscretisation, name::Symbol)
#     name in (:points, :N_f, :Q_f, :N_v, :N_start, :N_height, :Q_v, :Q_start, :Q_height) ? getfield(f, name) : getproperty(f.points, name)  
# end

# @inline function Base.setproperty!(f::VanedFlagellumNearestDiscretisation, name::Symbol, value)
#     name in (:points, :N_f, :Q_f, :N_v, :N_start, :N_height, :Q_v, :Q_start, :Q_height) ? setfield!(f, name, value) : setproperty!(f.points, name, value)
# end

# function VanedFlagellumNearestDiscretisation(
#     N_f::Int,
#     Q_f::Int, 
#     N_v::Int, 
#     N_start::Int, 
#     N_height::Int; 
#     location=SVector(0.,0.,0.), 
#     orientation=I3
# ) 
#     Q_v = floor(Int, (N_v / N_f) * Q_f) 
#     Q_start = ceil(Int, ((N_start-1) / (N_f-1)) * (Q_f-1))
#     Q_height = floor(Int, (N_height / N_f) * Q_f)

#     points = NearestDiscretisation(
#         N_f + N_height*N_v, Q_f + Q_height*Q_v;
#         location=location, orientation=orientation
#     )

#     VanedFlagellumNearestDiscretisation(points, N_f, Q_f, N_v, N_start, N_height, Q_v, Q_start, Q_height)
# end


struct TubeFlagellumNearestDiscretisation{T <: Number} <: Discretisation
    points::NearestDiscretisation

    N_cs::Int  # Number of cross section force points on the tube
    Q_cs::Int  # Number of cross section quadrature points on the tube
    radius::T  # cross-sectional radius of the tube
end

# Forward unknown properties to `points`
@inline function Base.getproperty(f::TubeFlagellumNearestDiscretisation, name::Symbol)
    name in (:points, :N_cs, :Q_cs, :radius) ? getfield(f, name) : getproperty(f.points, name)  
end

@inline function Base.setproperty!(f::TubeFlagellumNearestDiscretisation, name::Symbol, value)
    name in (:points, :N_cs, :Q_cs, :radius) ? setfield!(f, name, value) : setproperty!(f.points, name, value)
end


function TubeFlagellumNearestDiscretisation(
    N::Int,
    N_cs::Int, 
    Q::Int, 
    Q_cs::Int; 
    location=SVector(0.,0.,0.), 
    orientation=I3,
    radius=0.01
)

    points = NearestDiscretisation(N*N_cs, Q*Q_cs, location=location, orientation=orientation)
    TubeFlagellumNearestDiscretisation(points, N_cs, Q_cs, radius)
end


struct LineTubeFlagellumNearestDiscretisation{T <: Number} <: Discretisation
    points::NearestDiscretisation

    Q_cs::Int  # Number of cross section quadrature points on the tube
    radius::T  # cross-sectional radius of the tube
end

# Forward unknown properties to `points`
@inline function Base.getproperty(f::LineTubeFlagellumNearestDiscretisation, name::Symbol)
    name in (:points, :Q_cs, :radius) ? getfield(f, name) : getproperty(f.points, name)  
end

@inline function Base.setproperty!(f::LineTubeFlagellumNearestDiscretisation, name::Symbol, value)
    name in (:points, :Q_cs, :radius) ? setfield!(f, name, value) : setproperty!(f.points, name, value)
end


function LineTubeFlagellumNearestDiscretisation(
    N::Int, 
    Q::Int, 
    Q_cs::Int; 
    location=SVector(0.,0.,0.), 
    orientation=I3,
    radius=0.01
)
    points = NearestDiscretisation(N, Q*Q_cs, location=location, orientation=orientation)
    LineTubeFlagellumNearestDiscretisation(points, Q_cs, radius)
end




function nearest_neighbour(force_pts::AbstractMatrix{T}, quad_pts::AbstractMatrix{T}) where {T <: Number}
    nearest = Int[]
    for x in eachcol(quad_pts)
        d = vec(sum((force_pts .- x).^2, dims=1))
       push!(nearest, argmin(d))
    end
    nearest
end 

function nearest_neighbour(force_pts::Vector{SVector{3,T}}, quad_pts::Vector{SVector{3,T}}) where {T <: Number}
    nearest = Int[]
    for x in quad_pts
        d = norm.(force_pts .- Ref(x))
       push!(nearest, argmin(d))
    end
    nearest
end 

function nearest_neighbour(force_pts::AbstractVector{<:AbstractMatrix{T}}, quad_pts::AbstractVector{<:AbstractMatrix{T}}) where {T <: Number}
    nearest = Int[]
    idx = 0
    for (f_pts, q_pts) in zip(force_pts, quad_pts)
        nn = nearest_neighbour(f_pts, q_pts)
        append!(nearest, idx .+ nn)
        idx += size(f_pts, 2)
    end
    nearest
end



function nearest_neighbour!(points::Discretisation)
    points.nearest .= nearest_neighbour(points.force_pts, points.quad_pts)
end
