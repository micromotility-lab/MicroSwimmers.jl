abstract type Discretisation end

mutable struct NearestDiscretisation{T <: Number} <: Discretisation
    N::Int                          # number of force points
    Q::Int                          # number of quadrature points
    location::SVector{3,T}          # current location
    orientation::SMatrix{3,3,T}     # current orientation
    force_pts::Matrix{T}            # reference configuration force points
    velocity::Matrix{T}             # reference configuration velocities (at the force points)
    quad_pts::Matrix{T}             # reference configuration quadrature points
    nearest::Vector{Int}            # length Q
end

# convenience functions
traction_nodes(points::NearestDiscretisation)  = points.force_pts       # 3×N
n_traction_nodes(points) = size(traction_nodes(points), 2)
n_unknowns(points) = 3 * n_traction_nodes(points)


# Constructors
NearestDiscretisation(N::Int, Q::Int; location=SVector(0.,0.,0.), orientation=I3) = NearestDiscretisation(
    N,
    Q,
    SVector{3}(location),
    SMatrix{3,3}(orientation),
    zeros(3, N),
    zeros(3, N),
    zeros(3, Q),
    zeros(Int, Q)
)
    
NearestDiscretisation(force_pts::AbstractMatrix, quad_pts::AbstractMatrix; 
    location=SVector(0.,0.,0.), 
    orientation=I3
) = NearestDiscretisation(
    size(force_pts, 2),
    size(quad_pts, 2),
    location,
    orientation,
    force_pts,
    zeros(eltype(force_pts), size(force_pts)),
    quad_pts,
    nearest_neighbour(force_pts, quad_pts)
)

NearestDiscretisation(force_pts, quad_pts, nearest; 
    location=SVector(0.,0.,0.), 
    orientation=I3
) = NearestDiscretisation(
    size(force_pts, 2),
    size(quad_pts, 2),
    location,
    orientation,
    force_pts,
    zeros(eltype(force_pts), size(force_pts)),
    quad_pts,
    nearest
)

# spacing between force points, returns (min, med, max)
function hf(points::NearestDiscretisation)
    @unpack force_pts, N, = points
    dnn = [minimum(norm.(eachcol(force_pts .- force_pts[:, i]))[setdiff(1:N, i)]) for i in 1:N]
    @info "N=$N" findmin(dnn) median(dnn) maximum(dnn)
    (minimum(dnn), median(dnn), maximum(dnn))
end

# spacing between quadrature points, returns (min, med, max)
function hq(points::NearestDiscretisation)
    @unpack force_pts, quad_pts, N, Q, nearest = points
    dnn = Float64[]
    for i in 1:N
        patch_quad_pts = quad_pts[:, nearest .== i]
        Qp = size(patch_quad_pts, 2)
        append!(dnn, [minimum(norm.(eachcol(patch_quad_pts .- patch_quad_pts[:, j]))[setdiff(1:Qp, j)]) for j in 1:Qp])
    end
    @info "Q=$Q" minimum(dnn) median(dnn) findmax(dnn)
    (minimum(dnn), median(dnn), maximum(dnn))
end


# Calculate the nearest neighbour vector
function nearest_neighbour(force_pts::AbstractMatrix{T}, quad_pts::AbstractMatrix{T}) where {T <: Number}
    nearest = Int[]
    for x in eachcol(quad_pts)
        d = vec(sum((force_pts .- x).^2, dims=1))
       push!(nearest, argmin(d))
    end
    nearest
end 

# calculate the neighbour vector for a vector of parts
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
