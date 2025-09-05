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

NearestDiscretisation(force_pts, quad_pts; location=SVector(0.,0.,0.), orientation=I3) = NearestDiscretisation(
    size(force_pts, 2),
    size(quad_pts, 2),
    location,
    orientation,
    force_pts,
    zeros(eltype(force_pts), size(force_pts)),
    quad_pts,
    nearest_neighbour(force_pts, quad_pts)
)

NearestDiscretisation(force_pts, quad_pts, nearest; location=SVector(0.,0.,0.), orientation=I3) = NearestDiscretisation(
    size(force_pts, 2),
    size(quad_pts, 2),
    location,
    orientation,
    force_pts,
    zeros(eltype(force_pts), size(force_pts)),
    quad_pts,
    nearest
)

struct TubeFlagellumNearestDiscretisation{T <: Number} <: Discretisation
    points::NearestDiscretisation

    N_cs::Int  # Number of cross section force points on the tube
    Q_cs::Int  # Number of cross section quadrature points on the tube
    radius::T  # cross-sectional radius of the tube
end

# Forward unknown properties to `points`
@inline function Base.getproperty(f::TubeFlagellumNearestDiscretisation, name::Symbol)
    # Access wrapper's own fields normally
    if name in [:points, :N_cs, :Q_cs, :radius]
        return getfield(f, name)
    end
    # Otherwise, try the underlying pointsretisation
    d = getfield(f, :points)
    return getproperty(d, name)
end

function TubeFlagellumNearestDiscretisation(
    N::Int, 
    N_cs::Int, 
    Q::Int, 
    Q_cs::Int, 
    force_pts::Matrix{T}, 
    quad_pts::Matrix{T}; 
    location=SVector(0.,0.,0.), 
    orientation=I3,
    radius=0.01
) where {T <: Number} 

    points = NearestDiscretisation(
        N*N_cs,
        Q*Q_cs,
        location,
        orientation,
        force_pts,
        zeros(T, size(force_pts)),
        quad_pts,
        nearest_neighbour(force_pts, quad_pts)
    )

    TubeFlagellumNearestDiscretisation(points, N_cs, Q_cs, radius)
end


function nearest_neighbour(force_pts, quad_pts)
    nearest = Int[]
    for x in eachcol(quad_pts)
        d = vec(sum((force_pts .- x).^2, dims=1))
        push!(nearest, argmin(d))
    end
    nearest
end 


function nearest_neighbour!(points::Discretisation)
    points.nearest .= nearest_neighbour(points.force_pts, points.quad_pts)
end
