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

NearestDiscretisation(N::Int, Q::Int; location=SVector(0.,0.,0.), orientation=I3) = NearestDiscretisation(
    N,
    Q,
    location,
    orientation,
    zeros(3, N),
    zeros(3, N),
    zeros(3, Q),
    zeros(Int, Q)
)

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

struct VanedFlagellumNearestDiscretisation <: Discretisation
    points::NearestDiscretisation

    N_f::Int
    Q_f::Int 

    N_v::Int
    N_start::Int
    N_height::Int

    Q_v::Int
    Q_start::Int
    Q_height::Int
end

# Forward unknown properties to `points`
@inline function Base.getproperty(f::VanedFlagellumNearestDiscretisation, name::Symbol)
    name in (:points, :N_f, :Q_f, :N_v, :N_start, :N_height, :Q_v, :Q_start, :Q_height) ? getfield(f, name) : getproperty(f.points, name)  
end

@inline function Base.setproperty!(f::VanedFlagellumNearestDiscretisation, name::Symbol, value)
    name in (:points, :N_f, :Q_f, :N_v, :N_start, :N_height, :Q_v, :Q_start, :Q_height) ? setfield!(f, name, value) : setproperty!(f.points, name, value)
end

function VanedFlagellumNearestDiscretisation(
    N_f::Int,
    Q_f::Int, 
    N_v::Int, 
    N_start::Int, 
    N_height::Int; 
    location=SVector(0.,0.,0.), 
    orientation=I3
) 
    Q_v = floor(Int, (N_v / N_f) * Q_f) 
    Q_start = ceil(Int, ((N_start-1) / (N_f-1)) * (Q_f-1))
    Q_height = floor(Int, (N_height / N_f) * Q_f)

    points = NearestDiscretisation(
        N_f + N_height*N_v, Q_f + Q_height*Q_v;
        location=location, orientation=orientation
    )

    VanedFlagellumNearestDiscretisation(points, N_f, Q_f, N_v, N_start, N_height, Q_v, Q_start, Q_height)
end


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
    name in (:points, :N_cs, :Q_cs, :radius) ? setfield!(f, name, value) : setproperty!(f.points, name, value)
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
