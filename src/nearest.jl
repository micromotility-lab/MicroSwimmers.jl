mutable struct Configuration{T <: Number}
    location::SVector{3,T}             # current location
    orientation::SMatrix{3,3,T}          # current orientation
    force_pts::Matrix{T}            # reference configuration force points
    velocity::Matrix{T}             # reference configuration velocities (at the force points)
    quad_pts::Matrix{T}             # reference configuration quadrature points
    nearest::Vector{Int}              # length Q
end

# Configuration(location::SVector{3,T},
#               orientation::SMatrix{3,3,T},
#               force_pts::Matrix{T},
#               quad_pts::Matrix{T},
#               velocity::Matrix{T}) where {T} =
#     Configuration{T}(location, orientation, force_pts, quad_pts, velocity)

function nearest_neighbour(force_pts, quad_pts)
    N = Int[]
    for x in eachcol(quad_pts)
        d = vec(sum((force_pts .- x).^2, dims=1))
        push!(N, argmin(d))
    end
    N
end 