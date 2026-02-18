abstract type FluidBoundary end
update_boundary!(fb::FluidBoundary, t::T) where {T <: Number} = nothing  # In general case, the boundary persists through time without changing

function add_rigid_body_motion!(boundary::FluidBoundary, U::AbstractVector, Ω::AbstractVector)
    boundary.points.velocity .+= U .+ reduce(hcat, cross.(Ref(Ω), eachcol(boundary.points.force_pts)))
end

function add_velocity!(boundary::FluidBoundary, U::AbstractVector)
    boundary.points.velocity .+= U
end

function add_angular_velocity!(boundary::FluidBoundary, Ω::AbstractVector)
    boundary.points.velocity .+= reduce(hcat, cross.(Ref(Ω), eachcol(boundary.points.force_pts)))
end

function reset_velocity!(boundary::FluidBoundary)
    boundary.points.velocity .= 0.0
end

abstract type MicroSwimmer <: FluidBoundary end

# General functions for moving microswimmers

function move_boundary!(S::MicroSwimmer, x0::SVector{3,T}=SVector(0.,0.,0.), B::SMatrix{3,3,T}=I3, t::T=0.0) where {T <: Number}
    update_boundary!(S, t)
    S.points.location = x0
    S.points.orientation = B
end

function move_boundary!(S::MicroSwimmer, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3, T}, t::T) where {T <: Number}
    B = hcat(b1, b2, cross(b1, b2))
    move_boundary!(S, x0, B, t)
end

