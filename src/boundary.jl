abstract type FluidBoundary end
update_boundary!(fb::FluidBoundary, t::T) where {T <: Number} = nothing  # In general case, the boundary persists through time without changing

abstract type AbstractMicroSwimmer <: FluidBoundary end

# General functions for moving microswimmers
function move_boundary!(S::AbstractMicroSwimmer, x0=SVector(0.,0.,0.), B=I3, t=0.0)
    update_boundary!(S, t)
    S.points.location = SVector{3}(x0)
    S.points.orientation = SMatrix{3,3}(B)
end

function move_boundary!(S::AbstractMicroSwimmer, x0::AbstractVector, b1::AbstractVector, b2::AbstractVector, t=0.0)
    x0 = SVector{3}(x0)
    b1 = SVector{3}(b1)
    b2 = SVector{3}(b2)
    B = hcat(b1, b2, cross(b1, b2))
    move_boundary!(S, x0, B, t)
end

