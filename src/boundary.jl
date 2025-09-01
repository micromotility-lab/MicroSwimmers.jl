abstract type FluidBoundary end
update_boundary!(fb::FluidBoundary, t::T) where {T <: Number} = nothing  # In general case, the boundary persists through time without changing