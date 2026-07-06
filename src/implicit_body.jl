abstract type ImplicitBodyModel <: CellBodyModel end

function (m::ImplicitBodyModel)(disc::NearestDiscretisation)
      T = eltype(eltype(disc.force_pts))
      raymarch_cloud!(disc.force_pts, x -> implicit(m, x), bounding_radius(m))
      raymarch_cloud!(disc.quad_pts,  x -> implicit(m, x), bounding_radius(m))
      disc.velocity .= Ref(zero(SVector{3,T}))
end

mutable struct ImplicitEllipsoid{T <: Number} <: ImplicitBodyModel
    a::T
    b::T 
    c::T 
end

implicit(m::ImplicitEllipsoid, x::SVector{3,T}) where {T <: Number} = (x[1]/m.a)^2 + (x[2]/m.b)^2 + (x[3]/m.c)^2 - 1.0
bounding_radius(m::ImplicitEllipsoid) = 1.2 * maximum([m.a, m.b, m.c])
