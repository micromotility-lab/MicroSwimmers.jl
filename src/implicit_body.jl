abstract type ImplicitBodyModel <: CellBodyModel end

# Can this model supply surface quadrature weights? Ray-marching produces them as a
# by-product of the surface hit (see raymarch_cloud!); the Fibonacci-sampled models do not
# know their own area elements, so they can only use the absorbed convention.
supports_quadrature_weights(::Model)             = false
supports_quadrature_weights(::ImplicitBodyModel) = true

# Hook for refreshing quadrature weights after the geometry moves. A rigid Frame preserves
# area, so a rigid body never needs this; it exists so a deforming weighted model can opt in.
# Only called for parts that actually carry weights.
quadrature_weights!(::Model, ::Discretisation) = nothing

function (m::ImplicitBodyModel)(disc::NearestDiscretisation, N, Q; weighted=false)
      T = eltype(eltype(disc.force_pts))
      f = x -> implicit(m, x)
      R = bounding_radius(m)
      # Weights belong to the quadrature cloud only: the force points carry unknowns, not
      # area, so their ray-march weights are discarded.
      raymarch_cloud!(disc.force_pts, f, R, N; seed=seed(m))
      wts = T[]
      raymarch_cloud!(disc.quad_pts, wts, f, R, Q; seed=seed(m))
    #   relax_cloud!(disc.force_pts, f; area=sum(wts), iters=1000)
    #   relax_cloud!(disc.quad_pts, f; area=sum(wts), iters=1000)
      disc.quad_wts = weighted ? wts : nothing
      disc.velocity = [zero(SVector{3,T}) for _ in 1:length(disc.force_pts)]
end

## implicit equations for implicit bodies
ellipsoid(x::SVector{3}, a, b, c)= (x[1]/a)^2 + (x[2]/b)^2 + (x[3]/c)^2 - 1.0
function shifted_rotated_ellipsoid(x::SVector{3}, a, b, c, d::SVector{3}, R::SMatrix{3,3})
    y = R' * (x .- d)  # Apply inverse rotation and shift
    ellipsoid(y, a, b, c)
end

mutable struct ImplicitEllipsoid{T <: Number} <: ImplicitBodyModel
    a::T
    b::T 
    c::T 
end

implicit(m::ImplicitEllipsoid, x::SVector{3,T}) where {T <: Number} = ellipsoid(x, m.a, m.b, m.c)
bounding_radius(m::ImplicitEllipsoid) = 1.2 * maximum([m.a, m.b, m.c])
seed(m::ImplicitEllipsoid) = zero(SVector{3, eltype(m.a)})

mutable struct ImplicitGroovedEllipsoid{T <: Number} <: ImplicitBodyModel
    a::T
    b::T 
    c::T 
    g_a::T
    g_b::T
    g_c::T
    groove_center::SVector{3,T}
    orientation::SMatrix{3,3,T,9}
end

# mutable struct ImplicitGroovedEllipsoid{E <: ImplicitEllipsoid} <: ImplicitBodyModel
#     body::E
#     groove::E
#     groove_frame::Frame
# end

function ImplicitGroovedEllipsoid(a::T, b::T, c::T, g_a::T, g_b::T, g_c::T, groove_center; orientation=SMatrix{3,3,T,9}(I)) where {T <: Number}
    ImplicitGroovedEllipsoid{T}(a, b, c, g_a, g_b, g_c, SVector{3,T}(groove_center), SMatrix{3,3,T,9}(orientation))
end

implicit(m::ImplicitGroovedEllipsoid, x::SVector{3,T}; k=50) where {T <: Number} = smooth_max(ellipsoid(x, m.a, m.b, m.c), -shifted_rotated_ellipsoid(x, m.g_a, m.g_b, m.g_c, m.groove_center, m.orientation), k)
bounding_radius(m::ImplicitGroovedEllipsoid) = 1.2 * maximum([abs(m.a + m.g_a), abs(m.b + m.g_b), abs(m.c + m.g_c)])
seed(m::ImplicitGroovedEllipsoid) = m.groove_center
