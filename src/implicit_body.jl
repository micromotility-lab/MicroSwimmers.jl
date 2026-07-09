abstract type ImplicitBodyModel <: CellBodyModel end

function (m::ImplicitBodyModel)(disc::NearestDiscretisation, N, Q)
      T = eltype(eltype(disc.force_pts))
      raymarch_cloud!(disc.force_pts, x -> implicit(m, x), bounding_radius(m), N; seed=seed(m))
      raymarch_cloud!(disc.quad_pts,  x -> implicit(m, x), bounding_radius(m), Q; seed=seed(m))
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
    orientation::SMatrix{3,3,T}
end

mutable struct ImplicitGroovedEllipsoid{E <: ImplicitEllipsoid} <: ImplicitBodyModel
    body::E
    groove::E
    groove_frame::Frame
end

function ImplicitGroovedEllipsoid(a::T, b::T, c::T, g_a::T, g_b::T, g_c::T, groove_center; orientation=SMatrix{3,3,T}(I)) where {T <: Number}
    ImplicitGroovedEllipsoid{T}(a, b, c, g_a, g_b, g_c, SVector{3,T}(groove_center), SMatrix{3,3,T}(orientation))
end

implicit(m::ImplicitGroovedEllipsoid, x::SVector{3,T}; k=50) where {T <: Number} = smooth_max(ellipsoid(x, m.a, m.b, m.c), -shifted_rotated_ellipsoid(x, m.g_a, m.g_b, m.g_c, m.groove_center, m.orientation), k)
bounding_radius(m::ImplicitGroovedEllipsoid) = 1.2 * maximum([abs(m.a + m.g_a), abs(m.b + m.g_b), abs(m.c + m.g_c)])
seed(m::ImplicitGroovedEllipsoid) = m.groove_center
