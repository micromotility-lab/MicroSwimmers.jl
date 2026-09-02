abstract type Model end
abstract type CellBodyModel <: Model end

# Characteristic size, used to turn a target point spacing into a point count. Only models
# that define one of these can size their own discretisation; the rest have to be given N and Q.
_no_size(m, what) = throw(ArgumentError(
    "$(typeof(m)) does not define $what, so it cannot size its own discretisation. " *
    "Pass N and Q explicitly, as `Part(model, N, Q)`."))

surface_area(m::Model) = _no_size(m, "surface_area")
arclength(m::Model)    = _no_size(m, "arclength")

function (m::CellBodyModel)(disc::NearestDiscretisation)
      T = eltype(eltype(disc.force_pts))
      m(disc.force_pts, disc.velocity)
      m(disc.quad_pts)
end

# Nearest spacing helper functions for cell bodies

function hf(model::CellBodyModel; Ns=[113, 2*113+7, 4*113+7, 8*113+7])
    hfs = []
    for N in Ns
        body = CellBody(model, N, 4*N+7)
        push!(hfs, hf(body.points))
    end
    hfs
end

function hq(model::CellBodyModel; Qs=[(2^i)*400 + 7 for i in 1:6])
    hqs = []
    for Q in Qs
        body = CellBody(model, Q ÷ 5 - 7, Q)
        push!(hqs, hq(body.points))
    end
    hqs
end


mutable struct EllipsoidBody{T <: Number} <: CellBodyModel
    a::T
    b::T 
    c::T 
end

(m::EllipsoidBody)(N::Int; pts_fn=fibonacci_ellipsoid) = pts_fn(m.a, m.b, m.c, N)

surface_area(m::EllipsoidBody) = ellipsoid_area(m.a, m.b, m.c)


function (m::EllipsoidBody)(points::Vector{SVector{3,T}}) where {T <: Number}
    points .= m(length(points))
end

function (m::EllipsoidBody)(points::Vector{SVector{3,T}}, velocities::Vector{SVector{3,T}}) where {T <: Number}
    points .= m(length(points))
    velocities .= Ref(zero(SVector{3,T}))
end


mutable struct EllipsoidalGroovedBody{T <: Number} <: CellBodyModel
    a::T
    b::T
    c::T
    g_a::T
    g_b::T
    g_c::T
    groove_center::Vector{T}
    orientation::SMatrix{3,3,T,9}
end

# The groove only removes area, and N, Q are sampling budgets rather than final point counts
# for this model anyway (see SampledBodyModel), so the ungrooved ellipsoid is close enough to
# size a budget from.
surface_area(m::EllipsoidalGroovedBody) = ellipsoid_area(m.a, m.b, m.c)

EllipsoidalGroovedBody(a::T, b::T, c::T, groove_center::Vector{T}; orientation=SMatrix{3,3,T,9}(I)) where {T <: Number} = EllipsoidalGroovedBody(
    a, b, c,
    a, b, c,
    groove_center,
    orientation
)

EllipsoidalGroovedBody(a::T, b::T, c::T, g_a::T, g_b::T, g_c::T, groove_center::Vector{T}; orientation=SMatrix{3,3,T,9}(I)) where {T <: Number} = EllipsoidalGroovedBody(
    a, b, c,
    g_a, g_b, g_c,
    groove_center,
    orientation
)

function (m::EllipsoidalGroovedBody{T})(N::Int; tol=1e-8, pts_fn=fibonacci_ellipsoid) where {T <: Number}  # N is the number of points per ellipse, roughly the total
    centre = SVector{3,T}(m.groove_center)
    axis   = m.orientation * SVector{3,T}(ez)          # groove axis: points out of the mouth
    body_radii   = SVector{3,T}(m.a, m.b, m.c)
    groove_radii = SVector{3,T}(m.g_a, m.g_b, m.g_c)

    ell1 = pts_fn(m.a, m.b, m.c, N)
    ell2 = [centre + m.orientation * p for p in pts_fn(m.g_a, m.g_b, m.g_c, N)]

    body = filter(
        x -> !is_inside_ellipsoid(x, centre, groove_radii, orientation=m.orientation) &&
             dot(axis, x - centre) < tol,
        ell1
    )

    groove = filter(
        x -> is_inside_ellipsoid(x, zero(SVector{3,T}), body_radii) &&
             dot(axis, x - centre) < tol,
        ell2
    )

    [body, groove]
end

# The groove filter discards points, so the cloud size is only known once the surface has
# been sampled — size the discretisation from the clouds rather than from N and Q.
# `weighted` is accepted and ignored: this model filters a Fibonacci cloud and so has no
# area elements to offer. Part() rejects `weighted=true` for it up front via
# supports_quadrature_weights, so only the false case ever reaches here.
function (m::EllipsoidalGroovedBody)(disc::NearestDiscretisation, N::Int, Q::Int; weighted=false, kwargs...)
    T = eltype(eltype(disc.force_pts))
    disc.force_pts = reduce(vcat, m(N; kwargs...))
    disc.quad_pts  = reduce(vcat, m(Q; kwargs...))
    disc.velocity  = [zero(SVector{3,T}) for _ in 1:length(disc.force_pts)]
    disc.quad_wts  = nothing
    disc
end


mutable struct CylindricalGroovedBody{T <: Number} <: CellBodyModel
    a::T
    b::T
    c::T
    g_a::T
    g_b::T
    g_d::T
end


function (m::CylindricalGroovedBody)(N::Int)
    ell = fibonacci_ellipsoid(m.a, m.b, m.c, N)
    cyl = fibonacci_cylinder(m.g_a, m.g_b, m.g_d, N, N ÷ 2) .+ [0., 0., m.g_d]

    body = reduce(
        hcat, 
        filter(x -> !is_inside_cylinder(x, [0., 0., m.g_d],  [m.g_a, m.g_b, m.g_d]), eachcol(ell))
    )
    groove = reduce(
        hcat, 
        filter(x -> is_inside_ellipsoid(x, [0., 0., 0.], [m.a, m.b, m.c]), eachcol(cyl))
    )
    body, groove
end

mutable struct FlatGroovedBody{T <: Number} <: CellBodyModel
    a::T
    b::T
    c::T
    g_a::T
    g_b::T
    groove_floor_center::Vector{T}
end

function (m::FlatGroovedBody)(N::Int)
    ell = fibonacci_ellipsoid(m.a, m.b, m.c, N)

    cyl_sides, _, cyl_bottom = fibonacci_cylinder(m.g_a, m.g_b, 2m.c, N, N ÷ 2)
    cyl_sides .+= m.groove_floor_center
    cyl_bottom .+= m.groove_floor_center

    # cyl_height = m.c - g_depth
    # @info "" 0.5cyl_height
    cyl_center = m.groove_floor_center .+ [0., 0., m.c]
    # cyl_sides, _, cyl_bottom = fibonacci_cylinder(m.g_a, m.g_b, 0.5cyl_height, N, N ÷ 2)
    # @info "" size(cyl_sides)
    # cyl_sides .+= cyl_centerm.
    # @info "" size(cyl_sides)

    # cyl_bottom .+= cyl_bottom
    # return ell, cyl_sides, cyl_bottom, cyl_center
    body = reduce(
        hcat, 
        filter(x -> !is_inside_cylinder(x, cyl_center,  [m.g_a, m.g_b, m.c]), eachcol(ell))
    )

    groove_wall = reduce(
        hcat, 
        filter(x -> is_inside_ellipsoid(x, [0., 0., 0.], [m.a, m.b, m.c]), eachcol(cyl_sides)),
        init = zeros(eltype(cyl_sides), 3, 0)
    )

    groove_floor = reduce(
        hcat, 
        filter(x -> is_inside_ellipsoid(x, [0., 0., 0.], [m.a, m.b, m.c]), eachcol(cyl_bottom))
    )
    [body, groove_wall, groove_floor]
end


