abstract type CellBodyModel end

# calculate the rigid body velocity at pts due to translation U and rotation Ω
(m::CellBodyModel)(pts::AbstractMatrix, U::AbstractVector, Ω::AbstractVector)= U .+ reduce(hcat, cross.(Ref(Ω), eachcol(pts)))

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

(m::EllipsoidBody)(N::Int; pts_fn=fibonacci_ellipsoid) = [pts_fn(m.a, m.b, m.c, N)]

mutable struct EllipsoidalGroovedBody{T <: Number} <: CellBodyModel
    a::T
    b::T
    c::T
    g_a::T
    g_b::T
    g_c::T
    groove_center::Vector{T}
    orientation::SMatrix{3,3,T}
end

EllipsoidalGroovedBody(a::T, b::T, c::T, groove_center::Vector{T}; orientation=I3) where {T <: Number} = EllipsoidalGroovedBody(
    a, b, c,
    a, b, c,
    groove_center,
    orientation
)

EllipsoidalGroovedBody(a::T, b::T, c::T, g_a::T, g_b::T, g_c::T, groove_center::Vector{T}; orientation=I3) where {T <: Number} = EllipsoidalGroovedBody(
    a, b, c,
    g_a, g_b, g_c,
    groove_center,
    orientation
)

function (m::EllipsoidalGroovedBody)(N::Int; tol=1e-8, pts_fn=fibonacci_ellipsoid)  # N is the number of points per ellipse, roughly the total
    ell1 = fibonacci_ellipsoid(m.a, m.b, m.c, N)
    ell2 = m.groove_center .+ m.orientation*pts_fn(m.g_a, m.g_b, m.g_c, N)

    N = m.orientation * ez

    body = reduce(
        hcat, 
        filter(x -> !is_inside_ellipsoid(x, m.groove_center, [m.g_a; m.g_b; m.g_c], orientation=m.orientation) && dot(N, x - m.groove_center) < tol, eachcol(ell1))
    )

    groove = reduce(
        hcat, 
        filter(x -> is_inside_ellipsoid(x, zeros(3), [m.a; m.b; m.c]) && dot(N, x - m.groove_center) < tol, eachcol(ell2)), 
        init=zeros(3,0)
    )

    [body, groove]
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


