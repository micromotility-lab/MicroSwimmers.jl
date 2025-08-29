abstract type CellBody <: FluidBoundary end

mutable struct EllipsoidBody{T <: Number} <: CellBody
    a::T
    b::T 
    c::T 
    N::Int 
    Q::Int
    config::Configuration
    ϵ::T # regularisation parameter
end

function EllipsoidBody(
    a::T, 
    b::T, 
    c::T, 
    N::Int, 
    Q::Int, 
    ϵ::T;
    location=SVector{3,T}(zero(T),zero(T),zero(T)), 
    orientation=SMatrix{3,3,T}(I)
) where {T <: Number}
    force_pts = fibonacci_ellipsoid(a,b,c,N)
    quad_pts = fibonacci_ellipsoid(a,b,c,Q)
    nearest = nearest_neighbour(force_pts, quad_pts)

    EllipsoidBody(
        a,
        b,
        c,
        N,
        Q,
        Configuration(location, orientation, force_pts, zeros(T, 3, N), quad_pts, nearest),
        ϵ
    )
end

SphericalBody(;a=0.2, N=27, Q=99, ϵ=0.01) = EllipsoidBody(a, a, a, N, Q, ϵ)

mutable struct EllipsoidalGroovedBody{T <: Number} <: CellBody 
    a::T
    b::T
    c::T
    center_ell2::Vector{T}
    N::Int 
    Q::Int
    config::Configuration
    ϵ::T # regularisation parameter
end


function EllipsoidalGroovedBody(
    a::T, 
    b::T, 
    c::T, 
    center_ell2::Vector{T}, 
    N::Int, 
    Q::Int, 
    ϵ::T;
    location=SVector{3,T}(zero(T),zero(T),zero(T)), 
    orientation=SMatrix{3,3,T}(I)
) where {T <: Number}

    function body_without_groove(num_points)
        ell1 = fibonacci_ellipsoid(a, b, c, num_points)
        reduce(hcat, filter(x -> !is_inside_ellipsoid(x, center_ell2, [a; b; c]), eachcol(ell1)))
    end

    function groove_without_body(num_points)
        ell2 = center_ell2 .+ fibonacci_ellipsoid(a, b, c, num_points)
        reduce(hcat, filter(x -> is_inside_ellipsoid(x, zeros(T,3), [a; b; c]), eachcol(ell2)), init=zeros(T,3,0))
    end

    body_force_pts = body_without_groove(N)
    body_quad_pts  = body_without_groove(Q)
    nearest_body   = nearest_neighbour(body_force_pts, body_quad_pts)

    groove_force_pts = groove_without_body(N)
    groove_quad_pts  = groove_without_body(Q)
    nearest_groove   = nearest_neighbour(groove_force_pts, groove_quad_pts)

    force_pts = [body_force_pts groove_force_pts]
    quad_pts  = [body_quad_pts groove_quad_pts]
    nearest   = [nearest_body; size(body_force_pts, 2) .+ nearest_groove]
 
    EllipsoidalGroovedBody(
        a,
        b,
        c,
        center_ell2,
        size(force_pts, 2),
        size(quad_pts, 2),
        Configuration(
            location, 
            orientation, 
            force_pts, 
            zeros(T, 3, size(force_pts,2)), 
            quad_pts, 
            nearest
        ),
        ϵ
    )
end

