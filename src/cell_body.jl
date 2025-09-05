abstract type CellBody <: FluidBoundary end

mutable struct EllipsoidBody{T <: Number} <: CellBody
    a::T
    b::T 
    c::T 
    points::Discretisation
end

function EllipsoidBody(
    a::T, 
    b::T, 
    c::T,
    N::Int,
    Q::Int; 
    location::SVector{3,T}=SVector{3,T}(zero(T),zero(T),zero(T)), 
    orientation::SMatrix{3,3,T}=SMatrix{3,3,T}(I)
) where {T <: Number}

    force_pts = fibonacci_ellipsoid(a,b,c,N)
    quad_pts = fibonacci_ellipsoid(a,b,c,Q)

    points = NearestDiscretisation(force_pts, quad_pts; location=location, orientation=orientation)
    EllipsoidBody(a, b, c, points)
end

SphericalBody(;a=0.2, N=27, Q=99, ϵ=0.01) = EllipsoidBody(a, a, a, N, Q)

mutable struct EllipsoidalGroovedBody{T <: Number} <: CellBody 
    a::T
    b::T
    c::T
    groove_center::Vector{T}
    points::Discretisation
end


function EllipsoidalGroovedBody(
    a::T, 
    b::T, 
    c::T, 
    groove_center::Vector{T}, 
    N::Int, 
    Q::Int; 
    location=SVector{3,T}(zero(T),zero(T),zero(T)), 
    orientation=SMatrix{3,3,T}(I)
) where {T <: Number}

    function body_without_groove(num_points)
        ell1 = fibonacci_ellipsoid(a, b, c, num_points)
        reduce(hcat, filter(x -> !is_inside_ellipsoid(x, groove_center, [a; b; c]), eachcol(ell1)))
    end

    function groove_without_body(num_points)
        ell2 = groove_center .+ fibonacci_ellipsoid(a, b, c, num_points)
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

    points = NearestDiscretisation(force_pts, quad_pts, nearest; location=location, orientation=orientation)
 
    EllipsoidalGroovedBody(a, b, c, groove_center, points)
end

