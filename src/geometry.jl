### Ellipsoids

function fibonacci_ellipsoid(a::T, b::T, c::T, num_points::Int) where {T <: Number}
    points = Array{T}(undef, 3, num_points)
    phi = (sqrt(5) + 1) / 2 - 1  # Golden ratio minus 1
    ga = 2π * phi  

    for i in 0:num_points-1
        theta = ga * i  
        z = 1 - (2i + 1) / num_points  # height
        r = sqrt(1 - z^2)  # radius of the circle at height z

        x = a * cos(theta) * r
        y = b * sin(theta) * r
        z = c * z

        points[:, i+1] .= [x, y, z]  
    end
    points
end


is_inside_ellipsoid(x, center, radii; orientation=I3, tol=1e-8) = sum((orientation' * (x .- center) ./ radii) .^ 2) <= 1.0 + tol