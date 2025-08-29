### Ellipsoids


function fibonacci_ellipsoid(a::Float64, b::Float64, c::Float64, num_points::Int)
    points = Array{Float64}(undef, 3, num_points)
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

is_inside_ellipsoid(x, center, radii) = sum(((x .- center) ./ radii) .^ 2) <= 1


