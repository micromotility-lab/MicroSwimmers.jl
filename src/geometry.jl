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

# Cylinders
function fibonacci_cylinder(a::T, b::T, h::T, num_side::Int, num_cap::Int) where {T <: Number}
    # Golden angle
    phi = (sqrt(5) + 1) / 2 - 1
    ga = 2π * phi

    ### 1. Lateral surface: z ∈ (0, h) ###
    side = Array{T}(undef, 3, num_side)
    for i in 0:num_side - 1
        θ = ga * i
        z = h * (i + 0.5) / num_side          # in (0, h)
        x = a * cos(θ)
        y = b * sin(θ)
        side[:, i+1] = SVector(x, y, z)  # or [x, y, z]
    end

    ### 2. Top cap: z = h ###
    top = Array{T}(undef, 3, num_cap)
    for i in 0:num_cap - 1
        θ = ga * i
        r = sqrt((i + 0.5) / num_cap)
        x = a * r * cos(θ)
        y = b * r * sin(θ)
        z = h
        top[:, i+1] = [x, y, z]
    end

    ### 3. Bottom cap: z = 0 ###
    bottom = Array{T}(undef, 3, num_cap)
    for i in 0:num_cap - 1
        θ = ga * i
        r = sqrt((i + 0.5) / num_cap)
        x = a * r * cos(θ)
        y = b * r * sin(θ)
        z = zero(T)
        bottom[:, i+1] = [x, y, z]
    end

    return side, top, bottom
end


is_inside_cylinder(x, center, dims) = x[3] > center[3] - dims[3] && sum(((x[1:2] .- center[1:2]) ./ dims[1:2]).^2) < 1 

