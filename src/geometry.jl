### Ellipsoids

# function fibonacci_ellipsoid(a::T, b::T, c::T, num_points::Int) where {T <: Number}
#     points = Array{T}(undef, 3, num_points)
#     phi = (sqrt(5) + 1) / 2 - 1  # Golden ratio minus 1
#     ga = 2π * phi  

#     for i in 0:num_points-1
#         theta = ga * i  
#         z = 1 - (2i + 1) / num_points  # height
#         r = sqrt(1 - z^2)  # radius of the circle at height z

#         x = a * cos(theta) * r
#         y = b * sin(theta) * r
#         z = c * z

#         points[:, i+1] .= [x, y, z]  
#     end
#     points
# end

function fibonacci_ellipsoid(a::T, b::T, c::T, num_points::Int) where {T <: Number}
    points = Vector{SVector{3,T}}(undef, num_points)
    phi = (sqrt(5) + 1) / 2 - 1  # Golden ratio minus 1
    ga = 2π * phi  

    for i in 0:num_points-1
        theta = ga * i  
        z = 1 - (2i + 1) / num_points  # height
        r = sqrt(1 - z^2)  # radius of the circle at height z

        x = a * cos(theta) * r
        y = b * sin(theta) * r
        z = c * z

        points[i+1] = SVector(x, y, z)  
    end
    points
end


# rejection sampled version
function fibonacci_ellipsoid_rejection(a, b, c, num_points)
    # area element ||∂_θ r × ∂_z r|| as a function of (θ, z)
    function area_element(θ, z)
        r = sqrt(1 - z^2)
        sqrt(a^2*b^2*z^2/(1-z^2 + eps()) + c^2*(a^2*sin(θ)^2 + b^2*cos(θ)^2))
    end
    
    # maximum area element for rejection sampling
    M = maximum(area_element(θ, z) 
                for θ in range(0, 2π, 100), z in range(-0.99, 0.99, 100))
    
    points = zeros(3, num_points)
    i = 0
    while i < num_points
        θ = 2π * rand()
        z = 2*rand() - 1
        r = sqrt(1 - z^2)
        if rand() < area_element(θ, z) / M
            i += 1
            points[:, i] = [a*r*cos(θ), b*r*sin(θ), c*z]
        end
    end
    points
end

is_inside_ellipsoid(x, center, radii; orientation=I3, tol=1e-8) = sum((orientation' * (x .- center) ./ radii) .^ 2) <= 1.0 + tol

# Needs NonlinearSolve so I'm not sure whether to include just for this
#
# function ellipsoid_intersection(;x0=[-0.2, 0.05], p=[(0., 0.), (0.2, 0.1), (0., .1), (0.2 ,0.1)])
#     function f!(res, u, p)
#         X, A, Y, B = p
#         res .= [sum(((u .- X) ./ A).^2) - 1.,  sum(((u .- Y) ./ B).^2) - 1.]
#     end
#     prob = NonlinearProblem(f!, x0, p)
#     sol = solve(prob)
#     [sol[1]; 0.; sol[2]]
# end


## Cylinder

# function fibonacci_cylinder(a::T, b::T, c::T, num_side::Int, num_cap::Int) where {T <: Number}
#     @info "" a b c
#     total_points = num_side + 2 * num_cap
#     points = Array{T}(undef, 3, total_points)

#     # Golden angle
#     phi = (sqrt(5) + 1) / 2 - 1
#     ga = 2π * phi

#     idx = 1  # index into the array

#     ### 1. Lateral surface points ###
#     side = Array{T}(undef, 3, num_side)
#     for i in 0:num_side - 1
#         θ = ga * i
#         # z = c * (1 - (2 * i + 1) / num_side)  # z in [-c, +c]
#         z = c*(i + 0.5) / num_side
#         x = a * cos(θ)
#         y = b * sin(θ)
#         side[:,i+1] = [x, y, z]
#         # points[:, idx] = [x, y, z]
#         # idx += 1
#     end

#     ### 2. Top cap points (z = +c) ###
#     top = Array{T}(undef, 3, num_cap)
#     for i in 0:num_cap - 1
#         θ = ga * i
#         r = sqrt((i + 0.5) / num_cap)  # sqrt to ensure uniform area distribution
#         x = a * r * cos(θ)
#         y = b * r * sin(θ)
#         z = c  # top
#         top[:,i+1] = [x, y, z]
#         # points[:, idx] = [x, y, z]
#         # idx += 1
#     end

#     ### 3. Bottom cap points (z = -c) ###
#     bottom = Array{T}(undef, 3, num_cap)
#     for i in 0:num_cap - 1
#         θ = ga * i
#         r = sqrt((i + 0.5) / num_cap)
#         x = a * r * cos(θ)
#         y = b * r * sin(θ)
#         z = -c  # bottom
#         bottom[:,i+1] = [x, y, z]
#         # points[:, idx] = [x, y, z]
#         # idx += 1
#     end

#     side, top, bottom
# end

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

