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

is_inside_ellipsoid(x, center, radii) = sum(((x .- center) ./ radii) .^ 2) <= 1

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

