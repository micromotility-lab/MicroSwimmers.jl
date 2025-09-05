## CHECK WHETHER THE PROBLEM HAS BEEN SOLVED

total_force(forces, points::Discretisation) = sum(forces[:,n] for n in points.nearest)

function total_torque(forces, points::Discretisation)
    @unpack force_pts, quad_pts, nearest = points
    sum(cross(quad_pts[:,i], forces[:, nearest[i]]) for i in axes(quad_pts,2))
end

function total_force_and_torque(prob::SwimmingProblem)
    check_solved!(prob)
    forces = reshape(prob.force_vals, 3, :)
    points =  prob.points
    total_force(forces, points), total_torque(forces, points)
end

function total_power(forces, points::Discretisation)
    sum(dot(forces[:,n], points.velocity[:,n]) for n in points.nearest)
end

function total_power(prob::SwimmingProblem)
    check_solved!(prob) 
    total_power(reshape(prob.force_vals, 3, :), prob.points)
end

