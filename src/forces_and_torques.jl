total_force(forces, points::Discretisation) = sum(forces[:,n] for n in points.nearest)

function total_torque(forces, points::Discretisation)
    @unpack quad_pts, nearest = points
    sum(cross(quad_pts[:,i], forces[:, nearest[i]]) for i in axes(quad_pts,2))
end

function total_force_and_torque(prob::InstantaneousProblem)
    check_solved!(prob)
    forces = reshape(prob.force_vals, 3, :)
    points =  prob.points
    total_force(forces, points), total_torque(forces, points)
end


function stresslet_tensor(prob::InstantaneousProblem)
    check_solved!(prob)
    @unpack quad_pts, nearest = prob.points
    forces = get_forces(prob)
    S_raw = 0.5*sum(forces[nearest[i]] * quad_pts[:,i]' + quad_pts[:,i] * forces[nearest[i]]' for i in axes(quad_pts, 2))
    S_raw - (1/3)*tr(S_raw)*I
end

function average_stresslet_tensor(prob::InstantaneousProblem; period=1.0, num_ts=30)
    check_solved!(prob)
    Ss = []
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(prob, t)
        solve_problem!(prob)
        push!(Ss, stresslet_tensor(prob))
    end
    sum(Ss) / num_ts
end

function total_power(prob::InstantaneousProblem)
    check_solved!(prob) 
    forces = get_forces(prob)
    vels = get_velocities(prob)
    sum(dot(forces[n], vels[n]) for n in prob.points.nearest)
end

function total_energy_dissipated(prob::SwimmingTrajectoryProblem)
    check_solved!(prob)
    traj = prob.traj
    sprob = prob.swimming_problem

    Es = Float64[]

    for (i, t) in enumerate(traj.t)
        move_boundary!(sprob, traj.x[i], traj.b1[i], traj.b2[i], t) 
        solve_problem!(sprob)
        push!(Es, total_power(sprob))
    end
    # Trapezoidal integration
    sum(0.5 * (Es[i] + Es[i+1]) * (traj.t[i+1] - traj.t[i]) for i in 1:length(traj.t)-1)
end

