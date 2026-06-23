total_force(forces::AbstractVector{<:SVector}) = sum(forces)

function total_torque(forces::AbstractVector{<:SVector}, disc::NearestDiscretisation)
    sum(cross(disc.quad_pts[i], forces[disc.nearest[i]]) for i in eachindex(disc.quad_pts))
end

function total_torque(forces::AbstractVector{<:SVector}, disc::NystromDiscretisation)
    sum(cross(disc.force_pts[i], forces[i]) for i in eachindex(disc.force_pts))
end

function total_force_and_torque(prob::InstantaneousProblem)
    check_solved!(prob)
    forces = get_forces(prob)
    total_force(forces), total_torque(forces, prob.disc)
end


function stresslet_tensor(prob::InstantaneousProblem)
    check_solved!(prob)
    forces = get_forces(prob)
    _stresslet_tensor(forces, prob.disc)
end

function _stresslet_tensor(forces, disc::NearestDiscretisation)
    S_raw = 0.5 * sum(
        forces[disc.nearest[i]] * disc.quad_pts[i]' + disc.quad_pts[i] * forces[disc.nearest[i]]'
        for i in eachindex(disc.quad_pts)
    )
    S_raw - (1/3)*tr(S_raw)*I
end

function _stresslet_tensor(forces, disc::NystromDiscretisation)
    S_raw = 0.5 * sum(
        forces[i] * disc.force_pts[i]' + disc.force_pts[i] * forces[i]'
        for i in eachindex(disc.force_pts)
    )
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
    vels   = prob.disc.velocity
    sum(dot(forces[n], vels[n]) for n in eachindex(forces))
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
