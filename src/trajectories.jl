function average_swimming_velocity(traj::Trajectory; periods=100)
    if traj.periodic
        continue_periodic_trajectory!(traj, periods)
    end
    (traj.x[end] - traj.x[1]) / (traj.t[end] - traj.t[1])
end


function continue_periodic_trajectory!(traj::Trajectory, N_periods::Int=10)
    traj.periodic || error("Cannot continue a non-periodic trajectory")

    @unpack t, x, b1, b2 = traj

    Δx = x[end] - x[1]
    Δt = t[end]
    N = length(t)

    B = hcat(b1[end], b2[end], cross(b1[end], b2[end]))
    for i in 1:N_periods-1
        R = B^i

        for j in 2:N
            push!(t, t[j] + i*Δt)
            push!(x, R*x[j] + i*Δx)
            push!(b1, R*b1[j])
            push!(b2, R*b2[j])
        end
    end
end