struct Trajectory{T <: Number}
    t::Vector{T}
    x::Vector{SVector{3,T}}
    b1::Vector{SVector{3,T}}
    b2::Vector{SVector{3,T}}
    periodic::Bool
end

function Trajectory(prob::DynamicSwimmingProblem, periodic::Bool)
    if isnothing(prob.sol)
        solve_problem!(prob)
    end
    t = prob.sol.t
    u = prob.sol.u

    x  = [SVector{3}(u[i][1:3]) for i in eachindex(u)]
    b1 = [SVector{3}(u[i][4:6]) for i in eachindex(u)]
    b2 = [SVector{3}(u[i][7:9]) for i in eachindex(u)]
    Trajectory(t, x, b1, b2, periodic)
end

function swimming_velocity(traj::Trajectory)
    (traj.x[end] - traj.x[1]) / (traj.t[end] - traj.t[1])
end



function continue_periodic_trajectory!(traj::Trajectory; N_periods::Int=10)
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