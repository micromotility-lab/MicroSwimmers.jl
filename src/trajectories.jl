function average_swimming_velocity(traj::Trajectory )
    (traj.x[end] - traj.x[1]) / (traj.t[end] - traj.t[1])
end

# function continue_periodic_trajectory!(
#     t, x, b1, b2,          # new trajectory arrays to overwrite
#     t_T, x_T, b1_T, b2_T   # the initial period trajectory
# )
#     traj.periodic || error("Cannot continue a non-periodic trajectory")

#     for i in 1:N_periods-1
#         B = hcat(b1[end], b2[end], cross(b1[end], b2[end]))
        
#         ti = t_T .+ t_T[end]*i
#         xi = [B * x for x in x_T] .+ Ref(x[end])
#         b1_i = [B * b1 for b1 in b1_T]
#         b2_i = [B * b2 for b2 in b2_T]
        
#         append!(t, ti[2:end])
#         append!(x, xi[2:end])
#         append!(b1, b1_i[2:end])
#         append!(b2, b2_i[2:end])
#     end
# end


# continue_periodic_trajectory!(traj::Trajectory, N_periods::Int=10) = continue_periodic_trajectory!(
#     traj.t, traj.x, traj.b1, traj.b2,
#     traj.t, traj.x, traj.b
# )

#     traj.periodic || error("Cannot continue a non-periodic trajectory")

#     @unpack t, x, b1, b2 = traj

#     Δx = x[end] - x[1]
#     Δt = t[end]
#     N = length(t)

#     B = hcat(b1[end], b2[end], cross(b1[end], b2[end]))
#     for i in 1:N_periods-1
#         R = B^i

#         for j in 2:N
#             push!(t, t[j] + i*Δt)
#             push!(x, x[j] + R*Δx)
#             push!(b1, R*b1[j])
#             push!(b2, R*b2[j])
#         end
#     end
# end

function continue_periodic_trajectory(traj::Trajectory; N_periods=10)
    T = eltype(traj.t)
    new_ts = T[]
    new_xs = SVector{3,T}[]
    new_b1s = SVector{3,T}[]
    new_b2s = SVector{3,T}[]

    append!(new_ts, traj.t)
    append!(new_xs, traj.x)
    append!(new_b1s, traj.b1)
    append!(new_b2s, traj.b2)

    for i in 1:N_periods-1
        B = hcat(new_b1s[end], new_b2s[end], cross(new_b1s[end], new_b2s[end]))
        @info "" B
        new_t = traj.t .+ traj.t[end]*i
        new_x = [B * x for x in traj.x] .+ Ref(new_xs[end])
        # @info "" new_x new_xs[end]
        new_b1 = [B * b1 for b1 in traj.b1]
        new_b2 = [B * b2 for b2 in traj.b2]
        
        append!(new_ts, new_t[2:end])
        append!(new_xs, new_x[2:end])
        append!(new_b1s, new_b1[2:end])
        append!(new_b2s, new_b2[2:end])
    end
    Trajectory(
        new_ts,
        new_xs,
        new_b1s,
        new_b2s,
        true
    )
end