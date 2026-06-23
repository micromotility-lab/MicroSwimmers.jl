struct Trajectory{T<:Number}
    t::Vector{T}
    x::Vector{SVector{3,T}}
    b1::Vector{SVector{3,T}}
    b2::Vector{SVector{3,T}}
    periodic::Bool
end

# function centred_trajectory(traj::Trajectory)
#     c = mean(traj.x)
#     Trajectory(
#         traj.t,
#         traj.x .- Ref(c),
#         traj.b1, 
#         traj.b2,
#         traj.periodic
#     )
# end

function translate_trajectory(traj::Trajectory, x)
    Trajectory(
        traj.t,
        traj.x .- Ref(SVector{3}(x)),
        traj.b1, 
        traj.b2,
        traj.periodic
    ) 
end
   
centred_trajectory(traj::Trajectory) = translate_trajectory(traj, mean(traj.x))


function move_boundary!(S::AbstractMicroSwimmer, traj::Trajectory, t_ind::Int=1)
    b1, b2 = traj.b1[t_ind], traj.b2[t_ind]
    S.frame = Frame(traj.x[t_ind], SMatrix{3,3}(hcat(b1, b2, cross(b1, b2))))
    update_boundary!(S, traj.t[t_ind])
end

function continue_periodic_trajectory!(traj::Trajectory, N_periods=10)
    extended = continue_periodic_trajectory(traj, N_periods)
    n0 = length(traj.t)
    append!(traj.t,  @view extended.t[n0+1:end])
    append!(traj.x,  @view extended.x[n0+1:end])
    append!(traj.b1, @view extended.b1[n0+1:end])
    append!(traj.b2, @view extended.b2[n0+1:end])
    traj
end

function average_swimming_velocity(traj::Trajectory)
    (traj.x[end] - traj.x[1]) / (traj.t[end] - traj.t[1])
end

function continue_periodic_trajectory(traj::Trajectory, N_periods=10)
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

function running_mean(traj::Trajectory, radius::Int)
    x = traj.x
    b1 = traj.b1
    b2 = traj.b2

    n = length(traj.t)
    rmean(v) = [mean(v[max(1,i-radius):min(n,i+radius)]) for i in 1:n]
    Trajectory(traj.t, rmean(x), rmean(b1), rmean(b2), traj.periodic)
end
"""
Use principal components of the trajectory to construct some rough helix
parameters to use as an initial guess for a nonlinear least squares fit.
"""
function initial_helix_pars(traj::Trajectory)
    c = mean(traj.x)
    traj_moved = traj.x .- Ref(c)
    U, S, V = svd(reduce(hcat, traj_moved))
    if (S[1] - S[2])/S[1] < 0.1 || (S[1] - S[3])/S[1] < 0.1
        ax_ind=3
    else
        ax_ind=1
    end
    
    v = S[ax_ind]*(V[end,ax_ind] - V[1,ax_ind]) / traj.t[end]
    if v < 0
        @info "here"
        U[:,[ax_ind, 2]] *= -1
        V[:,[ax_ind, 2]] *= -1
    end
    a = U[:,ax_ind]
    @info ""
    u = abs(a[1]) > 0.9 ? ey : ex
    proj = u - dot(a, u)*a
    e1 = proj / norm(proj)
    e2 = cross(a, e1)

    
    X0 = c - dot(c, a)*a
    traj_centred = traj.x .- Ref(X0)
    e1_proj = dot.(Ref(e1), traj_centred)
    e2_proj = dot.(Ref(e2), traj_centred)
    phase = unwrap(atan.(e2_proj, e1_proj))
    A = [ones(length(traj.t)) traj.t]
    ψ, ω = A \ phase

    v = S[ax_ind]*(V[end,ax_ind] - V[1,ax_ind]) / traj.t[end]
    θ = acos(dot(a, ez)) # clamp(a[3], -1.0, 1.0))
    ϕ = sign(a[2])*acos(clamp(a[ax_ind]/(a[ax_ind]^2 + a[2]^2), -1.0, 1.0))
    r = S[2]*(maximum(V[:,2]) - minimum(V[:,2])) / 2
    [X0..., v, ω, θ, ϕ, r, ψ]
end

"""
Fit a helix of the form:

X(t) = [x0,y0,z0] + [sinθcosϕ, sinθsinϕ, cosθ]*v*t + r(cos(ωt + ψ)*e1 + sin(ωt + ψ)*e2)

where e1 and e2 form an orthonormal basis with the helical axis
"""
function fit_helix(traj::Trajectory; N=1, smooth=true, num_t=10)
    traj_extended = N > 1 ? continue_periodic_trajectory(traj, N) : traj
    traj_smooth = smooth ? running_mean(traj_extended, num_t) : traj_extended
    p0 = initial_helix_pars(traj_smooth)
    data = vec(reduce(hcat, traj_smooth.x))
    helix_flat(ts, p) = vec(reduce(hcat, helix(ts, p)))
    fit = curve_fit(helix_flat, traj_smooth.t, data, p0)
    if sqrt(mean(fit.resid .^2)) > 1e-2
        @info "Warning: poor fit with RMSE" sqrt(mean(fit.resid .^2))
    end
    Helix(fit.param...)
end
