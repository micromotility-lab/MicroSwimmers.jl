abstract type Problem end

abstract type InstantaneousProblem <: Problem end
abstract type DynamicProblem <: Problem end

###########################################################################################
### Instantaneous Problems ################################################################
###########################################################################################

function check_solved!(prob::InstantaneousProblem)
    if isnothing(prob.force_vals)
        @info "Solving problem"
        solve_problem!(prob)
    end
end

function time_collect!(prob::InstantaneousProblem,
    pre_transform!::Function, 
    output::Function, 
    t_final::T, 
    num_t::Int; 
    endpoint=false
) where {T <: Number}
    ts = endpoint ? range(0, t_final, num_t) : range(0, t_final, num_t)[1:end-1]
    X = Vector{Any}(undef, length(ts))
    for (i,t) in enumerate(ts)
        pre_transform!(prob, t)
        solve_problem!(prob)
        X[i] = output(prob)
    end
    X
end

function time_mean!(prob::InstantaneousProblem, pre_transform!::Function, output::Function, t_final::T, num_t::Int; endpoint=false) where {T <: Number}
    X = time_collect!(prob, pre_transform!, output, t_final, num_t; endpoint=endpoint)
    mean(X)
end

function time_mean_std(prob::InstantaneousProblem, pre_transform!::Function, output::Function, t_final::T, num_t::Int; endpoint=false) where {T <: Number}
    X = time_collect!(prob, pre_transform!, output, t_final, num_t; endpoint=endpoint)
    mean(X), std(X)
end

get_force_pts(prob::InstantaneousProblem) = [SVector{3}(pt) for pt in eachcol(prob.points.force_pts)]


function translate_problem!(prob::InstantaneousProblem, x0::AbstractVector)
    prob.points.force_pts .= x0 .+ prob.points.force_pts
    prob.points.quad_pts .= x0 .+ prob.points.quad_pts
    prob.microswimmer.points.location = prob.microswimmer.points.location + SVector{3}(x0)
end

function rotate_problem!(prob::InstantaneousProblem, B::AbstractMatrix)
    prob.points.force_pts .= B * prob.points.force_pts
    prob.points.velocity  .= B * prob.points.velocity
    prob.points.quad_pts  .= B * prob.points.quad_pts
    # prob.microswimmer.points.location = B* prob.microswimmer.points.location
    prob.microswimmer.points.orientation = B * prob.microswimmer.points.orientation
end

mutable struct SwimmingProblem{T<:Number} <: InstantaneousProblem
    microswimmer::MicroSwimmer
    points::Discretisation  # IS THIS NECESSARY?
    eps::T   # regularisation parameter
    mu::T    # viscosity

    lin_prob::LinearProblem
    force_vals::Union{Nothing,Vector{T}} #  keep track of the values of the forces at the force points
    wall::Bool
end

# function SwimmingProblem(
#     S::MicroSwimmer;
#     eps=0.01,
#     mu=1.0,
#     wall=false
# )
#     T = eltype(S.points.force_pts)

#     N = length(S.points.force_pts)
#     points = NearestDiscretisation(
#         zeros(T, 3, S.points.N),
#         zeros(T, 3, S.points.Q),
#         S.points.nearest
#     )
#     sp = SwimmingProblem(
#         S, points, T(eps), T(mu),
#         LinearProblem(zeros(T, N + 6, N + 6), zeros(T, N + 6)),
#         nothing,
#         wall
#     )
#     update_boundary!(sp, zero(T))   
#     sp
# end

_make_points(::Type{NearestDiscretisation}, S) = NearestDiscretisation(
    zeros(eltype(S.points.force_pts), 3, S.points.N),
    zeros(eltype(S.points.force_pts), 3, S.points.Q),
    S.points.nearest
)

_make_points(::Type{NystromDiscretisation}, S) = NystromDiscretisation(
    zeros(eltype(S.points.force_pts), 3, S.points.N),
    zeros(eltype(S.points.velocities), 3, S.points.N)
)

function SwimmingProblem(
    S::MicroSwimmer;
    discretisation::Type{D}=NearestDiscretisation,
    eps=0.01, mu=1.0, wall=false
) where {D}
    T = eltype(S.points.force_pts)
    points = _make_points(discretisation, S)
    N = n_unknowns(points)

    sp = SwimmingProblem(
        S, points, T(eps), T(mu),
        LinearProblem(zeros(T, N+6, N+6), zeros(T, N+6)),
        nothing,
        wall
    )
    update_boundary!(sp, zero(T))
    sp
end

function get_U(prob::SwimmingProblem)
    check_solved!(prob)
    SVector{3}(prob.force_vals[end-5:end-3])
end

function get_Ω(prob::SwimmingProblem)
    check_solved!(prob)
    SVector{3}(prob.force_vals[end-2:end])
end

function get_forces(prob::SwimmingProblem)
    check_solved!(prob)
    force_vectors = reshape(prob.force_vals[1:end-6], 3, :)
    [SVector{3}(f) for f in eachcol(force_vectors)]
end

# Get the total velocity including rigid body dynamics at the force points
function get_velocities(prob::SwimmingProblem)
    U = get_U(prob)
    Ω = get_Ω(prob)
    
    [U + SVector{3}(vel) for vel in eachcol(prob.points.velocity)] .+ cross.(Ref(Ω), get_force_pts(prob) .- Ref(prob.points.location))
end

# Get the total velocity including rigid body dynamics at the quad points
function get_quad_pt_velocities(prob::SwimmingProblem; t=0.0)
    pts = prob.points.quad_pts
    q_pts = zeros(size(pts))
    vs = zeros(size(pts))
    flgt = prob.microswimmer
    q_inds = flgt.quad_pt_indices

    for (i, flagellum) in enumerate(flgt.flagella) 
        Q = flagellum.points.Q
        fl_pts = @view q_pts[:,q_inds[i]:q_inds[i]+Q-1]
        fl_vs = @view vs[:,q_inds[i]:q_inds[i]+Q-1]
        flagellum.model(fl_pts, fl_vs, t)
        vs[:,q_inds[i]:q_inds[i]+Q-1] .= flgt.points.orientation*flagellum.points.orientation*vs[:,q_inds[i]:q_inds[i]+Q-1]
    end

    x0 = prob.microswimmer.points.location
    U = get_U(prob)
    Ω = get_Ω(prob)

    vs = [SVector{3}(v) for v in eachcol(vs)]
    rigid_body_vel = Ref(U) .+ cross.(Ref(Ω), eachcol(pts) .- Ref(x0))
    rigid_body_vel .+ vs
end

function update_boundary!(prob::SwimmingProblem, t::T) where {T <: Number}
    update_boundary!(prob.microswimmer, t)
    @unpack location, orientation, force_pts, quad_pts, velocity = prob.microswimmer.points

    @views begin
        prob.points.force_pts .= location .+ orientation * force_pts
        prob.points.velocity .= orientation * velocity
        prob.points.quad_pts .= location .+ orientation * quad_pts
    end
end

function move_boundary!(prob::SwimmingProblem, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::Number) where {T <: Number}
    tT = T(t)
    move_boundary!(prob.microswimmer, x0, B, tT)

    @unpack force_pts, quad_pts, velocity = prob.microswimmer.points
    @views begin
        prob.points.force_pts .= x0 .+ B * prob.microswimmer.points.force_pts
        prob.points.velocity .= B * prob.microswimmer.points.velocity
        prob.points.quad_pts .= x0 .+ B * prob.microswimmer.points.quad_pts
    end
end

function move_boundary!(prob::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3,T}, t::Number) where {T<:Number}
    tT = T(t)
    B = hcat(b1, b2, cross(b1, b2))
    move_boundary!(prob, x0, B, tT)
end

function solve_problem!(prob::SwimmingProblem)
    swimming_matrix!(
        prob.lin_prob.A,
        prob.microswimmer.points.location,
        # prob.points.force_pts,
        # prob.points.quad_pts,
        # prob.microswimmer.points.nearest,
        prob.points,
        prob.eps,
        μ=prob.mu,
        wall=prob.wall
    )

    @views prob.lin_prob.b[1:end-6] .= reshape(prob.points.velocity, :)
    prob.force_vals = solve(prob.lin_prob, MKLLUFactorization())
end

# Check body boundary conditions at quad_pts (fluid velocity should equal rigid body velocity)
function check_body_boundary_conditions(prob::SwimmingProblem)
    body_pts = prob.points.quad_pts[:,1:prob.microswimmer.body.points.Q]
    x0 = prob.microswimmer.points.location
    U = get_U(prob)
    Ω = get_Ω(prob)

    rigid_body_vel = Ref(U) .+ cross.(Ref(Ω), eachcol(body_pts) .- Ref(x0))
    u = FluidVelocity(prob)

    resid = norm.(u.(eachcol(body_pts)) .- rigid_body_vel)
    median(resid), maximum(resid)
end

# Check all boundary conditions at quad pts, using nearest to approximate the velocities at quad points
function check_boundary_conditions(prob::SwimmingProblem; t=0.0)
    update_boundary!(prob, t)
    solve_problem!(prob)
    pts = prob.points.quad_pts
    u = FluidVelocity(prob)

    q_vs = get_quad_pt_velocities(prob, t=t)
    V_scale = quantile(norm.(q_vs), 0.95) # median(norm.(q_vs))
    @info "" V_scale
    resid = norm.(u.(eachcol(pts)) .- q_vs) ./ V_scale
    resid
    # median(resid), findmax(resid)
end

mutable struct ResistanceProblem{T<:Number} <: InstantaneousProblem
    boundary::FluidBoundary
    points::Discretisation
    eps::T   # regularisation parameter
    mu::T    # viscosity

    lin_prob::LinearProblem
    force_vals::Union{Nothing,Vector{T}}
    wall::Bool
end


function ResistanceProblem(
    boundary::FluidBoundary;
    eps::T=0.01,
    mu::T=1.0,
    wall=false
) where {T<:Number}

    @unpack N, Q, force_pts, quad_pts, velocity, nearest, location, orientation = boundary.points

    points = NearestDiscretisation(
        N, Q,
        SVector(0.0, 0.0, 0.0), I3,
        zeros(T, size(force_pts)),
        zeros(T, size(force_pts)),
        zeros(T, size(quad_pts)),
        nearest
    )

    prob = ResistanceProblem(
        boundary, points, eps, mu,
        LinearProblem(zeros(T, 3N, 3N), zeros(T, 3N)),
        nothing,
        wall
    )

    update_boundary!(prob, 0.0)
    prob
end

get_velocities(prob::ResistanceProblem)  = [SVector{3}(vel) for vel in eachcol(prob.points.velocity)]

function get_forces(prob::ResistanceProblem)
    check_solved!(prob)
    force_vectors = reshape(prob.force_vals, 3, :)
    [SVector{3}(f) for f in eachcol(force_vectors)]
end

# Get velocities at the quad points
function get_quad_pt_velocities(prob::ResistanceProblem; t=0.0)
    pts = prob.points.quad_pts
    q_pts = zeros(size(pts))
    vs = zeros(size(pts))
    flgt = prob.boundary
    q_inds = flgt.quad_pt_indices

    for (i, flagellum) in enumerate(flgt.flagella) 
        Q = flagellum.points.Q
        fl_pts = @view q_pts[:,q_inds[i]:q_inds[i]+Q-1]
        fl_vs = @view vs[:,q_inds[i]:q_inds[i]+Q-1]
        flagellum.model(fl_pts, fl_vs, t)
        vs[:,q_inds[i]:q_inds[i]+Q-1] .= flgt.points.orientation*flagellum.points.orientation*vs[:,q_inds[i]:q_inds[i]+Q-1]
    end
    [SVector{3}(v) for v in eachcol(vs)]
end

function update_boundary!(prob::ResistanceProblem, t::T) where {T<:Number}
    update_boundary!(prob.boundary, t)
    @unpack location, orientation, force_pts, quad_pts, velocity = prob.boundary.points

    @views begin
        prob.points.force_pts .= location .+ orientation * force_pts
        prob.points.velocity .= orientation * velocity
        prob.points.quad_pts .= location .+ orientation * quad_pts
    end
end

function add_rigid_body_motion!(prob::ResistanceProblem, U::AbstractVector, Ω::AbstractVector)
    prob.points.velocity .+= U .+ reduce(hcat, cross.(Ref(Ω), eachcol(prob.points.force_pts)))
end

function solve_problem!(prob::ResistanceProblem)
    @unpack lin_prob, points, boundary, eps, mu = prob
    resistance_matrix!(
        lin_prob.A,
        points.force_pts,
        points.quad_pts,
        points.nearest,
        eps;
        μ=mu,
    )
    lin_prob.b .= reshape(points.velocity, :)
    prob.force_vals = solve(lin_prob, MKLLUFactorization())
end

function check_boundary_conditions(prob::ResistanceProblem; t=0.0)
    pts = prob.points.quad_pts
    vs = get_quad_pt_velocities(prob; t=t)
    V_scale = quantile(norm.(vs), 0.95) # median(norm.(q_vs))
    u = FluidVelocity(prob)
    resid = norm.(u.(eachcol(pts)) .- vs) ./ V_scale
    @info "" V_scale
    resid
    # median(resid), maximum(resid)
end


###########################################################################################
### Dynamic Problems ######################################################################
###########################################################################################


mutable struct SwimmingTrajectoryProblem <: DynamicProblem
    swimming_problem::SwimmingProblem
    ode_prob::ODEProblem
    traj::Union{Nothing,Trajectory}
end

function SwimmingTrajectoryProblem(
    S::MicroSwimmer;
    x0=SVector(0.0, 0.0, 0.0),
    B=I3,
    t_final=20.0,
    saveat=0.05,
    eps=0.01,
    mu=1.0,
    wall=false
)
    T = eltype(S.points.force_pts)

    sprob = SwimmingProblem(S; eps=T(eps), mu=T(mu), wall=wall)

    x0_0 = SVector{3,T}(x0)
    b1_0 = SVector{3,T}(B[:,1])
    b2_0 = SVector{3,T}(B[:,2])
    X0 = vcat(x0_0, b1_0, b2_0)

    function rhs(X, p, t)
        x0 = SVector{3,T}(X[1:3])
        b1 = SVector{3,T}(X[4:6])
        b2 = SVector{3,T}(X[7:9])

        move_boundary!(sprob, x0, b1, b2, t)
        solve_problem!(sprob)
        Ω = get_Ω(sprob)
        vcat(get_U(sprob), cross(Ω, b1), cross(Ω, b2))
    end

    SwimmingTrajectoryProblem(
        sprob,
        ODEProblem(rhs, X0, (T(0), T(t_final)), saveat=T(saveat)),
        nothing
    )
end

function solve_problem!(prob::SwimmingTrajectoryProblem; method=Tsit5(), periodic=false)
    sol = solve(prob.ode_prob, method)
    u = sol.u

    x = [SVector{3}(u[i][1:3]) for i in eachindex(u)]
    b1 = [SVector{3}(u[i][4:6]) for i in eachindex(u)]
    b2 = [SVector{3}(u[i][7:9]) for i in eachindex(u)]
    prob.traj = Trajectory(sol.t, x, b1, b2, periodic)
end

function get_sol!(prob::SwimmingTrajectoryProblem)
    solve(prob.ode_prob, Tsit5()).u
end

function check_solved!(prob::SwimmingTrajectoryProblem)
    if isnothing(prob.traj)
        @info "Solving swimming trajectory problem"
        solve_problem!(prob)
    end
end

function move_boundary!(prob::SwimmingTrajectoryProblem, t_ind::Int)
    check_solved!(prob)
    traj = prob.traj
    move_boundary!(prob.swimming_problem, traj.x[t_ind], traj.b1[t_ind], traj.b2[t_ind], traj.t[t_ind])
end

# Change to ParticleTrajectoryResistanceProblem and do the same for swimming

mutable struct ParticleTrajectoryProblem{T<:Number} <: Problem
    resistance_problem::ResistanceProblem{T}
    ode_prob::ODEProblem
    t::Union{Nothing,Vector{T}}
    trajectories::Union{Nothing,Matrix{T}}
end

function ParticleTrajectoryProblem(
    microswimmer::MicroSwimmer;
    x=-5.0,
    ys=range(-4.0, 4.0, 6),
    zs=range(0.2, 3.2, 6),
    t_final::T=20.0,
    saveat::T=0.05,
    eps=0.01,
    mu=1.0
) where {T<:Number}
    num_particles = length(ys)*length(zs)
    rprob = ResistanceProblem(microswimmer; eps=eps, mu=mu)
    A = zeros(3 * num_particles, 3rprob.points.N)

    function rhs!(dX, X, p, t)
        update_boundary!(rprob, t)
        solve_problem!(rprob)

        resistance_matrix!(
            A,
            reshape(X, 3, num_particles),
            rprob.points.quad_pts,
            rprob.points.nearest,
            eps;
            μ=mu
        )
        dX .= A * rprob.force_vals
    end

    tspan = (0.0, t_final)

    x0 = reduce(vcat, [[x, y, z] for y in ys, z in zs])

    ParticleTrajectoryProblem(
        rprob,
        ODEProblem(rhs!, x0, tspan, saveat=saveat),
        nothing,
        nothing
    )
end

function solve_problem!(prob::ParticleTrajectoryProblem; method=Tsit5())
    sol = solve(prob.ode_prob, method)
    prob.t = sol.t
    prob.trajectories = reduce(hcat, sol.u)
end

