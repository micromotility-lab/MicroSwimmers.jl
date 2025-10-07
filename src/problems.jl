abstract type Problem end

abstract type InstantaneousProblem <: Problem end
abstract type DynamicProblem <: Problem end


function check_solved!(prob::InstantaneousProblem)
    if isnothing(prob.force_vals)
        @info "Solving problem"
        solve_problem!(prob)
    end
end

get_force_pts(prob::InstantaneousProblem) = [SVector{3}(pt) for pt in eachcol(prob.points.force_pts)]


mutable struct SwimmingProblem{T<:Number} <: InstantaneousProblem
    microswimmer::MicroSwimmer
    points::Discretisation  # IS THIS NECESSARY?
    eps::T   # regularisation parameter
    mu::T    # viscosity

    lin_prob::LinearProblem
    force_vals::Union{Nothing,Vector{T}} #  keep track of the values of the forces at the force points
    wall::Bool
end

function SwimmingProblem(
    S::MicroSwimmer;
    x0=SVector(0.0, 0.0, 0.0),
    B=I3,
    eps=0.01,
    mu=1.0
)
    T = eltype(S.points.force_pts)
    x0 = SVector{3,T}(x0)
    B = SMatrix{3,3,T}(B)

    N = length(S.points.force_pts)
    points = NearestDiscretisation(
        zeros(T, 3, S.points.N),
        zeros(T, 3, S.points.Q),
        S.points.nearest,
        location=x0,
        orientation=B
    )
    sp = SwimmingProblem(
        S, points, T(eps), T(mu),
        LinearProblem(zeros(T, N + 6, N + 6), zeros(T, N + 6)),
        nothing,
        false
    )
    move_boundary!(sp, x0, B, zero(T))
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

# This gets the total velocity including rigid body dynamics at the force points
function get_velocities(prob::SwimmingProblem)
    U = get_U(prob)
    Ω = get_Ω(prob)

    [U + SVector{3}(vel) for vel in eachcol(prob.points.velocity)] .+ cross.(Ref(Ω), get_force_pts(prob))
end

function update_boundary!(prob::SwimmingProblem, t::T) where {T<:Number}
    update_boundary!(prob.microswimmer, t)
    @views begin
        prob.points.force_pts .= prob.microswimmer.points.force_pts
        prob.points.velocity .= prob.microswimmer.points.velocity
        prob.points.quad_pts .= prob.microswimmer.points.quad_pts
    end
end

function move_boundary!(sp::SwimmingProblem, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::Number) where {T<:Number}
    tT = T(t)
    move_boundary!(sp.microswimmer, x0, B, tT)

    @unpack force_pts, quad_pts, velocity = sp.microswimmer.points
    @views begin
        sp.points.force_pts .= x0 .+ B * sp.microswimmer.points.force_pts
        sp.points.velocity .= B * sp.microswimmer.points.velocity
        sp.points.quad_pts .= x0 .+ B * sp.microswimmer.points.quad_pts
    end
end

function move_boundary!(sp::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3,T}, t::Number) where {T<:Number}
    tT = T(t)
    B = hcat(b1, b2, cross(b1, b2))
    move_boundary!(sp, x0, B, tT)
end

function solve_problem!(sp::SwimmingProblem)
    swimming_matrix!(
        sp.lin_prob.A,
        sp.microswimmer.points.location,
        sp.points.force_pts,
        sp.points.quad_pts,
        sp.microswimmer.points.nearest,
        sp.eps,
        μ=sp.mu
    )

    @views sp.lin_prob.b[1:end-6] .= reshape(sp.points.velocity, :)
    sp.force_vals = solve(sp.lin_prob, MKLLUFactorization())
end

struct Trajectory{T<:Number}
    t::Vector{T}
    x::Vector{SVector{3,T}}
    b1::Vector{SVector{3,T}}
    b2::Vector{SVector{3,T}}
    periodic::Bool
end

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
    mu=1.0
)
    T = eltype(S.points.force_pts)
    x0 = SVector{3,T}(x0)
    B = SMatrix{3,3,T}(B)

    sprob = SwimmingProblem(S; x0=x0, B=B, eps=T(eps), mu=T(mu))

    x0_0 = SVector{3,T}(0, 0, 0)
    b1_0 = SVector{3,T}(1, 0, 0)
    b2_0 = SVector{3,T}(0, 1, 0)
    X0 = vcat(x0_0, b1_0, b2_0)

    function rhs(X, p, t)
        x0 = SVector{3,T}(X[1:3])
        b1 = SVector{3,T}(X[4:6])
        b2 = SVector{3,T}(X[7:9])

        move_boundary!(sprob, x0, b1, b2, t)
        solve_problem!(sprob)
        Ω = get_Ω(sprob)
        vcat(get_U(sprob), cross(Ω, b1), cross(Ω, b2))
        # @views begin
        #     dX[1:3] .= get_U(sprob)
        #     dX[4:6] .= cross(Ω, X[4:6])
        #     dX[7:9] .= cross(Ω, X[7:9])
        # end
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

mutable struct ResistanceProblem{T<:Number} <: InstantaneousProblem
    microswimmer::MicroSwimmer
    points::Discretisation
    eps::T   # regularisation parameter
    mu::T    # viscosity

    lin_prob::LinearProblem
    force_vals::Union{Nothing,Vector{T}}
    wall::Bool
end


function ResistanceProblem(
    microswimmer::MicroSwimmer;
    eps::T=0.01,
    mu::T=1.0,
    wall=false
) where {T<:Number}

    @unpack N, Q, force_pts, quad_pts, velocity, nearest, location, orientation = microswimmer.points

    points = NearestDiscretisation(
        N, Q,
        SVector(0.0, 0.0, 0.0), I3,
        zeros(T, size(force_pts)),
        zeros(T, size(force_pts)),
        zeros(T, size(quad_pts)),
        nearest
    )

    prob = ResistanceProblem(
        microswimmer, points, eps, mu,
        LinearProblem(zeros(T, 3N, 3N), zeros(T, 3N)),
        nothing,
        wall
    )

    update_boundary!(prob, 0.0)
    prob
end
# This gets the total velocity including rigid body dynamics at the force points
get_velocities(prob::ResistanceProblem)  = [SVector{3}(vel) for vel in eachcol(prob.points.velocity)]

function get_forces(prob::ResistanceProblem)
    check_solved!(prob)
    force_vectors = reshape(prob.force_vals, 3, :)
    [SVector{3}(f) for f in eachcol(force_vectors)]
end

function update_boundary!(prob::ResistanceProblem, t::T) where {T<:Number}
    update_boundary!(prob.microswimmer, t)
    @unpack location, orientation, force_pts, quad_pts, velocity = prob.microswimmer.points

    @views begin
        prob.points.force_pts .= location .+ orientation * force_pts
        prob.points.velocity .= orientation * velocity
        prob.points.quad_pts .= location .+ orientation * quad_pts
    end
end

function solve_problem!(prob::ResistanceProblem)
    @unpack lin_prob, points, microswimmer, eps, mu = prob
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

