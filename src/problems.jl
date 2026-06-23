abstract type Problem end

abstract type InstantaneousProblem <: Problem end
abstract type DynamicProblem <: Problem end

###########################################################################################
### Generic helpers #######################################################################
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
    for (i, t) in enumerate(ts)
        pre_transform!(prob, t)
        solve_problem!(prob)
        X[i] = output(prob)
    end
    X
end

function time_mean!(prob::InstantaneousProblem, pre_transform!::Function, output::Function, t_final::T, num_t::Int; endpoint=false) where {T <: Number}
    mean(time_collect!(prob, pre_transform!, output, t_final, num_t; endpoint=endpoint))
end

function time_mean_std(prob::InstantaneousProblem, pre_transform!::Function, output::Function, t_final::T, num_t::Int; endpoint=false) where {T <: Number}
    X = time_collect!(prob, pre_transform!, output, t_final, num_t; endpoint=endpoint)
    mean(X), std(X)
end

# # Old API helpers — tied to Matrix-based discretisation
# get_force_pts(prob::InstantaneousProblem) = [SVector{3}(pt) for pt in eachcol(prob.points.force_pts)]
#
# function translate_problem!(prob::InstantaneousProblem, x0::AbstractVector)
#     prob.points.force_pts .= x0 .+ prob.points.force_pts
#     prob.points.quad_pts  .= x0 .+ prob.points.quad_pts
#     prob.microswimmer.points.location = prob.microswimmer.points.location + SVector{3}(x0)
# end
#
# function rotate_problem!(prob::InstantaneousProblem, B::AbstractMatrix)
#     prob.points.force_pts .= B * prob.points.force_pts
#     prob.points.velocity  .= B * prob.points.velocity
#     prob.points.quad_pts  .= B * prob.points.quad_pts
#     prob.microswimmer.points.orientation = B * prob.microswimmer.points.orientation
# end

###########################################################################################
### SwimmingProblem #######################################################################
###########################################################################################

mutable struct SwimmingProblem{MS <: MicroSwimmer, D <: Discretisation, T <: Number, K <: Kernel} <: InstantaneousProblem
    microswimmer::MS
    disc::D
    force_pt_indices::Vector{Int}
    quad_pt_indices::Vector{Int}
    mu::T
    lin_prob::LinearProblem
    force_vals::Union{Nothing, Vector{T}}
    kernel::K
end

function SwimmingProblem(ms::MicroSwimmer{<:Part{<:Model, <:NearestDiscretisation}}; mu=1.0, eps=0.1)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    nq_sizes = [nq(p.disc) for p in ms.parts]
    N = sum(nf_sizes); Q = sum(nq_sizes)
    force_pt_indices = cumsum([1; nf_sizes[1:end-1]])
    quad_pt_indices  = cumsum([1; nq_sizes[1:end-1]])
    prob = SwimmingProblem(
        ms,
        NearestDiscretisation(N, Q),
        force_pt_indices,
        quad_pt_indices,
        Float64(mu),
        LinearProblem(zeros(3N+6, 3N+6), zeros(3N+6)),
        nothing,
        RegStokeslet(eps)
    )
    gather_nearest!(prob)
    update_boundary!(prob, 0.0)
    prob
end

function SwimmingProblem(ms::MicroSwimmer{<:Part{<:Model, <:NystromDiscretisation}}; mu=1.0, eps=0.1)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    N        = sum(nf_sizes)
    indices  = cumsum([1; nf_sizes[1:end-1]])
    prob = SwimmingProblem(
        ms,
        NystromDiscretisation(N),
        indices,
        indices,
        Float64(mu),
        LinearProblem(zeros(3N+6, 3N+6), zeros(3N+6)),
        nothing,
        RegStokeslet(eps)
    )
    update_boundary!(prob, 0.0)
    prob
end

function gather_nearest!(prob::SwimmingProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc, force_pt_indices, quad_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part   = microswimmer.parts[i]
        foff   = force_pt_indices[i] - 1
        qstart = quad_pt_indices[i]
        nq_i   = length(part.disc.nearest)
        @views disc.nearest[qstart:qstart+nq_i-1] .= part.disc.nearest .+ foff
    end
end

function gather!(prob::SwimmingProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc, force_pt_indices, quad_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part      = microswimmer.parts[i]
        lab_frame = microswimmer.frame * part.frame
        fstart    = force_pt_indices[i]
        qstart    = quad_pt_indices[i]
        nf_i      = nf(part.disc)
        nq_i      = nq(part.disc)
        @views disc.force_pts[fstart:fstart+nf_i-1] .= lab_frame.(part.disc.force_pts)
        @views disc.velocity[fstart:fstart+nf_i-1]  .= Ref(lab_frame.orientation) .* part.disc.velocity
        @views disc.quad_pts[qstart:qstart+nq_i-1]  .= lab_frame.(part.disc.quad_pts)
    end
end

function gather!(prob::SwimmingProblem{<:Any, <:NystromDiscretisation})
    @unpack microswimmer, disc, force_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part      = microswimmer.parts[i]
        lab_frame = microswimmer.frame * part.frame
        fstart    = force_pt_indices[i]
        nf_i      = nf(part.disc)
        @views disc.force_pts[fstart:fstart+nf_i-1] .= lab_frame.(part.disc.force_pts)
        @views disc.velocity[fstart:fstart+nf_i-1]  .= Ref(lab_frame.orientation) .* part.disc.velocity
    end
end

get_force_pts(prob::SwimmingProblem) = prob.disc.force_pts

function get_U(prob::SwimmingProblem)
    check_solved!(prob)
    fv = prob.force_vals
    SVector{3}(fv[end-5], fv[end-4], fv[end-3])
end

function get_Ω(prob::SwimmingProblem)
    check_solved!(prob)
    fv = prob.force_vals
    SVector{3}(fv[end-2], fv[end-1], fv[end])
end

function get_forces(prob::SwimmingProblem)
    check_solved!(prob)
    fv = prob.force_vals
    N  = nf(prob.disc)
    [SVector{3}(fv[3i-2], fv[3i-1], fv[3i]) for i in 1:N]
end

function update_boundary!(prob::SwimmingProblem, t::Number)
    update_boundary!(prob.microswimmer, t)
    gather!(prob)
end

function move_boundary!(prob::SwimmingProblem, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::Number) where T
    prob.microswimmer.frame = Frame(x0, B)
    update_boundary!(prob, T(t))
end

function move_boundary!(prob::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3,T}, t::Number) where T
    move_boundary!(prob, x0, SMatrix{3,3,T}(hcat(b1, b2, cross(b1, b2))), t)
end

function solve_problem!(prob::SwimmingProblem)
    @unpack lin_prob, disc, kernel, mu, microswimmer = prob
    assemble_swimming!(lin_prob.A, microswimmer.frame.location, disc, kernel; μ=mu)
    N3 = 3 * nf(disc)
    T  = eltype(eltype(disc.velocity))
    @views lin_prob.b[1:N3] .= reinterpret(T, disc.velocity)
    @views lin_prob.b[N3+1:end] .= zero(T)
    prob.force_vals = solve(lin_prob, MKLLUFactorization()).u
end

# # Old SwimmingProblem — Flagellate/CellBody/Flagellum API (Matrix-based discretisation)
#
# mutable struct SwimmingProblem{T<:Number} <: InstantaneousProblem
#     microswimmer::AbstractMicroSwimmer
#     points::Discretisation
#     eps::T
#     mu::T
#     lin_prob::LinearProblem
#     force_vals::Union{Nothing,Vector{T}}
#     wall::Bool
# end
#
# _make_points(::Type{NearestDiscretisation}, S) = NearestDiscretisation(
#     zeros(eltype(S.points.force_pts), 3, S.points.N),
#     zeros(eltype(S.points.force_pts), 3, S.points.Q),
#     S.points.nearest
# )
#
# _make_points(::Type{NystromDiscretisation}, S) = NystromDiscretisation(
#     zeros(eltype(S.points.force_pts), 3, S.points.N),
#     zeros(eltype(S.points.velocities), 3, S.points.N)
# )
#
# function SwimmingProblem(S::AbstractMicroSwimmer; discretisation=NearestDiscretisation,
#                          eps=0.01, mu=1.0, wall=false)
#     T = eltype(S.points.force_pts)
#     points = _make_points(discretisation, S)
#     N = n_unknowns(points)
#     sp = SwimmingProblem(S, points, T(eps), T(mu),
#                          LinearProblem(zeros(T, N+6, N+6), zeros(T, N+6)), nothing, wall)
#     update_boundary!(sp, zero(T))
#     sp
# end
#
# function get_U(prob::SwimmingProblem) ... end
# function get_Ω(prob::SwimmingProblem) ... end
# function get_forces(prob::SwimmingProblem) ... end
# function get_velocities(prob::SwimmingProblem) ... end
# function get_quad_pt_velocities(prob::SwimmingProblem; t=0.0) ... end
# function update_boundary!(prob::SwimmingProblem, t) ... end
# function move_boundary!(prob::SwimmingProblem, x0, B, t) ... end
# function move_boundary!(prob::SwimmingProblem, x0, b1, b2, t) ... end
# function solve_problem!(prob::SwimmingProblem) ... end
# function check_body_boundary_conditions(prob::SwimmingProblem) ... end
# function check_boundary_conditions(prob::SwimmingProblem; t=0.0) ... end

###########################################################################################
### ResistanceProblem #####################################################################
###########################################################################################

mutable struct ResistanceProblem{MS <: AbstractMicroSwimmer, D <: Discretisation, T <: Number, K <: Kernel, L <: LinearProblem} <: InstantaneousProblem
    microswimmer::MS
    disc::D
    force_pt_indices::Vector{Int}
    quad_pt_indices::Vector{Int}
    mu::T
    lin_prob::L
    force_vals::Union{Nothing, Vector{T}}
    kernel::K
end

function ResistanceProblem(ms::MicroSwimmer{<:Part{<:Model, <:NearestDiscretisation}}; mu=1.0, eps=0.1)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    nq_sizes = [nq(p.disc) for p in ms.parts]
    N = sum(nf_sizes); Q = sum(nq_sizes)
    force_pt_indices = cumsum([1; nf_sizes[1:end-1]])
    quad_pt_indices  = cumsum([1; nq_sizes[1:end-1]])
    prob = ResistanceProblem(
        ms,
        NearestDiscretisation(N, Q),
        force_pt_indices,
        quad_pt_indices,
        Float64(mu),
        LinearProblem(zeros(3N, 3N), zeros(3N)),
        nothing,
        RegStokeslet(eps)
    )
    gather_nearest!(prob)
    prob
end

function ResistanceProblem(ms::MicroSwimmer{<:Part{<:Model, <:NystromDiscretisation}}; mu=1.0, eps=0.1)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    N        = sum(nf_sizes)
    indices  = cumsum([1; nf_sizes[1:end-1]])
    ResistanceProblem(
        ms,
        NystromDiscretisation(N),
        indices,
        indices,
        Float64(mu),
        LinearProblem(zeros(3N, 3N), zeros(3N)),
        nothing,
        RegStokeslet(eps)
    )
end

function gather_nearest!(prob::ResistanceProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc, force_pt_indices, quad_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part   = microswimmer.parts[i]
        foff   = force_pt_indices[i] - 1
        qstart = quad_pt_indices[i]
        nq_i   = length(part.disc.nearest)
        @views disc.nearest[qstart:qstart+nq_i-1] .= part.disc.nearest .+ foff
    end
end

function gather!(prob::ResistanceProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc, force_pt_indices, quad_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part      = microswimmer.parts[i]
        lab_frame = microswimmer.frame * part.frame
        fstart    = force_pt_indices[i]
        qstart    = quad_pt_indices[i]
        nf_i      = nf(part.disc)
        nq_i      = nq(part.disc)
        @views disc.force_pts[fstart:fstart+nf_i-1] .= lab_frame.(part.disc.force_pts)
        @views disc.velocity[fstart:fstart+nf_i-1]  .= Ref(lab_frame.orientation) .* part.disc.velocity
        @views disc.quad_pts[qstart:qstart+nq_i-1]  .= lab_frame.(part.disc.quad_pts)
    end
end

function gather!(prob::ResistanceProblem{<:Any, <:NystromDiscretisation})
    @unpack microswimmer, disc, force_pt_indices = prob
    for i in eachindex(microswimmer.parts)
        part      = microswimmer.parts[i]
        lab_frame = microswimmer.frame * part.frame
        fstart    = force_pt_indices[i]
        nf_i      = nf(part.disc)
        @views disc.force_pts[fstart:fstart+nf_i-1] .= lab_frame.(part.disc.force_pts)
        @views disc.velocity[fstart:fstart+nf_i-1]  .= Ref(lab_frame.orientation) .* part.disc.velocity
    end
end

function add_rigid_body_motion!(prob::ResistanceProblem, U, Ω)
    prob.disc.velocity .= U .+ cross.(Ref(Ω), prob.disc.force_pts)
end

function get_forces(prob::ResistanceProblem)
    check_solved!(prob)
    fv = prob.force_vals
    N  = nf(prob.disc)
    [SVector{3}(fv[3i-2], fv[3i-1], fv[3i]) for i in 1:N]
end

function solve_problem!(prob::ResistanceProblem)
    @unpack lin_prob, disc, kernel, mu = prob
    gather!(prob)
    assemble!(lin_prob.A, disc, kernel; μ=mu)
    lin_prob.b .= reinterpret(eltype(eltype(disc.velocity)), disc.velocity)
    prob.force_vals = solve(lin_prob, MKLLUFactorization()).u
end

# # Old ResistanceProblem — Flagellate/CellBody/Flagellum API (Matrix-based discretisation)
#
# mutable struct ResistanceProblem{T<:Number} <: InstantaneousProblem
#     boundary::FluidBoundary
#     points::Discretisation
#     eps::T
#     mu::T
#     lin_prob::LinearProblem
#     force_vals::Union{Nothing,Vector{T}}
#     wall::Bool
# end
#
# function ResistanceProblem(boundary::FluidBoundary; eps=0.01, mu=1.0, wall=false)
#     @unpack N, Q, force_pts, quad_pts, velocity, nearest = boundary.points
#     points = NearestDiscretisation(N, Q, SVector(0.,0.,0.), I3,
#                                    zeros(T, size(force_pts)), zeros(T, size(force_pts)),
#                                    zeros(T, size(quad_pts)), nearest)
#     prob = ResistanceProblem(boundary, points, eps, mu,
#                              LinearProblem(zeros(T, 3N, 3N), zeros(T, 3N)), nothing, wall)
#     update_boundary!(prob, 0.0)
#     prob
# end
#
# get_velocities(prob::ResistanceProblem) = ...
# get_forces(prob::ResistanceProblem) = ...
# get_quad_pt_velocities(prob::ResistanceProblem; t=0.0) = ...
# update_boundary!(prob::ResistanceProblem, t) = ...
# add_rigid_body_motion!(prob::ResistanceProblem, U, Ω) = ...
# solve_problem!(prob::ResistanceProblem) = ...
# check_boundary_conditions(prob::ResistanceProblem; t=0.0) = ...

###########################################################################################
### Dynamic Problems ######################################################################
###########################################################################################

mutable struct SwimmingTrajectoryProblem <: DynamicProblem
    swimming_problem::SwimmingProblem
    ode_prob::ODEProblem
    traj::Union{Nothing, Trajectory}
end

function SwimmingTrajectoryProblem(
    ms::MicroSwimmer;
    x0=SVector(0.0, 0.0, 0.0),
    B=I3,
    t_final=20.0,
    saveat=0.05,
    eps=0.01,
    mu=1.0
)
    T = Float64
    sprob = SwimmingProblem(ms; eps=T(eps), mu=T(mu))
    @info "" typeof(sprob)

    x0_0 = SVector{3,T}(x0)
    b1_0 = SVector{3,T}(B[:,1])
    b2_0 = SVector{3,T}(B[:,2])
    X0   = vcat(x0_0, b1_0, b2_0)

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
    u   = sol.u
    x   = [SVector{3}(u[i][1:3]) for i in eachindex(u)]
    b1  = [SVector{3}(u[i][4:6]) for i in eachindex(u)]
    b2  = [SVector{3}(u[i][7:9]) for i in eachindex(u)]
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

# # ParticleTrajectoryProblem — tracks passive particles in the swimmer's flow field.
# # Needs updating to new ResistanceProblem interface before re-enabling.
#
# mutable struct ParticleTrajectoryProblem{T<:Number} <: Problem
#     resistance_problem::ResistanceProblem{T}
#     ode_prob::ODEProblem
#     t::Union{Nothing,Vector{T}}
#     trajectories::Union{Nothing,Matrix{T}}
# end
#
# function ParticleTrajectoryProblem(microswimmer; x=-5.0, ys=range(-4.,4.,6),
#                                    zs=range(0.2,3.2,6), t_final=20.0, saveat=0.05,
#                                    eps=0.01, mu=1.0)
#     ...
# end
#
# function solve_problem!(prob::ParticleTrajectoryProblem; method=Tsit5())
#     ...
# end


# only implemented for resistance problems currently
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

