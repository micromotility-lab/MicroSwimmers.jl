abstract type Problem end

abstract type InstantaneousProblem{MS <: MicroSwimmer, D <: Discretisation} <: Problem end
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

function translate_problem!(prob::InstantaneousProblem, x0::AbstractVector{T}) where {T <: Number}
    prob.disc.force_pts .= Ref(SVector{3,T}(x0)) .+ prob.disc.force_pts
    prob.disc.quad_pts  .= Ref(SVector{3,T}(x0)) .+ prob.disc.quad_pts
    fr = prob.microswimmer.frame
    prob.microswimmer.frame = Frame(fr.location + SVector{3,T}(x0), fr.orientation)
end

function rotate_problem!(prob::InstantaneousProblem, B::AbstractMatrix{T}) where {T <: Number}
    prob.disc.force_pts .= Ref(SMatrix{3,3,T,9}(B)) .* prob.disc.force_pts
    prob.disc.velocity  .= Ref(SMatrix{3,3,T,9}(B)) .* prob.disc.velocity
    prob.disc.quad_pts  .= Ref(SMatrix{3,3,T,9}(B)) .* prob.disc.quad_pts
    fr = prob.microswimmer.frame
    prob.microswimmer.frame = Frame(fr.location, SMatrix{3,3,T,9}(B) * fr.orientation)
end

function gather_nearest!(prob::InstantaneousProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc = prob
    for i in eachindex(microswimmer.parts)
        part  = microswimmer.parts[i]
        f_rng = disc.force_part_ranges[i]
        q_rng = disc.quad_part_ranges[i]
        foff  = first(f_rng) - 1
        @views disc.nearest[q_rng] .= part.disc.nearest .+ foff
    end
end

function gather!(prob::InstantaneousProblem{<:Any, <:NearestDiscretisation})
    @unpack microswimmer, disc = prob
    for (i, p) in enumerate(microswimmer.parts)
        gather_part!(disc, microswimmer.frame, p, i)   # <- barrier
    end
end
function gather_part!(disc, ms_frame::Frame, p::Part, i)
    lab_frame = ms_frame * p.frame
    f_rng     = disc.force_part_ranges[i]
    q_rng     = disc.quad_part_ranges[i]
    @views disc.force_pts[f_rng] .= lab_frame.(p.disc.force_pts)
    @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* p.disc.velocity
    @views disc.quad_pts[q_rng]  .= lab_frame.(p.disc.quad_pts)
    gather_part_weights!(disc, p, q_rng)
end

# A Frame is a rigid motion, so it preserves area: weights copy across untransformed. Parts
# with no weights of their own contribute w == 1, which is exactly the absorbed convention,
# so a weighted and an unweighted part can share one linear system. Re-copying on every
# gather (an O(Q) memcpy against an O(N*Q) assembly) keeps the part's own weights
# authoritative for models that update them.
function gather_part_weights!(disc::NearestDiscretisation, p::Part, q_rng)
    isnothing(disc.quad_wts) && return nothing
    if isnothing(p.disc.quad_wts)
        @views disc.quad_wts[q_rng] .= one(eltype(disc.quad_wts))
    else
        @views disc.quad_wts[q_rng] .= p.disc.quad_wts
    end
    nothing
end

# Give the global discretisation a weight vector only if at least one part supplies weights;
# otherwise leave it `nothing` so the whole solve takes the absorbed (zero-cost) path.
function init_quad_weights!(disc::NearestDiscretisation, ms::MicroSwimmer)
    any(p -> is_weighted(p.disc), ms.parts) || return disc
    disc.quad_wts = ones(eltype(eltype(disc.quad_pts)), nq(disc))
    disc
end

# function gather!(prob::InstantaneousProblem{<:Any, <:NearestDiscretisation})
#     @unpack microswimmer, disc = prob
#     for i in eachindex(microswimmer.parts)
#         part      = microswimmer.parts[i]
#         lab_frame = microswimmer.frame * part.frame
#         f_rng     = disc.force_part_ranges[i]
#         q_rng     = disc.quad_part_ranges[i]
#         @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
#         @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
#         @views disc.quad_pts[q_rng]  .= lab_frame.(part.disc.quad_pts)
#     end
# end

function gather!(prob::InstantaneousProblem{<:Any, <:NystromDiscretisation})
    @unpack microswimmer, disc = prob
    fstart = 1
    for part in microswimmer.parts
        lab_frame = microswimmer.frame * part.frame
        nf_i      = nf(part.disc)
        f_rng     = fstart:fstart+nf_i-1
        @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
        @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
        fstart += nf_i
    end
end

###########################################################################################
### Hybrid dense-LU / GMRES solve support #################################################
###########################################################################################
#
# For SwimmingTrajectoryProblem, consecutive timesteps perturb the geometry only slightly, so
# the previous step's dense LU factorisation is nearly as good a solver for the new step's
# matrix. swimming_function_operator wraps mul_swimming! (numerics.jl) as the matrix-free A
# used by GMRES on those in-between steps; PreconditionerBox holds the last dense
# factorisation to use as a left preconditioner, in a wrapper whose type stays fixed across
# repeated dense refactors (LinearSolve.LinearCache.Pl's type is fixed at `init` time, but
# each dense solve produces a new LU object).
#
# LinearProblem requires A to either be an AbstractArray or an AbstractSciMLOperator (it
# introspects call signatures on anything else, expecting an ODE-style `f(u,p,t)` function —
# a plain matrix-free struct doesn't satisfy that and errors). Subtyping AbstractSciMLOperator
# is what SciMLOperators.FunctionOperator does internally too, but FunctionOperator's own
# p/u/state-threading machinery (in the installed SciMLOperators version) does not reliably
# propagate a live mutable state object through to mul! — only size/eltype/mul! are actually
# required for SimpleGMRES, so implementing AbstractSciMLOperator directly avoids that
# indirection entirely: op.disc below is a plain field read on an object we construct and
# control, not something routed through a third party's internal calling convention.
struct SwimmingOperator{MS <: MicroSwimmer, D <: Discretisation, K <: Kernel, T <: Number} <:
       SciMLOperators.AbstractSciMLOperator{T}
    ms::MS
    disc::D
    kernel::K
    mu::T
    N::Int   # nf(disc); operator size is (3N+6, 3N+6)
end

Base.size(op::SwimmingOperator) = (3*op.N+6, 3*op.N+6)
Base.size(op::SwimmingOperator, i::Int) = size(op)[i]
Base.eltype(::SwimmingOperator{<:Any, <:Any, <:Any, T}) where {T} = T

LinearAlgebra.mul!(y::AbstractVector, op::SwimmingOperator, x::AbstractVector) =
    mul_swimming!(y, x, op.ms.frame.location, op.disc, op.kernel; μ=op.mu)

swimming_function_operator(ms::MicroSwimmer, disc::Discretisation, kernel::Kernel, mu, N) =
    SwimmingOperator(ms, disc, kernel, Float64(mu), N)

mutable struct PreconditionerBox{F}
    lu::F
end

# SimpleGMRES applies Pl via mul! (Pl stands in for the inverse operator M⁻¹, not M),
# so "multiplying" by this box means solving with the stored factorisation.
LinearAlgebra.mul!(y::AbstractVector, P::PreconditionerBox, x::AbstractVector) = ldiv!(y, P.lu, x)
LinearAlgebra.ldiv!(y::AbstractVector, P::PreconditionerBox, x::AbstractVector) = ldiv!(y, P.lu, x)
LinearAlgebra.ldiv!(P::PreconditionerBox, x::AbstractVector) = ldiv!(P.lu, x)

###########################################################################################
### SwimmingProblem #######################################################################
###########################################################################################

mutable struct SwimmingProblem{MS <: MicroSwimmer, D <: Discretisation, T <: Number, K <: Kernel, C <: LinearSolve.LinearCache, G <: Union{Nothing, LinearSolve.LinearCache}, P <: Union{Nothing, PreconditionerBox}} <: InstantaneousProblem{MS, D}
    microswimmer::MS
    disc::D
    mu::T
    cache::C
    force_vals::Union{Nothing, Vector{T}}
    kernel::K
    hybrid::Bool
    gmres_cache::G
    pl_box::P
    steps_since_refactor::Int
    refactor_interval::Int
    max_refactor_interval::Int
    gmres_maxiters::Int
end

function make_hybrid_cache(ms::MicroSwimmer, disc::Discretisation, kernel::Kernel, mu, N;
    gmres_reltol, gmres_maxiters
)
    op  = swimming_function_operator(ms, disc, kernel, Float64(mu), N)
    box = PreconditionerBox(lu(Matrix{Float64}(I, 3N+6, 3N+6)))
    # SimpleGMRES's own default_alias_A is false (unlike other Krylov algorithms), so
    # `init` would otherwise deepcopy `op` — freezing its captured `ms`/`disc` at their
    # current (pre-gather_nearest!) state instead of tracking the live problem. Force
    # aliasing so the operator keeps seeing the same, live, mutated-in-place objects.
    gmres_cache = init(
        LinearProblem(op, zeros(3N+6)), SimpleGMRES();
        alias=LinearSolve.SciMLBase.LinearAliasSpecifier(alias_A=true),
        Pl=box, u0=zeros(3N+6), reltol=gmres_reltol, maxiters=gmres_maxiters
    )
    gmres_cache, box
end

function SwimmingProblem(ms::MicroSwimmer{<:Part{<:Model, <:NearestDiscretisation}};
    mu=1.0, eps=0.1, wall=false, alg=LUFactorization(),
    hybrid=false, refactor_interval=20, max_refactor_interval=100,
    gmres_reltol=1e-10, gmres_maxiters=50
)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    nq_sizes = [nq(p.disc) for p in ms.parts]
    N = sum(nf_sizes); Q = sum(nq_sizes)
    disc = NearestDiscretisation(nf_sizes, nq_sizes)
    init_quad_weights!(disc, ms)
    kernel = wall ? RegBlakelet(eps) : RegStokeslet(eps)

    gmres_cache, pl_box = hybrid ?
        make_hybrid_cache(ms, disc, kernel, mu, N; gmres_reltol=gmres_reltol, gmres_maxiters=gmres_maxiters) :
        (nothing, nothing)

    prob = SwimmingProblem(
        ms,
        disc,
        Float64(mu),
        init(LinearProblem(zeros(3N+6, 3N+6), zeros(3N+6)), alg),
        nothing,
        kernel,
        hybrid,
        gmres_cache,
        pl_box,
        0,
        refactor_interval,
        max_refactor_interval,
        gmres_maxiters
    )
    gather_nearest!(prob)
    update_boundary!(prob, 0.0)
    prob
end

function SwimmingProblem(ms::MicroSwimmer{<:Part{<:Model, <:NystromDiscretisation}};
    mu=1.0, eps=0.1, alg=LUFactorization(),
    hybrid=false, refactor_interval=20, max_refactor_interval=100,
    gmres_reltol=1e-7, gmres_maxiters=5
)
    nf_sizes = [nf(p.disc) for p in ms.parts]
    N        = sum(nf_sizes)
    disc     = NystromDiscretisation(N)
    kernel   = RegStokeslet(eps)

    # make_hybrid_cache already returns a fully initialised PreconditionerBox, so the
    # non-hybrid branch has nothing to preallocate — it just carries `nothing`, exactly as
    # the NearestDiscretisation constructor above does.
    gmres_cache, pl_box = hybrid ?
        make_hybrid_cache(ms, disc, kernel, mu, N; gmres_reltol=gmres_reltol, gmres_maxiters=gmres_maxiters) :
        (nothing, nothing)

    prob = SwimmingProblem(
        ms,
        disc,
        Float64(mu),
        init(LinearProblem(zeros(3N+6, 3N+6), zeros(3N+6)), alg),
        nothing,
        kernel,
        hybrid,
        gmres_cache,
        pl_box,
        0,
        refactor_interval,
        max_refactor_interval,
        gmres_maxiters
    )
    update_boundary!(prob, 0.0)
    prob
end

# function gather_nearest!(prob::SwimmingProblem{<:Any, <:NearestDiscretisation})
#     @unpack microswimmer, disc = prob
#     for i in eachindex(microswimmer.parts)
#         part  = microswimmer.parts[i]
#         f_rng = disc.force_part_ranges[i]
#         q_rng = disc.quad_part_ranges[i]
#         foff  = first(f_rng) - 1
#         @views disc.nearest[q_rng] .= part.disc.nearest .+ foff
#     end
# end

# function gather!(prob::SwimmingProblem{<:Any, <:NearestDiscretisation})
#     @unpack microswimmer, disc = prob
#     for i in eachindex(microswimmer.parts)
#         part      = microswimmer.parts[i]
#         lab_frame = microswimmer.frame * part.frame
#         f_rng     = disc.force_part_ranges[i]
#         q_rng     = disc.quad_part_ranges[i]
#         @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
#         @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
#         @views disc.quad_pts[q_rng]  .= lab_frame.(part.disc.quad_pts)
#     end
# end

# function gather!(prob::SwimmingProblem{<:Any, <:NystromDiscretisation})
#     @unpack microswimmer, disc = prob
#     fstart = 1
#     for part in microswimmer.parts
#         lab_frame = microswimmer.frame * part.frame
#         nf_i      = nf(part.disc)
#         f_rng     = fstart:fstart+nf_i-1
#         @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
#         @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
#         fstart += nf_i
#     end
# end

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
    # gather!(prob)
end

function move_boundary!(prob::SwimmingProblem, x0::AbstractVector, B::AbstractMatrix, t::T) where {T <: Number} 
    prob.microswimmer.frame = Frame(SVector{3,T}(x0), SMatrix{3,3,T,9}(B))
    update_boundary!(prob, T(t))
end

function move_boundary!(prob::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3,T}, t::Number) where T
    move_boundary!(prob, x0, SMatrix{3,3,T,9}(hcat(b1, b2, cross(b1, b2))), t)
end

function fill_swimming_rhs!(b, disc)
    N3 = 3 * nf(disc)
    T  = eltype(eltype(disc.velocity))
    @views b[1:N3] .= reinterpret(T, disc.velocity)
    @views b[N3+1:end] .= zero(T)
end

function dense_solve!(prob::SwimmingProblem)
    @unpack cache, disc, kernel, mu, microswimmer, pl_box = prob
    assemble_swimming!(cache.A, microswimmer.frame.location, disc, kernel; μ=mu)
    cache.isfresh = true  # A was mutated in-place; tell LinearSolve to refactorise
    fill_swimming_rhs!(cache.b, disc)
    prob.force_vals = solve!(cache).u
    if prob.hybrid
        copyto!(pl_box.lu.factors, cache.cacheval.factors)
        copyto!(pl_box.lu.ipiv,    cache.cacheval.ipiv)
        pl_box.lu = LU(pl_box.lu.factors, pl_box.lu.ipiv, cache.cacheval.info)
        prob.steps_since_refactor = 0
    end

    prob
end

# Attempts a preconditioned-GMRES solve reusing the last dense factorisation as `Pl`.
# Returns true and updates prob.force_vals on convergence; returns false (leaving
# prob.force_vals untouched) if GMRES fails to converge, so the caller can fall back to a
# dense solve rather than accepting an unconverged answer.
function try_gmres_solve!(prob::SwimmingProblem)
    gc = prob.gmres_cache
    fill_swimming_rhs!(gc.b, prob.disc)
    # SimpleGMRES's internal cache snapshots A/b/x once and does not pick up new values on
    # repeat solve! calls unless isfresh is set — force it to rebuild from the current b
    # (preconditioned by the still-stale pl_box.lu) every hybrid step.
    gc.isfresh = true
    sol = solve!(gc)

    if sol.retcode != LinearSolve.ReturnCode.Success
        return false
    end
    prob.force_vals = sol.u
    prob.steps_since_refactor += 1

    # Adaptive refactor interval: the preconditioner is drifting stale if GMRES is working
    # hard to converge, so refactorise sooner; if it's converging easily, wait longer.
    # Symmetric ±1 (rather than halving down / incrementing up) so a handful of rough steps
    # don't collapse the interval to 1 and get stuck oscillating dense/GMRES/dense/GMRES —
    # which is worse than just staying dense, since a failed or maxed-out GMRES attempt gets
    # paid for on top of the dense solve that follows it.
    if sol.iters > 0.7 * prob.gmres_maxiters
        prob.refactor_interval = max(1, prob.refactor_interval - 1)
    elseif sol.iters < 0.3 * prob.gmres_maxiters
        prob.refactor_interval = min(prob.max_refactor_interval, prob.refactor_interval + 1)
    end
    true
end

function solve_problem!(prob::SwimmingProblem)
    gather!(prob)
    # @info "" prob.steps_since_refactor
    needs_dense = !prob.hybrid || isnothing(prob.force_vals) ||
                  prob.steps_since_refactor >= prob.refactor_interval

    if needs_dense
        dense_solve!(prob)
    elseif !try_gmres_solve!(prob)
        # Preconditioner too stale to converge — fall back to a fresh dense factorisation
        # rather than return an unconverged answer.
        dense_solve!(prob)
    end
    prob.force_vals
end


###########################################################################################
### ResistanceProblem #####################################################################
###########################################################################################

mutable struct ResistanceProblem{MS <: MicroSwimmer, D <: Discretisation, T <: Number, K <: Kernel, C <: LinearSolve.LinearCache} <: InstantaneousProblem{MS, D}
    microswimmer::MS
    disc::D
    mu::T
    cache::C
    force_vals::Union{Nothing, Vector{T}}
    kernel::K
end

function ResistanceProblem(ms::MicroSwimmer{<:Part{<:Model, <:NearestDiscretisation}}; mu=1.0, eps=0.1, wall=false, alg=LUFactorization())
    nf_sizes = [nf(p.disc) for p in ms.parts]
    nq_sizes = [nq(p.disc) for p in ms.parts]
    N = sum(nf_sizes); Q = sum(nq_sizes)
    disc = NearestDiscretisation(nf_sizes, nq_sizes)
    init_quad_weights!(disc, ms)
    prob = ResistanceProblem(
        ms,
        disc,
        Float64(mu),
        init(LinearProblem(zeros(3N, 3N), zeros(3N)), alg),
        nothing,
        wall ? RegBlakelet(eps) : RegStokeslet(eps)
    )
    gather_nearest!(prob)
    prob
end

# function ResistanceProblem(ms::MicroSwimmer{<:Part{<:Model, <:NystromDiscretisation}}; mu=1.0, eps=0.1)
#     nf_sizes = [nf(p.disc) for p in ms.parts]
#     N        = sum(nf_sizes)
#     ResistanceProblem(
#         ms,
#         NystromDiscretisation(N),
#         Float64(mu),
#         LinearProblem(zeros(3N, 3N), zeros(3N)),
#         nothing,
#         RegStokeslet(eps)
#     )
# end

update_boundary!(prob::ResistanceProblem, t::Number) = update_boundary!(prob.microswimmer, t)


# function gather_nearest!(prob::ResistanceProblem{<:Any, <:NearestDiscretisation})
#     @unpack microswimmer, disc = prob
#     for i in eachindex(microswimmer.parts)
#         part  = microswimmer.parts[i]
#         f_rng = disc.force_part_ranges[i]
#         q_rng = disc.quad_part_ranges[i]
#         foff  = first(f_rng) - 1
#         @views disc.nearest[q_rng] .= part.disc.nearest .+ foff
#     end
# end

# function gather!(prob::ResistanceProblem{<:Any, <:NearestDiscretisation})
#     @unpack microswimmer, disc = prob
#     for i in eachindex(microswimmer.parts)
#         part      = microswimmer.parts[i]
#         lab_frame = microswimmer.frame * part.frame
#         f_rng     = disc.force_part_ranges[i]
#         q_rng     = disc.quad_part_ranges[i]
#         @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
#         @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
#         @views disc.quad_pts[q_rng]  .= lab_frame.(part.disc.quad_pts)
#     end
# end

# function gather!(prob::ResistanceProblem{<:Any, <:NystromDiscretisation})
#     @unpack microswimmer, disc = prob
#     fstart = 1
#     for part in microswimmer.parts
#         lab_frame = microswimmer.frame * part.frame
#         nf_i      = nf(part.disc)
#         f_rng     = fstart:fstart+nf_i-1
#         @views disc.force_pts[f_rng] .= lab_frame.(part.disc.force_pts)
#         @views disc.velocity[f_rng]  .= Ref(lab_frame.orientation) .* part.disc.velocity
#         fstart += nf_i
#     end
# end

function get_forces(prob::ResistanceProblem)
    check_solved!(prob)
    fv = prob.force_vals
    N  = nf(prob.disc)
    [SVector{3}(fv[3i-2], fv[3i-1], fv[3i]) for i in 1:N]
end

function solve_problem!(prob::ResistanceProblem)
    @unpack cache, disc, kernel, mu = prob
    gather!(prob)
    assemble!(cache.A, disc, kernel; μ=mu)
    cache.isfresh = true  # A was mutated in-place; tell LinearSolve to refactorise
    cache.b .= reinterpret(eltype(eltype(disc.velocity)), disc.velocity)
    prob.force_vals = solve!(cache).u
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
    # x0=SVector(0.0, 0.0, 0.0),
    # B=I3,
    t_final=20.0,
    saveat=0.05,
    eps=0.1,
    mu=1.0,
    wall=false,
    alg=LUFactorization(),
    hybrid=false,
    refactor_interval=20,
    max_refactor_interval=100,
    gmres_reltol=1e-10,
    gmres_maxiters=50
)
    T = Float64
    sprob = SwimmingProblem(ms; eps=T(eps), mu=T(mu), wall=wall,
        hybrid=hybrid, alg=alg, refactor_interval=refactor_interval,
        max_refactor_interval=max_refactor_interval,
        gmres_reltol=gmres_reltol, gmres_maxiters=gmres_maxiters)

    X0 = SVector{9,T}(ms.frame.location..., ms.frame.orientation[:,1]..., ms.frame.orientation[:,2]...)
    # X0   = SVector{9,T}(x0..., B[:,1]..., B[:,2]...)

    function rhs(X, p, t)
        x0 = SVector(X[1],X[2],X[3])
        b1 = SVector(X[4],X[5],X[6])
        b2 = SVector(X[7],X[8],X[9])
        move_boundary!(sprob, x0, b1, b2, t)
        solve_problem!(sprob)
        Ω = get_Ω(sprob)
        SVector{9,T}(get_U(sprob)..., cross(Ω, b1)..., cross(Ω, b2)...)
    end

    SwimmingTrajectoryProblem(
        sprob,
        ODEProblem(rhs, X0, (T(0), T(t_final)), saveat=T(saveat)),
        nothing
    )
end

function solve_problem!(prob::SwimmingTrajectoryProblem; method=VCABM(), periodic=false)
    sol = solve(prob.ode_prob, method)
    u   = sol.u
    T = eltype(eltype(u))
    x   = [SVector{3,T}(u[i][1:3]) for i in eachindex(u)]
    b1  = [SVector{3,T}(u[i][4:6]) for i in eachindex(u)]
    b2  = [SVector{3,T}(u[i][7:9]) for i in eachindex(u)]
    prob.traj = Trajectory{T}(sol.t, x, b1, b2, periodic)
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
    N = nf(rprob.disc)
    A = zeros(3 * num_particles, 3N)

    function rhs!(dX, X, p, t)
        update_boundary!(rprob, t)
        solve_problem!(rprob)
        for i in 1:num_particles
            x_sv = SVector{3}(X[3i-2], X[3i-1], X[3i])
            Ai   = @view A[3i-2:3i, :]
            assemble!(Ai, [x_sv], rprob.disc.quad_pts, rprob.disc.nearest,
                      rprob.disc.quad_wts, rprob.kernel; μ=rprob.mu)
        end
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