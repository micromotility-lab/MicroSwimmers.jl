struct Part{M <: Model, D <: Discretisation} <: FluidBoundary
    model::M
    disc::D
    frame::Frame{Float64}
end

# function discretise(model::Model, disc::Discretisation, acc::Accessory..., location=zero(SVector{3,Float64}), orientation=I3)
#     part = Part(
#         model,
#         disc,
#         Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64,9}(orientation)),
#         acc
#     )
#     init_boundary!(part, )
# end

# `weighted=true` keeps this part's surface quadrature weights in the linear system, so its
# force unknowns become tractions rather than lumped forces. The default absorbs the weights
# into the unknowns, which is what every part did before this option existed.
function check_weighted(model::Model, weighted::Bool)
    weighted && !supports_quadrature_weights(model) && throw(ArgumentError(
        "quadrature weights are not implemented for $(typeof(model)); omit `weighted=true` " *
        "to absorb them into the force unknowns"))
    weighted
end

function discretise(model::Model, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3,
                    weighted=false, eps=DEFAULT_EPS)
    Part(model, N, Q; location=location, orientation=orientation, weighted=weighted, eps=eps)
end

function Part(model::Model, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3,
              weighted=false, eps=DEFAULT_EPS)
    check_weighted(model, weighted)
    check_tube(model, N, eps)
    part = Part(
        model,
        make_discretisation(model, N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64,9}(orientation))
    )
    init_boundary!(part, N, Q; weighted=weighted)
    nearest_neighbour!(part.disc)
    # after init_boundary!: sampled models only settle their region partition once the cloud
    # has been sampled, and set_eps! needs one value per region
    set_eps!(part.disc, eps)
    part
end

# Model-only form: size the discretisation from physical targets rather than making the caller
# invent N and Q. See src/defaults.jl for where the numbers come from.
#
# Two separate constraints set hq, and the binding one depends on eps. The quadrature has to
# resolve the blob (hq ≲ 2*eps), and it has to put several points in each force patch, or the
# piecewise-constant traction is not integrated at all — hence the hf/2 ceiling, which is what
# keeps Q ≥ 4N once eps grows past hf/4.
default_hq(hf, eps) = min(DEFAULT_HQ_FACTOR*minimum(eps), hf/2)

function Part(model::Model; hf=DEFAULT_HF, eps=DEFAULT_EPS, hq=default_hq(hf, eps), kwargs...)
    N = npoints_for_spacing(model, hf, false)
    Q = npoints_for_spacing(model, hq, true)
    check_default_sizing(model, N, Q, hf, hq, eps)
    Part(model, N, Q; eps=eps, kwargs...)
end

discretise(model::Model; kwargs...) = Part(model; kwargs...)

# The defaults are in microns, so they are wrong by orders of magnitude for a non-dimensional
# model — which is what most of the test suite uses. Say so rather than silently returning a
# two-point flagellum or a cloud that takes minutes to relax.
function check_default_sizing(model, N, Q, hf, hq, eps)
    N < 8 && @warn "Part($(nameof(typeof(model)))) sized to only N=$N force points at hf=$hf. " *
                   "The defaults assume microns; pass `hf` explicitly for a non-dimensional model." maxlog=1
    Q > 200_000 && @warn "Part($(nameof(typeof(model)))) needs Q=$Q quadrature points at hq=$hq. " *
                         "Raising `eps` above $eps would cut this quadratically." maxlog=1
    nothing
end

# Force points use the open rule (hf = L/(N+1)) and quadrature points the closed one
# (hq = L/(Q-1)) — see get_s0_and_ds in flagellum_models.jl.
npoints_for_spacing(m::FlagellumModel, h, include_endpoints) =
    include_endpoints ? max(round(Int, arclength(m)/h) + 1, 3) :
                        max(round(Int, arclength(m)/h) - 1, 2)

# Everything else is a surface. Its samplers spread N points over the whole area, so
# h = sqrt(area/N) — the same convention relax_cloud! uses to pick its target spacing. A model
# with no surface_area gets a clear error from that function's fallback.
npoints_for_spacing(m::Model, h, _) = max(ceil(Int, surface_area(m) / h^2), 8)

set_eps!(p::Part, ε) = (set_eps!(p.disc, ε); p)
get_eps(p::Part)     = get_eps(p.disc)

function add_rigid_body_motion!(part::Part, U, Ω)
    part.disc.velocity .= Ref(SVector{3}(U)) .+ cross.(Ref(SVector{3}(Ω)), part.disc.force_pts)
end

# models whose surface is sampled then filtered (raymarched, or carved by a groove): the
# requested N, Q are only sampling budgets, so the final cloud size is not known upfront
const SampledBodyModel = Union{ImplicitBodyModel, EllipsoidalGroovedBody}

# fixed-cloud models: N,Q are final point counts → size the arrays
make_discretisation(::Model, N, Q)             = NearestDiscretisation(N, Q)
# sampled models: cloud size unknown → start empty and size it from the sampled cloud
make_discretisation(::SampledBodyModel, N, Q)  = NearestDiscretisation()
function make_discretisation(m::PlanarVanedFlagellum, N, Q)
    N_v = vane_npoints(m.vane, N, m.flagellum.L, false)
    Q_v = vane_npoints(m.vane, Q, m.flagellum.L, true)
    NearestDiscretisation([N, N_v], [Q, Q_v])
end
# Tubes: N and Q are *station* counts along the centreline, which is exactly what the inherited
# arclength rule in npoints_for_spacing produces, so `Part(model)` sizes them correctly with no
# extra method. The cross-section counts come from the model's own fields.
make_discretisation(m::LineTubeFlagellum, N, Q)    = NearestDiscretisation(N, Q*m.Q_cs)
make_discretisation(m::SurfaceTubeFlagellum, N, Q) = NearestDiscretisation(N*m.N_cs, Q*m.Q_cs)


init_boundary!(part::Part)       = init_boundary!(part.model, part.disc)
init_boundary!(part::Part, N, Q; weighted=false) = init_boundary!(part.model, part.disc)
init_boundary!(part::Part{<:SampledBodyModel}, N, Q; weighted=false) =
    init_boundary!(part.model, part.disc, N, Q; weighted=weighted)
        
init_boundary!(m::FlagellumModel, disc)     = m(disc, 0.0)     
init_boundary!(m::CellBodyModel, disc)      = m(disc)         
# init_boundary!(m::PlanarVanedFlagellum, disc) = nothing


function init_boundary!(m::SampledBodyModel, disc, N, Q; weighted=false)
    m(disc, N, Q; weighted=weighted)
    disc.nearest            = zeros(Int, length(disc.quad_pts))
    disc.force_part_ranges  = [1:length(disc.force_pts)]
    disc.quad_part_ranges   = [1:length(disc.quad_pts)]
end

function update_boundary!(part::Part, t::T) where {T <: Number}
    update_boundary!(part.model, part.disc, t)
    # only weighted parts can have stale weights, and only deforming models define the hook
    is_weighted(part.disc) && quadrature_weights!(part.model, part.disc)
    nothing
end
update_boundary!(::Model, disc::Discretisation, t::T) where {T <: Number} = nothing          # static default
update_boundary!(m::FlagellumModel, disc::Discretisation, t::T) where {T <: Number} = m(disc, t)       # deforming opts in

mutable struct MicroSwimmer{P <: Part} <: AbstractMicroSwimmer
    parts::Vector{P}
    frame::Frame{Float64}
end

MicroSwimmer(parts::Vector{P}; location=zero(SVector{3}), orientation=I3) where {P <: Part} = MicroSwimmer(parts, Frame(location, orientation))

update_boundary!(ms::MicroSwimmer, t::T) where {T <: Number} = foreach(p -> update_boundary!(p, t), ms.parts)

add_rigid_body_motion!(ms::MicroSwimmer, U, Ω) = foreach(p -> add_rigid_body_motion!(p, U, Ω), ms.parts)

function grand_resistance_matrix(ms::MicroSwimmer; eps=nothing, alg=LUFactorization())
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez])
        prob = ResistanceProblem(ms, eps=eps, alg=alg)
        add_rigid_body_motion!(prob.microswimmer, n, zero(SVector{3,Float64}))
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, i] .= F
        R[4:6, i] .= T
        
        add_rigid_body_motion!(prob.microswimmer, zero(SVector{3,Float64}), n)
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, 3+i] .= F
        R[4:6, 3+i] .= T
    end
    R
end