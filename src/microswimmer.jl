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

function discretise(model::Model, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3, weighted=false)
    check_weighted(model, weighted)
    part = Part(
        model,
        make_discretisation(model, N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64,9}(orientation))
    )
    init_boundary!(part, N, Q; weighted=weighted)
    nearest_neighbour!(part.disc)
    part
end

function Part(model::Model, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3, weighted=false)
    check_weighted(model, weighted)
    part = Part(
        model,
        make_discretisation(model, N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64,9}(orientation))
    )
    init_boundary!(part, N, Q; weighted=weighted)
    nearest_neighbour!(part.disc)
    part
end

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

add_rigid_body_motion!!(ms::MicroSwimmer, U, Ω) = foreach(p -> add_rigid_body_motion!(p, U, Ω), ms.parts)

function grand_resistance_matrix(ms::MicroSwimmer; eps=0.1, alg=LUFactorization())
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez])
        prob = ResistanceProblem(ms, eps=eps, alg=alg)
        [add_rigid_body_motion!(part, n, zero(SVector{3,Float64})) for part in prob.microswimmer.parts]
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, i] .= F
        R[4:6, i] .= T
        
        [add_rigid_body_motion!(part, zero(SVector{3,Float64}), n) for part in prob.microswimmer.parts]
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, 3+i] .= F
        R[4:6, 3+i] .= T
    end
    R
end