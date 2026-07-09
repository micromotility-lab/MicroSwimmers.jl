struct Part{M <: Model, D <: Discretisation} <: FluidBoundary
    model::M
    disc::D
    frame::Frame{Float64}
end


function Part(model::Model, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3)
    part = Part(
        model,
        make_discretisation(model, N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64}(orientation))
    )
    init_boundary!(part, N, Q)
    nearest_neighbour!(part.disc)
    part
end

function add_rigid_body_motion!(part::Part, U, Ω)
    part.disc.velocity .= Ref(SVector{3}(U)) .+ cross.(Ref(SVector{3}(Ω)), part.disc.force_pts)
end

# fixed-cloud models: N,Q are final point counts → size the arrays
make_discretisation(::Model, N, Q)             = NearestDiscretisation(N, Q)
# raymarched models: N is a ray budget, cloud size unknown → start empty
make_discretisation(::ImplicitBodyModel, N, Q) = NearestDiscretisation()
    
init_boundary!(part::Part)       = init_boundary!(part.model, part.disc)
init_boundary!(part::Part, N, Q) = init_boundary!(part.model, part.disc)
init_boundary!(part::Part{<:ImplicitBodyModel}, N, Q) =  init_boundary!(part.model, part.disc, N, Q)
        
init_boundary!(m::FlagellumModel, disc)     = m(disc, 0.0)     # place at t=0
init_boundary!(m::CellBodyModel, disc)      = m(disc)          # fixed cloud
function init_boundary!(m::ImplicitBodyModel, disc, N, Q)
    m(disc, N, Q)       
    disc.nearest = zeros(Int, length(disc.quad_pts))
end

update_boundary!(part::Part, t::T) where {T <: Number} = update_boundary!(part.model, part.disc, t)
update_boundary!(::Model, disc::Discretisation, t::T) where {T <: Number} = nothing          # static default
update_boundary!(m::FlagellumModel, disc::Discretisation, t::T) where {T <: Number} = m(disc, t)       # deforming opts in

mutable struct MicroSwimmer{P <: Part} <: AbstractMicroSwimmer
    parts::Vector{P}
    frame::Frame{Float64}
end

MicroSwimmer(parts::Vector{P}) where {P <: Part} = MicroSwimmer(parts, Frame(zero(SVector{3}), I3))

update_boundary!(ms::MicroSwimmer, t::T) where {T <: Number} = foreach(p -> update_boundary!(p, t), ms.parts)

function grand_resistance_matrix(ms::MicroSwimmer; eps=0.1)
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez]) 
        prob = ResistanceProblem(ms, eps=eps)
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