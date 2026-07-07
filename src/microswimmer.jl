struct Part{M <: Model, D <: Discretisation} <: FluidBoundary
    model::M
    disc::D
    frame::Frame{Float64}
end


function Part(model::Model, N::Int, Q::Int; location=zero(SVector{3, Float64}), orientation=I3)
    part = Part(
        model,
        NearestDiscretisation(N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64}(orientation))
    )
    init_boundary!(part, N, Q)
    nearest_neighbour!(part.disc)
    part
end

# function Part(model::ImplicitBodyModel, N::Int, Q::Int; location=zero(SVector{3,Float64}), orientation=I3)
#     part = Part(
    #         model,
    #         NearestDiscretisation(),
    #         Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64}(orientation))
    #     )
    #     init_boundary!(model, part.disc, N, Q)
    #     part.disc.nearest = zeros(Int, length(part.disc.quad_pts))
    #     nearest_neighbour!(part.disc)
    #     part
    # end
    
init_boundary!(part::Part)       = init_boundary!(part.model, part.disc)
init_boundary!(part::Part, N, Q) = init_boundary!(part.model)
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
    prob = ResistanceProblem(cell, eps=eps)
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez]) 
        prob = ResistanceProblem(cell, eps=eps)
        add_rigid_body_motion!(prob, n, [0.0, 0.0, 0.0])
        F, T = total_force_and_torque(prob)
        R[1:3, i] .= F
        R[4:6, i] .= T
        
        add_rigid_body_motion!(prob, [0.0, 0.0, 0.0], n)
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, 3+i] .= F
        R[4:6, 3+i] .= T
    end
    R
end