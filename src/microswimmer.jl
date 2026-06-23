struct Part{M <: Model, D <: Discretisation}
    model::M
    disc::D
    frame::Frame{Float64}
end

function Part(model::Model, N, Q; location=zero(SVector{3, Float64}), orientation=I3)
    part = Part(
        model,
        NearestDiscretisation(N, Q),
        Frame(SVector{3,Float64}(location), SMatrix{3,3,Float64}(orientation))
    )
    update_boundary!(part, 0.0)
    nearest_neighbour!(part.disc)
    part
end

update_boundary!(part::Part, t) = part.model(part.disc, t)

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