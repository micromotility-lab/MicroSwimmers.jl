mutable struct CellBody <: MicroSwimmer
    model::CellBodyModel
    points::NearestDiscretisation
end

function CellBody(
    model::CellBodyModel,
    N::Int,
    Q::Int; 
    location=SVector(0., 0., 0.),
    orientation=I3
)
    force_pts = model(N)
    quad_pts  = model(Q)

    # remove empty parts
    force_pts_clean = [force_pts[i] for i in eachindex(force_pts) if size(force_pts[i], 2) > 0]
    quad_pts_clean = [quad_pts[i] for i in eachindex(force_pts) if size(force_pts[i], 2) > 0]

    nearest = nearest_neighbour(force_pts_clean, quad_pts_clean)
    points = NearestDiscretisation(
        reduce(hcat, force_pts_clean), 
        reduce(hcat, quad_pts_clean), 
        nearest; 
        location=SVector{3}(location), 
        orientation=SMatrix{3,3}(orientation)
    )
    CellBody(model, points)
end

SphericalBody(a=0.2; N=27, Q=99) = CellBody(EllipsoidBody(a, a, a), N, Q)


function grand_resistance_matrix(cell::CellBody; eps=0.1)
    prob = ResistanceProblem(cell, eps=eps)
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez]) 
        prob = ResistanceProblem(cell, eps=eps)
        add_rigid_body_motion!(prob, n, 0.0)
        F, T = total_force_and_torque(prob)
        R[1:3, i] .= F
        R[4:6, i] .= T
        
        add_rigid_body_motion!(prob, 0.0, n)
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, 3+i] .= F
        R[4:6, 3+i] .= T
    end
    R
end


