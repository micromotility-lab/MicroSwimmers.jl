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

# function CellBody(
#     model::Union{EllipsoidalGroovedBody, CylindricalGroovedBody, FlatGroovedBody},
#     N::Int, 
#     Q::Int;
#     location=SVector(0., 0., 0.),
#     orientation=I3
# )
#     body_force_pts, groove_force_pts = model(N)
#     body_quad_pts, groove_quad_pts   = model(Q)
    
#     force_pts = hcat(body_force_pts, groove_force_pts)
#     quad_pts =  hcat(body_quad_pts, groove_quad_pts)
#     nearest = nearest_neighbour(force_pts, quad_pts)
#     points = NearestDiscretisation(force_pts, quad_pts, nearest; location=SVector{3}(location), orientation=SMatrix{3,3}(orientation))
#     CellBody(model, points)
# end

SphericalBody(a=0.2; N=27, Q=99) = CellBody(EllipsoidBody(a, a, a), N, Q)


function RigidMotionBody(
    model::CellBodyModel,
    N::Int,
    Q::Int,
    U::AbstractVector,
    Ω::AbstractVector
)
    body = CellBody(model, N, Q)
    body.points.velocity .= model(body.points.force_pts, U, Ω)
    body
end

function grand_resistance_matrix(cell::CellBody; eps=0.1)
    prob = ResistanceProblem(cell, eps=eps)
    R = zeros(6,6)

    for (i, n) in enumerate([ex, ey, ez]) 
        reset_velocity!(cell)  
        add_velocity!(cell, n)
        prob = ResistanceProblem(cell, eps=eps)
        F, T = total_force_and_torque(prob)
        R[1:3, i] .= F
        R[4:6, i] .= T
        
        reset_velocity!(cell) 
        add_angular_velocity!(cell, n)
        prob = ResistanceProblem(cell, eps=eps)
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        R[1:3, 3+i] .= F
        R[4:6, 3+i] .= T
    end
    R
end


