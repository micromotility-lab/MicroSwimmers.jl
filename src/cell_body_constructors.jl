mutable struct CellBody <: AbstractMicroSwimmer
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
    quad_pts_clean = [quad_pts[i] for i in eachindex(quad_pts) if size(quad_pts[i], 2) > 0]

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