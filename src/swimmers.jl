abstract type Swimmer <: FluidBoundary end

function move_boundary!(S::Swimmer, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::T) where {T <: Number}
    update_boundary!(S, t)
    S.points.location = x0
    S.points.orientation = B
end

struct Flagellum{M <: FlagellumModel} <: Swimmer   # subtypes Swimmer for isolated flagella
    model::M
    points::Discretisation
end

function Flagellum( 
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, -2π, 2π),
    N=23, 
    Q=127; 
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = NearestDiscretisation(
        zeros(3,N), zeros(3,Q), zeros(Int, Q); 
        location=location, orientation=orientation
    )
    f = Flagellum(model, points)

    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end

function TubeFlagellum( 
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, -2π, 2π),
    N=23, 
    N_cs=5,
    Q=127,
    Q_cs=5; 
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = TubeFlagellumNearestDiscretisation(N, N_cs, Q, Q_cs,
        zeros(3,N*N_cs), zeros(3,Q*Q_cs); 
        location=location, orientation=orientation
    )
    f = Flagellum(model, points)

    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end

function update_boundary!(f::Flagellum, t::T) where {T <: Number}
    f.model(f.points, t)
    # @unpack force_pts, quad_pts, velocity = f.points
    # f.model(force_pts, velocity, t)
    # f.model(quad_pts, t)
end


struct Flagellate{F <: Flagellum} <: Swimmer
    body::CellBody
    flagella::Vector{F}
    points::Discretisation
    force_pt_indices::Vector{Int} # locate the force points for each flagellum in the points structure
    quad_pt_indices::Vector{Int}  # locate the quad points 
end

function Flagellate(
    body::CellBody, 
    flagella::Vector{F};
    location=SVector(0., 0., 0.),
    orientation=I3,
) where {F <: Flagellum}
    force_pt_indices = [body.points.N + 1]
    quad_pt_indices  = [body.points.Q + 1]
    for i in eachindex(flagella[1:end-1])
        push!(force_pt_indices, body.points.N + 1 + sum(flagella[j].points.N - 1 for j = 1:i))
        push!(quad_pt_indices, body.points.Q + 1 + sum(flagella[j].points.Q - 1 for j in 1:i))
    end

    flagella_nearest = []
    for (i, flagellum) in enumerate(flagella)
        if i == 1
            append!(flagella_nearest, flagellum.points.nearest[2:end] .- 1)
        else
            append!(flagella_nearest, flagellum.points.nearest[2:end] .- 1 .+ sum(f.N - 1 for f in flagella[1:i-1]))
        end
    end

    nearest = [body.points.nearest; body.points.N .+ flagella_nearest]
    
    points = NearestDiscretisation(
        zeros(3, body.points.N + sum(f.points.N - 1 for f in flagella)),
        zeros(3, body.points.Q + sum(f.points.Q - 1 for f in flagella)),
        nearest;
        location=location,
        orientation=orientation
    )

    flgt = Flagellate(body, flagella, points, force_pt_indices, quad_pt_indices)

    @views begin
        flgt.points.force_pts[:, 1:body.points.N] .= body.points.force_pts
        flgt.points.quad_pts[:,  1:body.points.Q] .= body.points.quad_pts
    end
    update_boundary!(flgt, 0.)
    flgt
end

UniFlagellate(
    body::CellBody,
    flagellum::Flagellum;
    location=SVector(0., 0., 0.),
    orientation=I3,
) = Flagellate(body, [flagellum]; location=location, orientation=orientation)


function update_boundary!(flgt::Flagellate, t::T) where {T <: Number}
    for (i, flagellum) in enumerate(flgt.flagella) 
        update_boundary!(flagellum, t)
    
        f_start = flgt.force_pt_indices[i]
        q_start = flgt.quad_pt_indices[i]
        
        @unpack force_pts, quad_pts, velocity, location, orientation = flagellum.points
        @views begin
            flgt.points.force_pts[:, f_start:f_start + flagellum.points.N-2] .= location .+ orientation * force_pts[:,2:end]
            flgt.points.quad_pts[:,  q_start:q_start + flagellum.points.Q-2] .= location .+ orientation * quad_pts[:,2:end]
            flgt.points.velocity[:,  f_start:f_start + flagellum.points.N-2] .= orientation * velocity[:,2:end]
        end 
    end
end
