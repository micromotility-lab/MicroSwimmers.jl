struct Flagellate{F <: Flagellum} <: MicroSwimmer
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
    T = eltype(flagella[1].points.force_pts)
    location = SVector{3,T}(location)
    orientation = SMatrix{3,3,T}(orientation)

    force_pt_indices = [body.points.N + 1]
    quad_pt_indices  = [body.points.Q + 1]
    for i in eachindex(flagella[1:end-1])
        push!(force_pt_indices, body.points.N + 1 + sum(flagella[j].points.N  for j = 1:i))
        push!(quad_pt_indices, body.points.Q + 1 + sum(flagella[j].points.Q for j in 1:i))
    end

    flagella_nearest = []
    for (i, flagellum) in enumerate(flagella)
        if i == 1
            append!(flagella_nearest, flagellum.points.nearest)
        else
            append!(flagella_nearest, flagellum.points.nearest .+ sum(f.points.N for f in flagella[1:i-1]))
        end
    end

    nearest = [body.points.nearest; body.points.N .+ flagella_nearest]
    
    points = NearestDiscretisation(
        zeros(T, 3, body.points.N + sum(f.points.N for f in flagella)),
        zeros(T, 3, body.points.Q + sum(f.points.Q for f in flagella)),
        nearest;
        location=location,
        orientation=orientation
    )

    flgt = Flagellate(body, flagella, points, force_pt_indices, quad_pt_indices)

    @views begin
        flgt.points.force_pts[:, 1:body.points.N] .= body.points.force_pts
        flgt.points.quad_pts[:,  1:body.points.Q] .= body.points.quad_pts
    end

    update_boundary!(flgt, zero(T))
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
            flgt.points.force_pts[:, f_start:f_start + flagellum.points.N-1] .= location .+ orientation * force_pts
            flgt.points.quad_pts[:,  q_start:q_start + flagellum.points.Q-1] .= location .+ orientation * quad_pts
            flgt.points.velocity[:,  f_start:f_start + flagellum.points.N-1] .= orientation * velocity
        end 
    end
end
