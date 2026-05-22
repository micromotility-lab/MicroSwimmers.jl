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

function discretisation(flgt::Flagellate)
    @info "Discretisation" flgt.body.points.N flgt.body.points.Q flgt.flagella[1].points.N flgt.flagella[1].points.Q
end

struct Colony{F <: Flagellate} <: MicroSwimmer
    members::Vector{F}
    points::Discretisation
    member_force_indices::Vector{Int}  # start index of each member's force pts
    member_quad_indices::Vector{Int}   # start index of each member's quad pts
end

function Colony(members::Vector{F};
    location=SVector(0., 0., 0.),
    orientation=I3,
) where {F <: Flagellate}
    T = eltype(members[1].points.force_pts)
    location = SVector{3,T}(location)
    orientation = SMatrix{3,3,T}(orientation)

    # Build index offsets exactly as Flagellate does for its flagella
    member_force_indices = [1]
    member_quad_indices  = [1]
    for i in eachindex(members[1:end-1])
        push!(member_force_indices, 1 + sum(m.points.N for m in members[1:i]))
        push!(member_quad_indices,  1 + sum(m.points.Q for m in members[1:i]))
    end

    # Concatenate nearest neighbour indices with offsets
    colony_nearest = Int[]
    for (i, member) in enumerate(members)
        offset = i == 1 ? 0 : sum(m.points.N for m in members[1:i-1])
        append!(colony_nearest, member.points.nearest .+ offset)
    end

    total_N = sum(m.points.N for m in members)
    total_Q = sum(m.points.Q for m in members)

    points = NearestDiscretisation(
        zeros(T, 3, total_N),
        zeros(T, 3, total_Q),
        colony_nearest;
        location=location,
        orientation=orientation
    )

    colony = Colony(members, points, member_force_indices, member_quad_indices)
    update_boundary!(colony, zero(T))
    colony
end

function update_boundary!(colony::Colony, t::T) where {T <: Number}
    for (i, member) in enumerate(colony.members)
        @unpack force_pts, quad_pts, velocity, location, orientation = member.points
        update_boundary!(member, t)

        f_start = colony.member_force_indices[i]
        q_start = colony.member_quad_indices[i]

        @views begin
            colony.points.force_pts[:, f_start:f_start + member.points.N - 1] .= location .+ orientation * force_pts
            colony.points.quad_pts[:,  q_start:q_start + member.points.Q - 1] .= location .+ orientation * quad_pts
            colony.points.velocity[:,  f_start:f_start + member.points.N - 1] .= orientation * velocity
        end
    end
end