abstract type Swimmer <: FluidBoundary end

function move!(S::Swimmer, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::T) where {T <: Number}
    update_boundary!(S, t)
    S.config.location = x0
    S.config.orientation = B
end

mutable struct UniFlagellate{T <: Number} <: Swimmer
    body::CellBody
    flagellum::Flagellum
    config::Configuration
    ϵ::T
end

function UniFlagellate(
    body::CellBody=SphericalBody(a=0.2, N=53, Q=253), 
    flagellum::Flagellum=SymmetricFlagellum(N=27, Q=127, location=Vec3(0.2, 0., 0));
    location=SVector(0., 0., 0.),
    orientation=I3,
    ϵ=0.01
)
    nearest = [body.nearest; body.N .+ flagellum.nearest]
    uf = UniFlagellate(
        body,
        flagellum,
        Configuration(
            location,
            orientation,
            zeros(3, size(body.config.force_pts, 2) + size(flagellum.config.force_pts,2)),
            zeros(3, body.N + size(flagellum.config.velocity, 2)),
            zeros(3, size(body.config.quad_pts, 2)  + size(flagellum.config.quad_pts, 2)),
            nearest
        ),
        ϵ
    )

    @views begin
        uf.config.force_pts[:, 1:body.N] .= body.config.force_pts
        uf.config.quad_pts[:,  1:body.Q] .= body.config.quad_pts
    end
    update_boundary!(uf, 0.)
    uf
end

function update_boundary!(uf::UniFlagellate, t::T) where {T <: Number}
    update_boundary!(uf.flagellum, t)
    
    # unpack flagellum configuration
    @unpack force_pts, quad_pts, velocity, location, orientation = uf.flagellum.config

    # don't include the connection point which is already part of the body
    @views begin
        uf.config.force_pts[:, uf.body.N+1:end] .= location .+ orientation * force_pts
        # uf.config.force_pts[:, uf.body.N+1:end] .= location .+ orientation * force_pts # [:,2:end]
        uf.config.quad_pts[:,  uf.body.Q+1:end] .= location .+ orientation * quad_pts # [:, uf.flagellum.Q_cs+1:end]
        uf.config.velocity[:,  uf.body.N+1:end] .= orientation * velocity # [:, 2:end]
    end 
end

struct Flagellate{T <: Number} <: Swimmer
    body::CellBody
    flagella::Vector{Flagellum}
    config::Configuration
    ϵ::T

    force_pt_indices::Vector{Int} # locate the force points for each flagellum in the config structure
    quad_pt_indices::Vector{Int}  # locate the quad points 
end

function Flagellate(
    body::CellBody, 
    flagella::Vector{<:Flagellum};
    location=SVector(0., 0., 0.),
    orientation=I3,
    ϵ=0.01
)
    force_pt_indices = [body.N + 1]
    quad_pt_indices  = [body.Q + 1]
    for i in eachindex(flagella[1:end-1])
        push!(force_pt_indices, body.N + 1 + sum(flagella[j].N - 1 for j = 1:i))
        push!(quad_pt_indices, body.Q + 1 + sum(flagella[j].Q - 1 for j in 1:i))
    end

    flagella_nearest = []
    for (i, flagellum) in enumerate(flagella)
        if i == 1
            append!(flagella_nearest, flagellum.config.nearest[2:end] .- 1)
        else
            append!(flagella_nearest, flagellum.config.nearest[2:end] .- 1 .+ sum(f.N - 1 for f in flagella[1:i-1]))
        end
    end

    nearest = [body.config.nearest; body.N .+ flagella_nearest]
    
    f = Flagellate(
        body,
        flagella,
        Configuration(
            location,
            orientation,
            zeros(3, body.N + sum(f.N - 1 for f in flagella)),
            zeros(3, body.Q + sum(f.Q - 1 for f in flagella)),
            zeros(3, body.N + sum(f.N - 1 for f in flagella)),
            nearest
        ),
        ϵ,
        force_pt_indices,
        quad_pt_indices
    )

    @views begin
        f.config.force_pts[:, 1:body.N] .= body.config.force_pts
        f.config.quad_pts[:,  1:body.Q] .= body.config.quad_pts
    end
    update_boundary!(f, 0.)
    f
end

function update_boundary!(f::Flagellate{T}, t::T) where {T <: Number}
    for (i, flagellum) in enumerate(f.flagella) 
        update_boundary!(flagellum, t)
    
        f_start = f.force_pt_indices[i]
        q_start = f.quad_pt_indices[i]
        
        @unpack force_pts, quad_pts, velocity, location, orientation = flagellum.config
        @views begin
            f.config.force_pts[:, f_start:f_start + flagellum.N-2] .= location .+ orientation * force_pts[:,2:end]
            f.config.quad_pts[:,  q_start:q_start + flagellum.Q-2] .= location .+ orientation * quad_pts[:,2:end]
            f.config.velocity[:,  f_start:f_start + flagellum.N-2] .= orientation * velocity[:,2:end]
        end 
    end
end
