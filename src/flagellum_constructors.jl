abstract type Flagellum <: MicroSwimmer end

struct BareFlagellum{M <: FlagellumModel, D <: Discretisation} <: Flagellum
    model::M
    points::D
end

function update_boundary!(f::BareFlagellum, t::T) where {T <: Number}
    f.model(f.points, t)
end

# Provide N and Q for the nearest discretisation
function Flagellum(
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, 2π, 2π, 0.0),
    N::Int=23, 
    Q::Int=127; 
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = NearestDiscretisation(N, Q; location=SVector{3}(location), orientation=orientation)
    
    f = BareFlagellum(model, points)
    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end


struct VanedFlagellum{M <: FlagellumModel, D <: Discretisation} <: Flagellum
    model::M
    points::D  # flagellum points

    N_f::Int
    Q_f::Int

    N_v::Int
    N_start::Int
    N_height::Int

    Q_v::Int
    Q_start::Int
    Q_height::Int
end

function VanedFlagellum(
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, 2π, 2π, 0.0),
    N_f=23, 
    Q_f=127,
    N_v=10,
    N_start=5,
    N_height=3;
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    Q_v = floor(Int, (N_v / N_f) * Q_f) 
    Q_start = ceil(Int, ((N_start-1) / (N_f-1)) * (Q_f-1))
    Q_height = floor(Int, (N_height / N_f) * Q_f)

    points = NearestDiscretisation(
        N_f + N_height*N_v, Q_f + Q_height*Q_v;
        location=SVector{3}(location), orientation=orientation
    )

    vf = VanedFlagellum(
        model, points, 
        N_f, Q_f,
        N_v, N_start, N_height, 
        Q_v, Q_start, Q_height
    )

    update_boundary!(vf, 0.)
    nearest_neighbour!(vf.points)
    vf
end

function update_boundary!(vf::VanedFlagellum, t::T) where {T <: Number}
    @unpack N_f, N_v, N_start, N_height, Q_f, Q_v, Q_start, Q_height = vf
    vf.model(vf.points.force_pts, vf.points.velocity, N_f, N_v, N_start, N_height, t)
    vf.model(vf.points.quad_pts, Q_f, Q_v, Q_start, Q_height, t)
end

get_vane_pts(vf::VanedFlagellum) = hcat(
    vf.points.force_pts[:, vf.N_start:vf.N_start + vf.N_v - 1],  
    vf.points.force_pts[:, vf.N_f+1:end]
)
