function Flagellum(
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, 2π, 2π, 0.0),
    N=23, 
    Q=127; 
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = NearestDiscretisation(N, Q; location=SVector{3}(location), orientation=orientation)
    
    f = Flagellum(model, points)
    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end

function Flagellum(::Type{T},
    model,
    N, 
    Q; 
    location=SVector{3,T}(0, 0, 0), orientation=SMatrix{3,3,T}(I)
) where {T <: Number}
    points = NearestDiscretisation(
        T, N, Q; 
        location=SVector{3,T}(location), orientation=SMatrix{3,3,T}(orientation)
    )
    
    f = Flagellum(model, points)
    update_boundary!(f, T(0.0))
    nearest_neighbour!(f.points)
    f
end

function VanedFlagellum(
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, 2π, 2π, 0.0),
    N=23, 
    Q=127,
    N_v=10,
    N_start=5,
    N_height=3;
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = VanedFlagellumNearestDiscretisation(
        N, Q, N_v, N_start, N_height, 
        location=SVector{3}(location), orientation=orientation
    )

    f = Flagellum(model, points)
    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end


function TubeFlagellum( 
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, -2π, 2π, 0.0),
    N=23, 
    N_cs=5,
    Q=127,
    Q_cs=5; 
    location=SVector(0., 0., 0.),
    orientation=I3,
    radius=0.1
)
    points = TubeFlagellumNearestDiscretisation(
        N, N_cs, Q, Q_cs; 
        location=SVector{3}(location), 
        orientation=orientation,
        radius=radius
    )
    f = Flagellum(model, points)

    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end

function TubeFlagellum(::Type{T}, 
    model,  
    N, 
    N_cs,
    Q,
    Q_cs; 
    location=SVector{3,T}(0, 0, 0), orientation=SMatrix{3,3,T}(I),
    radius=T(0.1)
) where {T <: Number}
    points = TubeFlagellumNearestDiscretisation(
        N, N_cs, Q, Q_cs; 
        location=SVector{3, T}(location), 
        orientation=SMatrix{3,3,T}(orientation),
        radius=radius
    )
    f = Flagellum(model, points)

    update_boundary!(f, zero(T))
    nearest_neighbour!(f.points)
    f
end

function LineTubeFlagellum( 
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, -2π, 2π, 0.0),
    N=23, 
    Q=127,
    Q_cs=5; 
    location=SVector(0., 0., 0.),
    orientation=I3,
)
    points = LineTubeFlagellumNearestDiscretisation(N, Q, Q_cs; 
        location=SVector{3}(location), orientation=orientation
    )
    f = Flagellum(model, points)

    update_boundary!(f, 0.)
    nearest_neighbour!(f.points)
    f
end