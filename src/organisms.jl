function Chlamydomonas(
    a=5.,
    b=4.,
    c=4.,
    L=10.,
    C=-2.5,
    α=π/8,   # attachment angle of flagella in x-y plane relative to x-axis
    N_body=117,
    Q_body=913,
    N_f=17,
    Q_f=111
)

    body = EllipsoidBody(a, b, c, N_body, Q_body)
    model = PlanarFlagellum(L, C, 0.3, 0.15, 2π, 2π, 2π)

    flagella = [
        Flagellum(
            model, N_f, Q_f,
            location=SVector(a*cos(α), b*sin(α), 0),
            orientation=rotation_matrix([1.,0.,0.], 1.0π)
        ), 
        Flagellum(
            model, N_f, Q_f,
            location=SVector(a*cos(α), -b*sin(α), 0)
        )
    ]

    Flagellate(
        body,
        flagella
    )
end
