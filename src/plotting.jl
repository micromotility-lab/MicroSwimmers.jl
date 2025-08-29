function ellipsoid(p, a, b, c)
    x, y, z = p
    x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
end

function shifted_ellipsoid(p, a, b, c, d)
    x, y, z = p .- d  # subtract the shift vector d = (dx, dy, dz)
    x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
end

function grooved_ellipsoid(p, a, b, c, d, k=50)
    smooth_max(ellipsoid(p, a, b, c), -shifted_ellipsoid(p, a, b, c, d), k)
end

function gen_mesh(body::EllipsoidBody)
    @unpack a, b, c = body
    Mesh(
        p -> sum(x -> x^2, p ./ [a, b, c]) - 1, 
        origin=Vec3f(-a, -b, -c),
        widths=Vec3f(2a, 2b, 2c),
        samples=(150,150,150)
    )
end 

function gen_mesh(body::EllipsoidalGroovedBody)
    @unpack a, b, c, center_ell2 = body
    Mesh(
        p -> grooved_ellipsoid(p, a, b, c, center_ell2), 
        origin=Vec3f(-a, -b, -c), 
        widths=Vec3f(2a, 2b, 2c),   
        samples=(150,150,150)
    )
end 