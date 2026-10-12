const I3 = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]
const ex = @SVector [1., 0., 0.]
const ey = @SVector [0., 1., 0.]
const ez = @SVector [0., 0., 1.]

function rotation_matrix(axis::AbstractVector, angle)
    axis = normalize(axis)  # Ensure it's a unit vector
    x, y, z = axis
    c = cos(angle)
    s = sin(angle)
    C = 1 - c

    @SMatrix [
        c + x^2*C     x*y*C - z*s   x*z*C + y*s;
        y*x*C + z*s   c + y^2*C     y*z*C - x*s;
        z*x*C - y*s   z*y*C + x*s   c + z^2*C
    ]
end

skew_symmetric_static(x::T) where {T} = @SMatrix [0.0   -x[3]   x[2];
                                                 x[3]   0.0   -x[1];
                                                -x[2]   x[1]   0.0]

struct Helix{T <: Number}
    x0::T
    y0::T
    z0::T
    v::T
    omega::T
    theta::T
    phi::T
    r::T
    psi::T
end

function Base.Tuple(h::Helix)
    (h.x0, h.y0, h.z0, h.v, h.omega, h.theta, h.phi, h.r, h.psi)
end

function helix(ts, p)
    x0, y0, z0, v, omega, theta, phi, r, psi = p
    X0 = @SVector [x0, y0, z0]
    a = @SVector [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    u = abs(a[1]) > 0.9 ? ey : ex

    proj = u - dot(a, u)*a
    e1 = proj / norm(proj)
    e2 = cross(a, e1)

    [X0 + v*t*a + r*(cos(omega*t + psi)*e1 + sin(omega*t + psi)*e2) for t in ts]
end

(h::Helix)(ts) = helix(ts, Tuple(h))

translate_helix(h::Helix, x) = Helix(
    h.x0 + x[1], h.y0 + x[2], h.z0 + x[3],
    h.v, h.omega, h.theta, h.phi, h.r, h.psi
)

# Helix quantities
radius(h::Helix) = h.r
axis_polar_angle(h::Helix) = mod(h.theta, 2pi)
axis_azimuthal_angle(h::Helix) = mod(h.phi + pi, 2pi) - pi
axis_velocity(h::Helix) = h.v
axis_angular_velocity(h::Helix) = h.omega
pitch(h::Helix) = 2pi * h.v / h.omega                     # P
pitch_angle(h::Helix) = atan(h.v/(h.omega*h.r))     # alpha
curvature(h::Helix) =   h.r / (h.r^2 + (h.v/h.omega)^2)
torsion(h::Helix) = (h.v/h.omega) / (h.r^2 + (h.v/h.omega)^2)
initial_point(h::Helix) = @SVector [h.x0, h.y0, h.z0]
axis_vector(h::Helix) = @SVector [sin(h.theta)*cos(h.phi), sin(h.theta)*sin(h.phi), cos(h.theta)]
chirality_sign(h::Helix) = sign(h.omega)             # +1 right-handed, -1 left-handed
