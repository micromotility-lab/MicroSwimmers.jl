const I3 = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]
const ex = @SVector [1., 0., 0.]
const ey = @SVector [0., 1., 0.]
const ez = @SVector [0., 0., 1.]

function smooth_max(a, b, k)
    m = max(a, b)
    m + log(exp(k*(a - m)) + exp(k*(b - m))) / k
end

function rotation_matrix(axis::AbstractVector{T}, angle::T) where T
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

function rotation_matrix(axis::AbstractVector{Float64}, angle::T) where T
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

skew_symmetric_static(x::T) where{T} = @SMatrix [0.0   -x[3]   x[2];
                                                 x[3]   0.0   -x[1];
                                                -x[2]   x[1]   0.0]

# Cumulative trapezoid integral: returns I[i] = ∫_{x[1]}^{x[i]} y(s) ds
function cumulative_trapz!(
    out::AbstractVector{T},
    x::AbstractVector{T},
    y::AbstractVector{T}
) where {T}
    @assert length(out) == length(x) == length(y)
    @inbounds begin
        out[1] = zero(T)
        for i in eachindex(x)[2:end]
            dx = x[i] - x[i-1]
            out[i] = out[i-1] + (y[i-1] + y[i]) * (dx/2)
        end
    end
    out
end

struct Helix{T <: Number}
    x0::T
    y0::T
    z0::T
    v::T
    ω::T
    θ::T
    ϕ::T
    r::T
    ψ::T
end

function Base.Tuple(h::Helix)
    (h.x0, h.y0, h.z0, h.v, h.ω, h.θ, h.ϕ, h.r, h.ψ)
end

function helix(ts, p)
    x0, y0, z0, v, ω, θ, ϕ, r, ψ = p
    X0 = @SVector [x0, y0, z0]
    a = @SVector [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    u = abs(a[1]) > 0.9 ? ey : ex

    proj = u - dot(a, u)*a
    e1 = proj / norm(proj)
    e2 = cross(a, e1)

    [X0 + v*t*a + r*(cos(ω*t + ψ)*e1 + sin(ω*t + ψ)*e2) for t in ts]
end

(h::Helix)(ts) = helix(ts, Tuple(h))

translate_helix(h::Helix, x) = Helix(
    h.x0 + x[1], h.y0 + x[2], h.z0 + x[3],
    h.v, h.ω, h.θ, h.ϕ, h.r, h.ψ
)

# Helix quantities
radius(h::Helix) = h.r
axis_polar_angle(h::Helix) = mod(h.θ, 2π)
axis_azimuthal_angle(h::Helix) = mod(h.ϕ + π, 2π) - π
axis_velocity(h::Helix) = h.v
axis_angular_velocity(h::Helix) = h.ω
pitch(h::Helix) = 2π * h.v / h.ω                     # P
pitch_angle(h::Helix) = atan(h.v/(h.ω*h.r))     # α
curvature(h::Helix) =   h.r / (h.r^2 + (h.v/h.ω)^2)
torsion(h::Helix) = (h.v/h.ω) / (h.r^2 + (h.v/h.ω)^2)
initial_point(h::Helix) = @SVector [h.x0, h.y0, h.z0]
axis_vector(h::Helix) = @SVector [sin(h.θ)*cos(h.ϕ), sin(h.θ)*sin(h.ϕ), cos(h.θ)]
chirality_sign(h::Helix) = sign(h.ω)             # +1 right-handed, -1 left-handed
