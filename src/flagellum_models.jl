abstract type FlagellumModel end

## 2D flagella models

struct PlanarFlagellum{T <: Number} <: FlagellumModel
    L::T
    C::T
    R₀::T
    R₁::T
    k::T
    ϕ::T
    ω::T
end

@inline function get_sincos(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ₁ = m.R₀ + m.R₁*sin(m.k*s)
    θ = m.C*s + θ₁*cos(m.ω*t + m.ϕ*s)
    (sin(θ), cos(θ))
end

@inline function get_sincosω(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ₁ = m.R₀ + m.R₁*sin(m.k*s)
    θ = m.C*s + θ₁*cos(m.ω*t + m.ϕ*s)
    (sin(θ), cos(θ), m.ω*θ₁*sin(m.ω*t + m.ϕ*s))
end


function (m::PlanarFlagellum)(pts::Matrix{T}, t::T) where {T <: Number}
    N = size(pts,2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev = get_sincos(s_prev, t, m)
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ = get_sincos(s, t, m)

        pts[1,i] = pts[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        pts[2,i] = pts[2,i-1] + (sinθ_prev + sinθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
    end
end


function (m::PlanarFlagellum)(pts::Matrix{T}, velocity::Matrix{T}, t::T) where {T <: Number}
    N = size(pts,2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev, ωθ₁sin_prev = get_sincosω(s_prev, t, m)  
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, ωθ₁sin = get_sincosω(s, t, m)  

        pts[1,i] = pts[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        pts[2,i] = pts[2,i-1] + (sinθ_prev + sinθ) * half_L_ds
        velocity[1,i] = velocity[1,i-1] + (ωθ₁sin_prev*sinθ_prev + ωθ₁sin*sinθ) * half_L_ds
        velocity[2,i] = velocity[2,i-1] - (ωθ₁sin_prev*cosθ_prev + ωθ₁sin*cosθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
        ωθ₁sin_prev = ωθ₁sin
    end
end

function (m::PlanarFlagellum)(tube_pts::Matrix{T}, tube_vel::Matrix{T}, M::Int, t::T; radius::T=0.01) where {T <: Number}
    N = size(tube_pts, 2) ÷ M
    ds = 1 / (N - 1)
    half_L_ds = 0.5 * m.L * ds

    # Precompute circle angles
    angles = range(0, 2π, length=M+1)[1:end-1]  # avoid duplication at 2π
    cosα = cos.(angles)
    sinα = sin.(angles)

    # Initialise centerline position and velocity
    x, y = zero(T), zero(T)
    vx, vy = zero(T), zero(T)

    s_prev = 0.0
    sinθ_prev, cosθ_prev, ωθ₁sin_prev = get_sincosω(s_prev, t, m)

    @inbounds for i in 1:N
        s = (i - 1) * ds
        sinθ, cosθ, ωθ₁sin = get_sincosω(s, t, m)

        # Update centerline position (x, y)
        if i > 1
            dx = (cosθ_prev + cosθ) * half_L_ds
            dy = (sinθ_prev + sinθ) * half_L_ds
            x += dx
            y += dy

            # Update centerline velocity (vx, vy)
            dvx = (ωθ₁sin_prev*sinθ_prev + ωθ₁sin*sinθ) * half_L_ds
            dvy = -(ωθ₁sin_prev*cosθ_prev + ωθ₁sin*cosθ) * half_L_ds
            vx += dvx
            vy += dvy
        end

        # Frame vectors (fixed in the x-y plane)
        normal  = SVector(-sinθ, cosθ, 0.0)
        binorm  = SVector(0.0, 0.0, 1.0)
        center  = SVector(x, y, 0.0)
        v_center = SVector(vx, vy, 0.0)

        # Angular velocity vector
        ω_vec = ωθ₁sin * binorm  # binorm is z-hat

        # Fill tube points and velocities
        for j in 1:M
            offset = radius * (cosα[j] * normal + sinα[j] * binorm)
            idx = (i - 1) * M + j
            tube_pts[:, idx] = center + offset
            tube_vel[:, idx] = v_center + cross(ω_vec, offset)
        end

        # Save for next step
        s_prev = s
        sinθ_prev = sinθ
        cosθ_prev = cosθ
        ωθ₁sin_prev = ωθ₁sin
    end
end

function (m::PlanarFlagellum)(tube_pts::Matrix{T}, M::Int, t::T; radius::T=0.01) where {T <: Number}
    N = size(tube_pts, 2) ÷ M
    ds = 1 / (N - 1)
    half_L_ds = 0.5 * m.L * ds

    # Precompute circle angles
    angles = range(0, 2π, length=M+1)[1:end-1]  # avoid duplication at 2π
    cosα = cos.(angles)
    sinα = sin.(angles)

    # Set initial centerline point
    x, y = zero(T), zero(T)
    s_prev = 0.0
    sinθ_prev, cosθ_prev = get_sincos(s_prev, t, m)

    @inbounds for i in 1:N
        # Arclength
        s = (i - 1) * ds
        sinθ, cosθ = get_sincos(s, t, m)

        if i > 1
            x += (cosθ_prev + cosθ) * half_L_ds
            y += (sinθ_prev + sinθ) * half_L_ds
        end

        # Normal and binormal (fixed in plane)
        normal = SVector(-sinθ, cosθ, 0.0)
        binorm = SVector(0.0, 0.0, 1.0)

        # Center point of this cross-section
        center = SVector(x, y, 0.0)

        # Fill M points around the circle
        for j in 1:M
            offset =  radius * (cosα[j] * normal + sinα[j] * binorm)
            idx = (i - 1) * M + j
            tube_pts[:, idx] = center + offset
        end

        # Save for next iteration
        s_prev = s
        sinθ_prev = sinθ
        cosθ_prev = cosθ
    end
end

(m::PlanarFlagellum)(points::NearestDiscretisation, t::T) where {T <: Number} = m.model(points.force_pts, points.velocity, t)
function (m::PlanarFlagellum)(points::TubeFlagellumNearestDiscretisation, t::T) where {T <: Number}
    @unpack force_pts, velocity, quad_pts, N_cs, Q_cs, radius = points
    m(force_pts, velocity, N_cs, t, radius=radius)
    m(quad_pts, Q_cs, t, radius=radius)
end


## 3D flagella models

struct QuasiPlanarFlagellum{T <: Number} <: FlagellumModel
    L::T      # Flagellum length
    ω::T      # Beat frequency
    A::T     # Amplitude
    δ::T      # Amplitude modulation
    λ::T      # Wavelength
    C::T      # Static curvature
    C_z::T    # Out-of-plane modulation
end

@inline function get_sincosinvsqrt(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    θ = m.A * (1 - exp(-s*m.L/m.δ)) * sin(m.ω*t - 2π*s*m.L/m.λ) + m.C*s*m.L
    (sin(θ), cos(θ), 1 / sqrt(1 + (s*m.C_z)^2))
end

@inline function get_sincosdotinvsqrt(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    θ = m.A * (1 - exp(-s*m.L/m.δ)) * sin(m.ω*t - 2π*s*m.L/m.λ) + m.C*s*m.L
    θdot = m.ω*m.A*(1 - exp(-s*m.L/m.δ)) * cos(m.ω*t - 2π*s*m.L/m.λ)
    (sin(θ), cos(θ), θdot, 1 / sqrt(1 + (s*m.C_z)^2))
end

function (m::QuasiPlanarFlagellum)(pts::Matrix{T}, t::T) where {T <: Number}
    N = size(pts, 2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev, invsqrt_prev = get_sincosinvsqrt(s_prev, t, m)

    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, invsqrt = get_sincosinvsqrt(s, t, m)

        pts[1,i] = pts[1,i-1] + (invsqrt_prev*cosθ_prev + invsqrt*cosθ) * half_L_ds
        pts[2,i] = pts[2,i-1] + (invsqrt_prev*sinθ_prev + invsqrt*sinθ) * half_L_ds
        pts[3,i] = pts[3,i-1] + m.C_z*s*(invsqrt_prev + invsqrt) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
        invsqrt_prev = invsqrt
    end
end
