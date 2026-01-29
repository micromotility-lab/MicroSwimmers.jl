abstract type FlagellumModel end

function (m::FlagellumModel)(points::NearestDiscretisation, t::T) where {T <: Number}
    m(points.force_pts, points.velocity, t)
    m(points.quad_pts, t)
end

function (m::FlagellumModel)(points::TubeFlagellumNearestDiscretisation, t::T) where {T <: Number}
    @unpack force_pts, velocity, quad_pts, N_cs, Q_cs, radius = points
    m(force_pts, velocity, N_cs, t, radius=radius)
    m(quad_pts, Q_cs, t, radius=radius)
end

function (m::FlagellumModel)(points::LineTubeFlagellumNearestDiscretisation, t::T) where {T <: Number}
    @unpack force_pts, velocity, quad_pts, Q_cs, radius = points
    m(force_pts, velocity, t)
    m(quad_pts, Q_cs, t, radius=radius)
end

# function(m::FlagellumModel)(points::VanedFlagellumNearestDiscretisation, t::T) where {T <: Number}
#     @unpack force_pts, velocity, quad_pts, N_f, Q_f, N_v, N_start, N_height, Q_v, Q_start, Q_height = points

#     m(force_pts, velocity, N_f, N_v, N_start, N_height, t)
#     m(quad_pts, Q_f, Q_v, Q_start, Q_height, t)
# end

"""Flagellum with a vane (only extends in the z direction currently)"""
function (m::FlagellumModel)(
    points::AbstractMatrix,
    N_f::Int, N_v::Int, start::Int, height::Int, 
    t::Number
)
    f_points = @view points[:,1:N_f]
    m(f_points, t)  # Fill flagellum points and velocities
    # Fill vane
    for i in 1:height
        cstart = N_f + (i-1)*N_v + 1
        cend   = N_f + i*N_v
        
        @views begin
            points[1:2, cstart:cend] .= points[1:2, start:start+N_v-1]
            points[3, cstart:cend]   .= -i*m.L / N_f
        end
    end
end

function (m::FlagellumModel)(
    points::AbstractMatrix, velocities::AbstractMatrix, 
    N_f::Int, N_v::Int, start::Int, height::Int, 
    t::Number
)
    f_points = @view points[:,1:N_f]
    f_velocities = @view velocities[:,1:N_f]
    m(f_points, f_velocities, t)  # Fill flagellum points and velocities
    
    # Fill vane
    for i in 1:height
        cstart = N_f + (i-1)*N_v + 1
        cend   = N_f + i*N_v
        
        @views begin
            points[1:2, cstart:cend] .= points[1:2, start:start+N_v-1]
            points[3, cstart:cend]   .= -i*m.L / N_f
            velocities[1:2, cstart:cend]  .= velocities[1:2, start:start+N_v-1]
        end
    end
end



## 2D flagella models

mutable struct PlanarFlagellum{T <: Number} <: FlagellumModel
    L::T
    C::T
    R₀::T
    R₁::T
    k::T
    ϕ::T
    ω::T
    δ::T
end

@inline function tangent_angle(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ₁ = m.R₀ + m.R₁*sin(m.k*s)
    m.C*s + θ₁*cos(m.ω*t - m.ϕ*s + m.δ)
end

@inline function tangent_angle_and_velocity(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ₁ = m.R₀ + m.R₁*sin(m.k*s)
    m.C*s + θ₁*cos(m.ω*t - m.ϕ*s + m.δ), m.ω*θ₁*sin(m.ω*t - m.ϕ*s + m.δ)
end


@inline function orientation_integrands(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ = tangent_angle(s, t, m)
    (sin(θ), cos(θ))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    θ, θdot = tangent_angle_and_velocity(s, t, m)
    (sin(θ), cos(θ), θdot)
end

(m::PlanarFlagellum)(points::AbstractVector{T}, t::T) where {T <: Number} = tangent_angle(range(0, L, size(points,2)), t, m)

function (m::PlanarFlagellum)(points::AbstractMatrix{T}, t::T) where {T <: Number}
    N = size(points,2)
    ds = T(1/(N-1))
    half_L_ds = 0.5*m.L*ds

    s_prev = T(0.0)
    sinθ_prev, cosθ_prev = orientation_integrands(s_prev, t, m)
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ = orientation_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (sinθ_prev + sinθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
    end
end


function (m::PlanarFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T) where {T <: Number}
    # tT = T(t) # promote t for autodiff
    N = size(points,2)
    ds = T(1/(N-1))
    half_L_ds = 0.5*m.L*ds

    s_prev = T(0.0)
    sinθ_prev, cosθ_prev, ωθ₁sin_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, ωθ₁sin = orientation_and_velocity_integrands(s, t, m)  

        points[1,i] = points[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (sinθ_prev + sinθ) * half_L_ds
        velocities[1,i] = velocities[1,i-1] + (ωθ₁sin_prev*sinθ_prev + ωθ₁sin*sinθ) * half_L_ds
        velocities[2,i] = velocities[2,i-1] - (ωθ₁sin_prev*cosθ_prev + ωθ₁sin*cosθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
        ωθ₁sin_prev = ωθ₁sin
    end
end


function (m::PlanarFlagellum)(tube_points::AbstractMatrix{T}, M::Int, t::T; radius::T=0.01) where {T <: Number}
    N = size(tube_points, 2) ÷ M
    ds = 1 / (N - 1)
    half_L_ds = 0.5 * m.L * ds

    # Precompute circle angles
    angles = range(0, 2π, length=M+1)[1:end-1]  # avoid duplication at 2π
    cosα = cos.(angles)
    sinα = sin.(angles)

    # Set initial centerline point
    x, y = zero(T), zero(T)
    s_prev = 0.0
    sinθ_prev, cosθ_prev = orientation_integrands(s_prev, t, m)

    @inbounds for i in 1:N
        # Arclength
        s = (i - 1) * ds
        sinθ, cosθ = orientation_integrands(s, t, m)

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
            tube_points[:, idx] = center + offset
        end

        # Save for next iteration
        s_prev = s
        sinθ_prev = sinθ
        cosθ_prev = cosθ
    end
end

function (m::PlanarFlagellum)(tube_points::AbstractMatrix{T}, tube_vel::AbstractMatrix{T}, M::Int, t::T; radius::T=0.01) where {T <: Number}
    N = size(tube_points, 2) ÷ M
    ds = 1 / (N - 1)
    half_L_ds = 0.5 * m.L * ds

    # Precompute circle angles
    angles = range(0, 2π, length=M+1)[1:end-1]  # avoid duplication at 2π
    cosα = cos.(angles)
    sinα = sin.(angles)

    # Initialise centerline position and velocities
    x, y = zero(T), zero(T)
    vx, vy = zero(T), zero(T)

    s_prev = 0.0
    sinθ_prev, cosθ_prev, ωθ₁sin_prev = orientation_and_velocity_integrands(s_prev, t, m)

    @inbounds for i in 1:N
        s = (i - 1) * ds
        sinθ, cosθ, ωθ₁sin = orientation_and_velocity_integrands(s, t, m)

        # Update centerline position (x, y)
        if i > 1
            dx = (cosθ_prev + cosθ) * half_L_ds
            dy = (sinθ_prev + sinθ) * half_L_ds
            x += dx
            y += dy

            # Update centerline velocities (vx, vy)
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
            tube_points[:, idx] = center + offset
            tube_vel[:, idx] = v_center + cross(ω_vec, offset)
        end

        # Save for next step
        s_prev = s
        sinθ_prev = sinθ
        cosθ_prev = cosθ
        ωθ₁sin_prev = ωθ₁sin
    end
end

mutable struct StandingWaveFlagellum{T <: Number} <: FlagellumModel
    # Arclength discretization
    L::T
    A01::T
    ϕ01::T
    A11::T
    ϕ11::T
    A21::T
    ϕ21::T
    A31::T
    ϕ31::T
    ω::T
end

@inline function tangent_angle(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ_space = m.A01*exp(1.0im*m.ϕ01)*sin(π*s/2) + m.A11*exp(1.0im*m.ϕ11)*sin(3π*s/2) + m.A21*exp(1.0im*m.ϕ21)*sin(5π*s/2) + m.A31*exp(1.0im*m.ϕ31)*sin(7π*s/2)
    real(exp(1im*m.ω*t)*θ_space + exp(-1im*m.ω*t)*conj(θ_space))
end

@inline function tangent_angle_and_velocity(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ_space = m.A01*exp(1.0im*m.ϕ01)*sin(π*s/2) + m.A11*exp(1.0im*m.ϕ11)*sin(3π*s/2) + m.A21*exp(1.0im*m.ϕ21)*sin(5π*s/2) + m.A31*exp(1.0im*m.ϕ31)*sin(7π*s/2)
    real(exp(1im*m.ω*t)*θ_space + exp(-1im*m.ω*t)*conj(θ_space)), real(1im*m.ω*exp(1im*m.ω*t)*θ_space - 1im*m.ω*exp(-1im*m.ω*t)*conj(θ_space)) 
end


@inline function orientation_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ = tangent_angle(s, t, m)
    (sin(θ), cos(θ))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ, θdot = tangent_angle_and_velocity(s, t, m)
    (sin(θ), cos(θ), θdot)
end

function (m::StandingWaveFlagellum)(points::AbstractVector{T}, t::T) where {T <: Number}
    s = range(0, 1, length(points))
    for i in eachindex(points)
        points[i] = tangent_angle(s[i], t, m)
    end
end

function (m::StandingWaveFlagellum)(points::AbstractMatrix{T}, t::T) where {T <: Number}
    N = size(points,2)
    ds = T(1/(N-1))
    half_L_ds = 0.5*m.L*ds

    s_prev = T(0.0)
    sinθ_prev, cosθ_prev = orientation_integrands(s_prev, t, m)
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ = orientation_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (sinθ_prev + sinθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
    end
end


function (m::StandingWaveFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T) where {T <: Number}
    # tT = T(t) # promote t for autodiff
    N = size(points,2)
    ds = T(1/(N-1))
    half_L_ds = 0.5*m.L*ds

    s_prev = T(0.0)
    sinθ_prev, cosθ_prev, θdot_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, θdot = orientation_and_velocity_integrands(s, t, m)  

        points[1,i] = points[1,i-1] + (cosθ_prev + cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (sinθ_prev + sinθ) * half_L_ds
        velocities[1,i] = velocities[1,i-1] - (θdot_prev*sinθ_prev + θdot*sinθ) * half_L_ds
        velocities[2,i] = velocities[2,i-1] + (θdot_prev*cosθ_prev + θdot*cosθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
        θdot_prev = θdot
    end
end



## 3D flagella models

mutable struct QuasiPlanarFlagellum{T <: Number} <: FlagellumModel
    L::T      # Flagellum length
    ω::T      # Beat frequency
    A::T     # Amplitude
    δ::T      # Amplitude modulation
    λ::T      # Wavelength
    C::T      # Static curvature
    C_z::T    # Out-of-plane modulation
end

@inline function orientation_integrands(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    θ = m.A * (1 - exp(-s*m.L/m.δ)) * sin(m.ω*t - 2π*s*m.L/m.λ) + m.C*s*m.L
    (sin(θ), cos(θ), 1 / sqrt(1 + (s*m.C_z)^2))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    θ = m.A * (1 - exp(-s*m.L/m.δ)) * sin(m.ω*t - 2π*s*m.L/m.λ) + m.C*s*m.L
    θdot = m.ω*m.A*(1 - exp(-s*m.L/m.δ)) * cos(m.ω*t - 2π*s*m.L/m.λ)
    (sin(θ), cos(θ), θdot, 1 / sqrt(1 + (s*m.C_z)^2))
end

function (m::QuasiPlanarFlagellum)(points::AbstractMatrix{T}, t::T) where {T <: Number}
    N = size(points, 2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev, invsqrt_prev = orientation_integrands(s_prev, t, m)
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, invsqrt = orientation_integrands(s, t, m)
        
        points[1,i] = points[1,i-1] + (invsqrt_prev*cosθ_prev + invsqrt*cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (invsqrt_prev*sinθ_prev + invsqrt*sinθ) * half_L_ds
        points[3,i] = points[3,i-1] + m.C_z*s*(invsqrt_prev + invsqrt) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
        invsqrt_prev = invsqrt
    end
end

function (m::QuasiPlanarFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T) where {T <: Number}
    N = size(points, 2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev, θdot_prev,  invsqrt_prev = orientation_and_velocity_integrands(s_prev, t, m)

    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, θdot, invsqrt = orientation_and_velocity_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (invsqrt_prev*cosθ_prev + invsqrt*cosθ) * half_L_ds
        points[2,i] = points[2,i-1] + (invsqrt_prev*sinθ_prev + invsqrt*sinθ) * half_L_ds
        points[3,i] = points[3,i-1] + m.C_z*s*(invsqrt_prev + invsqrt) * half_L_ds

        velocities[1,i] = velocities[1,i-1] - (invsqrt_prev*θdot_prev*sinθ_prev + invsqrt*θdot*sinθ) * half_L_ds
        velocities[2,i] = velocities[2,i-1] + (invsqrt_prev*θdot_prev*cosθ_prev + invsqrt*θdot*cosθ) * half_L_ds

        s_prev       = s
        sinθ_prev    = sinθ 
        cosθ_prev    = cosθ
        θdot_prev    = θdot
        invsqrt_prev = invsqrt
    end
end

# Three dimensional flagellum based on the model used in Suzuki-Tellier et. al 2024

mutable struct ThreeDimensionalFlagellum{T <: Number}
    ## Waveform parameters
    L::Float64      # Flagellum length
    fᵩ::Float64      # Azimuthal beat frequency
    Aᵩ::Float64     # Azimuthal amplitude
    δᵩ::Float64     # Azimuthal amplitude modulation
    λᵩ::Float64     # Azimuthal Wavelength
    Cᵩ::Float64     # Azimuthal static curvature

    f_θ::Float64    # Polar beat frequency
    A_θ::Float64    # Polar amplitude
    δ_θ::Float64    # Polar amplitude modulation
    λ_θ::Float64    # Polar wavelength
    C_θ::Float64    # Polar static curvature

    γ::Float64      # overall phase 
    Δγ::Float64     # relative phase   
end

@inline function get_integrands(s::T, t::T, m::ThreeDimensionalFlagellum) where {T <: Number}
    θ = m.A_θ * (1 - exp(-s*L/m.δ_θ))*sin(2π*m.f_θ*t - 2π*s*L/m.λ_θ + m.γ + m.Δγ) + m.C_θ*s
    ϕ = m.Aᵩ  * (1 - exp(-s*L/m.δᵩ))*sin(2π*m.fᵩ*t - 2π*s*L/m.λᵩ + m.γ) + m.C_θ*s
    (sin(θ), cos(θ), sin(ϕ), cos(ϕ))
end


function (m::ThreeDimensionalFlagellum)(points::AbstractMatrix{T}, t::T) where {T <: Number}
    N = size(points,2)
    ds = 1/(N-1)
    half_L_ds = 0.5*m.L*ds

    s_prev = 0.0
    sinθ_prev, cosθ_prev, sinϕ_prev, cosϕ_prev  = get_integrands(s_prev, t, m)
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sinθ, cosθ, sinϕ_prev, cosϕ_prev = get_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (cosθ_prev*cosϕ_prev + cosθ*cosϕ) * half_L_ds
        points[2,i] = points[2,i-1] + (cosθ_prev*sinϕ_prev + cosθ*sinϕ) * half_L_ds
        points[3,i] = points[3,i-1] + (sinθ_prev + sinθ) * half_L_ds

        s_prev      = s
        sinθ_prev   = sinθ
        cosθ_prev   = cosθ
    end
end

# function (m::ThreeDimensionalFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T)
# end


    
