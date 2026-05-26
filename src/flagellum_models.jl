abstract type FlagellumModel end

function (m::FlagellumModel)(points::NearestDiscretisation, t::T) where {T <: Number}
    m(points.force_pts, points.velocity, t; include_endpoints=false)
    m(points.quad_pts, t; include_endpoints=true)
end

# If you also want quad point velocities input a pre-allocated matrix quad_velocities
function (m::FlagellumModel)(points::NearestDiscretisation, quad_velocities::Matrix{T}, t::T) where {T <: Number}
    m(points.force_pts, points.velocity, t; include_endpoints=false)
    m(points.quad_pts, quad_velocities, t; include_endpoints=true)
end

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
    t::Number; include_endpoints=false
)
    f_points = @view points[:,1:N_f]
    f_velocities = @view velocities[:,1:N_f]
    m(f_points, f_velocities, t, include_endpoints=include_endpoints)  # Fill flagellum points and velocities
    
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

function get_s0_and_ds(T::Type, N::Int, include_endpoints::Bool)
    if include_endpoints
        return (T(0),  T(1) / T(N-1))
    else
        return (T(1) / T(N+1), T(1) / T(N+1))
    end
end  

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

# an input vector is filled with the tangent angle
(m::PlanarFlagellum)(points::AbstractVector{T}, t::T) where {T <: Number} = tangent_angle(range(0, L, size(points,2)), t, m)

# an input matrix is filled with cartesian points
function (m::PlanarFlagellum)(points::AbstractMatrix{T}, t::T; include_endpoints::Bool=false) where {T <: Number}
    N = size(points, 2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 
    half_L_ds = 0.5*m.L*ds
    
    sinθ_prev, cosθ_prev = orientation_integrands(s_prev, t, m)

    if include_endpoints == false
        sinθ0, cosθ0 = orientation_integrands(zero(T), t, m)
        points[1,1] = (cosθ0 + cosθ_prev) * half_L_ds
        points[2,1] = (sinθ0 + sinθ_prev) * half_L_ds
    end

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


function (m::PlanarFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T; include_endpoints::Bool=false) where {T <: Number}
    N = size(points,2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 

    half_L_ds = 0.5*m.L*ds

    # s_prev = T(0.0)
    sinθ_prev, cosθ_prev, ωθ₁sin_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    if include_endpoints == false
        sinθ0, cosθ0, ωθ₁sin0 = orientation_and_velocity_integrands(zero(T), t, m)
        points[1,1] = (cosθ0 + cosθ_prev) * half_L_ds
        points[2,1] = (sinθ0 + sinθ_prev) * half_L_ds
        velocities[1,1] = (ωθ₁sin_prev*sinθ_prev + ωθ₁sin0*sinθ0) * half_L_ds
        velocities[1,2] =  - (ωθ₁sin_prev*cosθ_prev + ωθ₁sin0*cosθ0) * half_L_ds
    end
    
    
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


mutable struct StandingWaveFlagellum{T <: Number} <: FlagellumModel
    # Arclength discretization
    L::T
    C::T
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
    m.C*s + real(exp(1im*m.ω*t)*θ_space + exp(-1im*m.ω*t)*conj(θ_space))
end

@inline function tangent_angle_and_velocity(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ_space = m.A01*exp(1.0im*m.ϕ01)*sin(π*s/2) + m.A11*exp(1.0im*m.ϕ11)*sin(3π*s/2) + m.A21*exp(1.0im*m.ϕ21)*sin(5π*s/2) + m.A31*exp(1.0im*m.ϕ31)*sin(7π*s/2)
    m.C*s + real(exp(1im*m.ω*t)*θ_space + exp(-1im*m.ω*t)*conj(θ_space)), real(1im*m.ω*exp(1im*m.ω*t)*θ_space - 1im*m.ω*exp(-1im*m.ω*t)*conj(θ_space)) 
end


@inline function orientation_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ = tangent_angle(s, t, m)
    (sin(θ), cos(θ))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    θ, θdot = tangent_angle_and_velocity(s, t, m)
    (sin(θ), cos(θ), θdot)
end

# an input vector is filled with the tangent angle of the centreline
function (m::StandingWaveFlagellum)(points::AbstractVector{T}, t::T) where {T <: Number}
    s = range(0, 1, length(points))
    for i in eachindex(points)
        points[i] = tangent_angle(s[i], t, m)
    end
end

# a 3xN input matrix is filled with the cartesian coordinates of the centerline
function (m::StandingWaveFlagellum)(points::AbstractMatrix{T}, t::T; include_endpoints=false) where {T <: Number}
    N = size(points, 2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 
    half_L_ds = 0.5*m.L*ds
    
    # s_prev = T(0.0)
    sinθ_prev, cosθ_prev = orientation_integrands(s_prev, t, m)

    if include_endpoints == false
        sinθ0, cosθ0 = orientation_integrands(zero(T), t, m)
        points[1,1] = (cosθ0 + cosθ_prev) * half_L_ds
        points[2,1] = (sinθ0 + sinθ_prev) * half_L_ds
    end
    half_L_ds = 0.5*m.L*ds

    # s_prev = T(0.0)
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


function (m::StandingWaveFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T; include_endpoints=false) where {T <: Number}
    N = size(points,2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 

    half_L_ds = 0.5*m.L*ds

    sinθ_prev, cosθ_prev, θdot_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    if include_endpoints == false
        sinθ0, cosθ0, θdot0 = orientation_and_velocity_integrands(zero(T), t, m)
        points[1,1] = (cosθ0 + cosθ_prev) * half_L_ds
        points[2,1] = (sinθ0 + sinθ_prev) * half_L_ds
        velocities[1,1] = (θdot_prev*sinθ_prev + θdot0*sinθ0) * half_L_ds
        velocities[1,2] =  - (θdot_prev*cosθ_prev + θdot0*cosθ0) * half_L_ds
    end
    
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
