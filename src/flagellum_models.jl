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
    R0::T
    R1::T
    k::T
    phi::T
    omega::T
    delta::T
end

@inline function tangent_angle(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    theta1 = m.R0 + m.R1*sin(m.k*s)
    m.C*s + theta1*cos(m.omega*t - m.phi*s + m.delta)
end

@inline function tangent_angle_and_velocity(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    theta1 = m.R0 + m.R1*sin(m.k*s)
    m.C*s + theta1*cos(m.omega*t - m.phi*s + m.delta), m.omega*theta1*sin(m.omega*t - m.phi*s + m.delta)
end


@inline function orientation_integrands(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    theta = tangent_angle(s, t, m)
    (sin(theta), cos(theta))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    theta, thetadot = tangent_angle_and_velocity(s, t, m)
    (sin(theta), cos(theta), thetadot)
end

# an input vector is filled with the tangent angle
(m::PlanarFlagellum)(points::AbstractVector{T}, t::T) where {T <: Number} = tangent_angle(range(0, L, size(points,2)), t, m)

# an input matrix is filled with cartesian points
function (m::PlanarFlagellum)(points::AbstractMatrix{T}, t::T; include_endpoints::Bool=false) where {T <: Number}
    N = size(points, 2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 
    half_L_ds = 0.5*m.L*ds
    
    sintheta_prev, costheta_prev = orientation_integrands(s_prev, t, m)

    if include_endpoints == false
        sintheta0, costheta0 = orientation_integrands(zero(T), t, m)
        points[1,1] = (costheta0 + costheta_prev) * half_L_ds
        points[2,1] = (sintheta0 + sintheta_prev) * half_L_ds
    end

    @inbounds for i in 2:N
        s = s_prev + ds
        sintheta, costheta = orientation_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (costheta_prev + costheta) * half_L_ds
        points[2,i] = points[2,i-1] + (sintheta_prev + sintheta) * half_L_ds

        s_prev      = s
        sintheta_prev   = sintheta
        costheta_prev   = costheta
    end
end


function (m::PlanarFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T; include_endpoints::Bool=false) where {T <: Number}
    N = size(points,2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 

    half_L_ds = 0.5*m.L*ds

    # s_prev = T(0.0)
    sintheta_prev, costheta_prev, omega_theta1_sin_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    if include_endpoints == false
        sintheta0, costheta0, omega_theta1_sin0 = orientation_and_velocity_integrands(zero(T), t, m)
        points[1,1] = (costheta0 + costheta_prev) * half_L_ds
        points[2,1] = (sintheta0 + sintheta_prev) * half_L_ds
        velocities[1,1] = (omega_theta1_sin_prev*sintheta_prev + omega_theta1_sin0*sintheta0) * half_L_ds
        velocities[2,1] =  - (omega_theta1_sin_prev*costheta_prev + omega_theta1_sin0*costheta0) * half_L_ds
    end
    
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sintheta, costheta, omega_theta1_sin = orientation_and_velocity_integrands(s, t, m)  

        points[1,i] = points[1,i-1] + (costheta_prev + costheta) * half_L_ds
        points[2,i] = points[2,i-1] + (sintheta_prev + sintheta) * half_L_ds
        velocities[1,i] = velocities[1,i-1] + (omega_theta1_sin_prev*sintheta_prev + omega_theta1_sin*sintheta) * half_L_ds
        velocities[2,i] = velocities[2,i-1] - (omega_theta1_sin_prev*costheta_prev + omega_theta1_sin*costheta) * half_L_ds

        s_prev      = s
        sintheta_prev   = sintheta
        costheta_prev   = costheta
        omega_theta1_sin_prev = omega_theta1_sin
    end
end


mutable struct StandingWaveFlagellum{T <: Number} <: FlagellumModel
    # Arclength discretization
    L::T
    C::T
    A01::T
    phi01::T
    A11::T
    phi11::T
    A21::T
    phi21::T
    A31::T
    phi31::T
    omega::T
end

@inline function tangent_angle(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    theta_space = m.A01*exp(1.0im*m.phi01)*sin(pi*s/2) + m.A11*exp(1.0im*m.phi11)*sin(3pi*s/2) + m.A21*exp(1.0im*m.phi21)*sin(5pi*s/2) + m.A31*exp(1.0im*m.phi31)*sin(7pi*s/2)
    m.C*s + real(exp(1im*m.omega*t)*theta_space + exp(-1im*m.omega*t)*conj(theta_space))
end

@inline function tangent_angle_and_velocity(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    theta_space = m.A01*exp(1.0im*m.phi01)*sin(pi*s/2) + m.A11*exp(1.0im*m.phi11)*sin(3pi*s/2) + m.A21*exp(1.0im*m.phi21)*sin(5pi*s/2) + m.A31*exp(1.0im*m.phi31)*sin(7pi*s/2)
    m.C*s + real(exp(1im*m.omega*t)*theta_space + exp(-1im*m.omega*t)*conj(theta_space)), real(1im*m.omega*exp(1im*m.omega*t)*theta_space - 1im*m.omega*exp(-1im*m.omega*t)*conj(theta_space)) 
end


@inline function orientation_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    theta = tangent_angle(s, t, m)
    (sin(theta), cos(theta))
end

@inline function orientation_and_velocity_integrands(s::T, t::T, m::StandingWaveFlagellum) where {T <: Number}
    theta, thetadot = tangent_angle_and_velocity(s, t, m)
    (sin(theta), cos(theta), thetadot)
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
    sintheta_prev, costheta_prev = orientation_integrands(s_prev, t, m)

    if include_endpoints == false
        sintheta0, costheta0 = orientation_integrands(zero(T), t, m)
        points[1,1] = (costheta0 + costheta_prev) * half_L_ds
        points[2,1] = (sintheta0 + sintheta_prev) * half_L_ds
    end
    half_L_ds = 0.5*m.L*ds

    # s_prev = T(0.0)
    sintheta_prev, costheta_prev = orientation_integrands(s_prev, t, m)
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sintheta, costheta = orientation_integrands(s, t, m)

        points[1,i] = points[1,i-1] + (costheta_prev + costheta) * half_L_ds
        points[2,i] = points[2,i-1] + (sintheta_prev + sintheta) * half_L_ds

        s_prev      = s
        sintheta_prev   = sintheta
        costheta_prev   = costheta
    end
end


function (m::StandingWaveFlagellum)(points::AbstractMatrix{T}, velocities::AbstractMatrix{T}, t::T; include_endpoints=false) where {T <: Number}
    N = size(points,2)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints) 

    half_L_ds = 0.5*m.L*ds

    sintheta_prev, costheta_prev, thetadot_prev = orientation_and_velocity_integrands(s_prev, t, m)  
    if include_endpoints == false
        sintheta0, costheta0, thetadot0 = orientation_and_velocity_integrands(zero(T), t, m)
        points[1,1] = (costheta0 + costheta_prev) * half_L_ds
        points[2,1] = (sintheta0 + sintheta_prev) * half_L_ds
        velocities[1,1] = - (thetadot_prev*sintheta_prev + thetadot0*sintheta0) * half_L_ds
        velocities[2,1] = (thetadot_prev*costheta_prev + thetadot0*costheta0) * half_L_ds
    end
    
    @inbounds for i in 2:N
        s = s_prev + ds
        sintheta, costheta, thetadot = orientation_and_velocity_integrands(s, t, m)  

        points[1,i] = points[1,i-1] + (costheta_prev + costheta) * half_L_ds
        points[2,i] = points[2,i-1] + (sintheta_prev + sintheta) * half_L_ds
        velocities[1,i] = velocities[1,i-1] - (thetadot_prev*sintheta_prev + thetadot*sintheta) * half_L_ds
        velocities[2,i] = velocities[2,i-1] + (thetadot_prev*costheta_prev + thetadot*costheta) * half_L_ds

        s_prev      = s
        sintheta_prev   = sintheta
        costheta_prev   = costheta
        thetadot_prev = thetadot
    end
end
