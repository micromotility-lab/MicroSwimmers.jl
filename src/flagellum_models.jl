abstract type FlagellumModel <: Model end

# ============================================================================
#  Flagellum waveform models
#
#  Each concrete model supplies ONLY the geometry through two methods:
#
#      unit_tangent(s, t, m)        -> SVector{3,T}                  (t̂)
#      unit_tangent_and_dt(s, t, m) -> (SVector{3,T}, SVector{3,T})  (t̂, ∂t̂/∂t)
#
#  with s the fractional arclength on [0,1] and t̂ a UNIT tangent (so the
#  centreline stays arclength-parameterised: |t̂| ≡ 1 ⇒ BEM nodes evenly
#  spaced).  The shared integrator below turns those into positions and
#  velocities by cumulative trapezoid, handling the base node and the
#  force-point (open) vs quad-point (closed) endpoint rules in one place.
#
# ============================================================================

# ---------------------------------------------------------------------------
#  Shared discretisation + integrator
# ---------------------------------------------------------------------------

# Starting station and step for the cumulative trapezoid.
#   include_endpoints=true  : closed rule, nodes at s = 0, 1/(N-1), …, 1   (quad pts)
#   include_endpoints=false : open rule,   nodes at s = ds, 2ds, …,        (force pts)
function get_s0_and_ds(T::Type, N::Int, include_endpoints::Bool)
    include_endpoints ? (T(0), T(1) / T(N - 1)) : (T(1) / T(N + 1), T(1) / T(N + 1))
end

# position only
@inline function integrate_centreline!(points::AbstractVector{SVector{3,T}},
                                       m::FlagellumModel, t::T;
                                       include_endpoints::Bool) where {T <: Number}
    N = length(points)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints)
    half_L_ds  = T(0.5) * m.L * ds

    τ_prev = unit_tangent(s_prev, t, m)

    if include_endpoints
        points[1] = zero(SVector{3,T})
    else
        τ0 = unit_tangent(zero(T), t, m)              # integrate the first panel from s=0
        points[1] = (τ0 + τ_prev) * half_L_ds
    end

    @inbounds for i in 2:N
        s = s_prev + ds
        τ = unit_tangent(s, t, m)
        points[i] = points[i-1] + (τ_prev + τ) * half_L_ds
        s_prev, τ_prev = s, τ
    end
    nothing
end

# position + velocity
@inline function integrate_centreline!(points::AbstractVector{SVector{3,T}},
                                       velocities::AbstractVector{SVector{3,T}},
                                       m::FlagellumModel, t::T;
                                       include_endpoints::Bool) where {T <: Number}
    N = length(points)
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints)
    half_L_ds  = T(0.5) * m.L * ds

    τ_prev, τdot_prev = unit_tangent_and_dt(s_prev, t, m)

    if include_endpoints
        points[1]     = zero(SVector{3,T})
        velocities[1] = zero(SVector{3,T})
    else
        τ0, τdot0 = unit_tangent_and_dt(zero(T), t, m)
        points[1]     = (τ0  + τ_prev)  * half_L_ds
        velocities[1] = (τdot0 + τdot_prev) * half_L_ds
    end

    @inbounds for i in 2:N
        s = s_prev + ds
        τ, τdot = unit_tangent_and_dt(s, t, m)
        points[i]     = points[i-1]     + (τ_prev  + τ)  * half_L_ds
        velocities[i] = velocities[i-1] + (τdot_prev + τdot) * half_L_ds
        s_prev, τ_prev, τdot_prev= s, τ, τdot
    end
    nothing
end


@inline (m::FlagellumModel)(points::AbstractVector{SVector{3,T}}, t::T;
                            include_endpoints::Bool=false) where {T<:Number} =
    integrate_centreline!(points, m, t; include_endpoints)


@inline (m::FlagellumModel)(points::AbstractVector{SVector{3,T}}, velocities::AbstractVector{SVector{3,T}},
                            t::T; include_endpoints::Bool=false) where {T<:Number} =
    integrate_centreline!(points, velocities, m, t; include_endpoints)

# return a matrix of points on an (s,t) grid
function (m::FlagellumModel)(s_vals::AbstractVector{T}, t_vals::AbstractVector{T}) where {T <: Number}
    pts = Matrix{SVector{3,T}}(undef, length(s_vals), length(t_vals))
    for (n,t) in enumerate(t_vals)
        integrate_centreline!(view(pts, :, n), m, t; include_endpoints=true)
    end
    pts
end


# Discretisation glue (force pts → open rule, quad pts → closed rule).
function (m::FlagellumModel)(disc::NearestDiscretisation, t::T) where {T <: Number}
    m(disc.force_pts, disc.velocity, t; include_endpoints=false)
    m(disc.quad_pts,                 t; include_endpoints=true)
end

# for use if quadrature point velocities are also required
function (m::FlagellumModel)(points::NearestDiscretisation, quad_velocities::Matrix{T}, t::T) where {T <: Number}
    m(points.force_pts, points.velocity, t; include_endpoints=false)
    m(points.quad_pts, quad_velocities, t;  include_endpoints=true)
end

# ===========================================================================
#  1. PlanarFlagellum  — planar travelling wave
#       θ(s,t) = C·s + A(s)·cos(ωt − ϕs + δ),   A(s) = R₀ + R₁·sin(k s)
# ===========================================================================
mutable struct PlanarFlagellum{T <: Number} <: FlagellumModel
    L::T
    C::T        # static curvature (accumulated over fractional s)
    R₀::T       # base amplitude
    R₁::T       # amplitude modulation
    k::T        # amplitude wavenumber
    ϕ::T        # phase gradient (sets wavelength)
    ω::T        # beat frequency
    δ::T        # phase offset
end

@inline function unit_tangent(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    A = m.R₀ + m.R₁*sin(m.k*s)
    θ = m.C*s + A*cos(m.ω*t - m.ϕ*s + m.δ)
    SVector(cos(θ), sin(θ), zero(T))
end

@inline function unit_tangent_and_dt(s::T, t::T, m::PlanarFlagellum) where {T <: Number}
    A = m.R₀ + m.R₁*sin(m.k*s)
    φ = m.ω*t - m.ϕ*s + m.δ
    θ, θdot = m.C*s + A*cos(φ),  -m.ω*A*sin(φ) 
    (SVector(cos(θ), sin(θ), zero(T)), θdot * SVector(-sin(θ), cos(θ), zero(T)))
end

# ===========================================================================
#  2. QuasiPlanarFlagellum — planar travelling beat with a STATIC out-of-plane
#     lean.  t̂ = w·(cosθ, sinθ, C_z·s),  w = 1/√(1+(C_z s)²)  (keeps |t̂|=1).
#     z carries no t ⇒ no z-velocity.  Physical units: δ,λ are lengths.
#       θ(s,t) = A·(1 − e^{−sL/δ})·sin(ωt − 2π sL/λ) + C·sL
# ===========================================================================
mutable struct QuasiPlanarFlagellum{T <: Number} <: FlagellumModel
    L::T
    ω::T        # beat frequency
    A::T        # angular amplitude
    δ::T        # base ramp length
    λ::T        # wavelength
    C::T        # static curvature (1/length)
    C_z::T      # out-of-plane lean (tan of tip elevation)
end

@inline function unit_tangent(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    sL    = s*m.L
    θ     = m.A*(1 - exp(-sL/m.δ))*sin(m.ω*t - 2*T(π)*sL/m.λ) + m.C*sL
    scale = one(T)/sqrt(one(T) + (s*m.C_z)^2)
    scale * SVector(cos(θ), sin(θ), m.C_z*s)
end

@inline function unit_tangent_and_dt(s::T, t::T, m::QuasiPlanarFlagellum) where {T <: Number}
    sL    = s*m.L
    φ     = m.ω*t - 2*T(π)*sL/m.λ
    env   = m.A*(1 - exp(-sL/m.δ))
    θ     = env*sin(φ) + m.C*sL
    θdot   = m.ω*env*cos(φ)                          # true ∂θ/∂t
    scale = one(T)/sqrt(one(T) + (s*m.C_z)^2)
    (scale * SVector(cos(θ), sin(θ), m.C_z*s),
     scale*θdot* SVector(-sin(θ), cos(θ), zero(T)))
end

# ===========================================================================
#  3. ThreeDimensionalFlagellum — full 3D travelling wave.
#     Spherical tangent  t̂ = (cosθ cosϕ, cosθ sinϕ, sinθ)  (unit-norm).
#       ϕ : azimuth (in xy-plane),  θ : elevation above it.
#       each angle = base-ramped travelling wave + static curvature.
#       Δγ (elevation − azimuth phase) sets helicity: 0 planar … π/2 conical.
#     NOTE: this model is parameterised by FREQUENCY f (with explicit 2πf);
#           the others use angular ω. 
# ===========================================================================
mutable struct ThreeDimensionalFlagellum{T <: Number} <: FlagellumModel
    L::T
    fᵩ::T;  Aᵩ::T;  δᵩ::T;  λᵩ::T;  Cᵩ::T          # azimuthal bank
    f_θ::T; A_θ::T; δ_θ::T; λ_θ::T; C_θ::T          # elevation bank
    γ::T                                            # overall phase
    Δγ::T                                           # relative phase (elev − azim)
end

@inline function unit_tangent(s::T, t::T, m::ThreeDimensionalFlagellum) where {T <: Number}
    sL = s*m.L
    θ  = m.A_θ*(1 - exp(-sL/m.δ_θ))*sin(2*T(π)*m.f_θ*t - 2*T(π)*sL/m.λ_θ + m.γ + m.Δγ) + m.C_θ*sL
    ϕ  = m.Aᵩ *(1 - exp(-sL/m.δᵩ ))*sin(2*T(π)*m.fᵩ *t - 2*T(π)*sL/m.λᵩ + m.γ)         + m.Cᵩ*sL
    sθ, cθ = sincos(θ); sϕ, cϕ = sincos(ϕ)
    SVector(cθ*cϕ, cθ*sϕ, sθ)
end

@inline function unit_tangent_and_dt(s::T, t::T, m::ThreeDimensionalFlagellum) where {T <: Number}
    sL   = s*m.L
    envθ = m.A_θ*(1 - exp(-sL/m.δ_θ))
    envϕ = m.Aᵩ *(1 - exp(-sL/m.δᵩ ))
    φθ   = 2*T(π)*m.f_θ*t - 2*T(π)*sL/m.λ_θ + m.γ + m.Δγ
    φϕ   = 2*T(π)*m.fᵩ *t - 2*T(π)*sL/m.λᵩ + m.γ
    θ    = envθ*sin(φθ) + m.C_θ*sL
    ϕ    = envϕ*sin(φϕ) + m.Cᵩ*sL
    θdot  = envθ * 2*T(π)*m.f_θ * cos(φθ)            # true ∂θ/∂t
    ϕ̇   = envϕ * 2*T(π)*m.fᵩ  * cos(φϕ)            # true ∂ϕ/∂t

    sθ, cθ = sincos(θ); sϕ, cϕ = sincos(ϕ)
    τ  = SVector(cθ*cϕ, cθ*sϕ, sθ)
    τdot = θdot * SVector(-sθ*cϕ, -sθ*sϕ, cθ) +         # ∂t̂/∂θ
         ϕ̇ * SVector(-cθ*sϕ,  cθ*cϕ, zero(T))      # ∂t̂/∂ϕ
    (τ, τdot)
end

# ===========================================================================
#  Standing-wave mode machinery (shared by the two standing-wave models)
#
#  One angle bank, equivalent to the complex form
#     Re[e^{iωt}θ_space + e^{-iωt}conj(θ_space)],  θ_space = Σ Aₙ e^{iφₙ} bₙ(s):
#       angle = C·s + 2 Σ Aₙ cos(ωt+φₙ) bₙ(s)
#       rate  =      −2ω Σ Aₙ sin(ωt+φₙ) bₙ(s)
#  with the clamped-base spatial modes bₙ(s) = sin((2n−1)π s/2), n = 1..4.
# ===========================================================================
@inline _modes(s::T) where {T<:Number} =
    SVector(sin(T(π)*s/2), sin(3*T(π)*s/2), sin(5*T(π)*s/2), sin(7*T(π)*s/2))

# @inline function _bank_angle(s::T, t::T, ω::T, C::T,
#                              A::SVector{4,T}, φ::SVector{4,T}) where {T<:Number}
#     b  = _modes(s)
#     ωt = ω*t
#     c  = SVector(cos(ωt+φ[1]), cos(ωt+φ[2]), cos(ωt+φ[3]), cos(ωt+φ[4]))
#     C*s + 2*sum(A .* c .* b)
# end


@inline function get_angle(s::T, t::T, ω::T, C::T,
                                R::SVector{4,T}, I::SVector{4,T}) where {T<:Number}
    b  = _modes(s)
    ωt = ω*t
    C*s + 2*(cos(ωt)*dot(R, b) - sin(ωt)*dot(I, b))
end

# @inline function get_angle_and_rate(s::T, t::T, ω::T, C::T,
#                                       A::SVector{4,T}, φ::SVector{4,T}) where {T<:Number}
#     b  = _modes(s)
#     ωt = ω*t
#     c  = SVector(cos(ωt+φ[1]), cos(ωt+φ[2]), cos(ωt+φ[3]), cos(ωt+φ[4]))
#     sn = SVector(sin(ωt+φ[1]), sin(ωt+φ[2]), sin(ωt+φ[3]), sin(ωt+φ[4]))
#     ( C*s + 2*sum(A .* c .* b),  -2*ω*sum(A .* sn .* b) )
# end

@inline function get_angle_and_rate(s::T, t::T, ω::T, C::T,
                                      R::SVector{4,T}, I::SVector{4,T}) where {T<:Number}
    b  = _modes(s)
    ωt = ω*t
    (C*s + 2*(cos(ωt)*dot(R, b) - sin(ωt)*dot(I, b)), -2ω*(sin(ωt)*dot(R, b) + cos(ωt)*dot(I, b)))
end
# ===========================================================================
#  4. PlanarStandingWaveFlagellum — single modal bank, in-plane (z ≡ 0).
# ===========================================================================
mutable struct PlanarStandingWaveFlagellum{T <: Number} <: FlagellumModel
    L::T
    ω::T
    C::T                       # static curvature
    R::SVector{4,T}            # real part of amplitudes
    I::SVector{4,T}            # imaginary part of amplitudes
end

PlanarStandingWaveFlagellum(L, ω, C, R, I) = PlanarStandingWaveFlagellum(
    Float64(L), Float64(ω), Float64(C), SVector{4,Float64}(R), SVector{4,Float64}(I)
)

(m::PlanarStandingWaveFlagellum)(s::Float64, t::Float64) = get_angle(s, t, m.ω, m.C, m.R, m.I)

@inline function unit_tangent(s::T, t::T, m::PlanarStandingWaveFlagellum) where {T<:Number}
    θ = get_angle(s, t, m.ω, m.C, m.R, m.I)
    SVector(cos(θ), sin(θ), zero(T))
end

@inline function unit_tangent_and_dt(s::T, t::T, m::PlanarStandingWaveFlagellum) where {T<:Number}
    θ, θdot = get_angle_and_rate(s, t, m.ω, m.C, m.R, m.I)
    (SVector(cos(θ), sin(θ), zero(T)), θdot * SVector(-sin(θ), cos(θ), zero(T)))
end

# ===========================================================================
#  5. ThreeDimensionalStandingWaveFlagellum — two independent modal banks.
#     Single beat frequency ω; helicity set by the phase offset between the
#     elevation (θ) and azimuth (ϕ) banks.  t̂ spherical ⇒ unit-norm.
# ===========================================================================
mutable struct ThreeDimensionalStandingWaveFlagellum{T <: Number} <: FlagellumModel
    L::T
    ω::T                       # single beat frequency
    C_θ::T                     # elevation static curvature
    R_θ::SVector{4,T}          # elevation mode real part
    I_θ::SVector{4,T}          # elevation mode imaginary part
    Cᵩ::T                      # azimuth static curvature
    Rᵩ::SVector{4,T}           # azimuth mode real part
    Iᵩ::SVector{4,T}           # azimuth mode imaginary part
end

ThreeDimensionalStandingWaveFlagellum(L, ω, C_θ, R_θ, I_θ, Cᵩ, Rᵩ, Iᵩ) = ThreeDimensionalStandingWaveFlagellum(
    Float64(L), Float64(ω), Float64(C_θ), SVector{4,Float64}(R_θ), SVector{4,Float64}(I_θ),
    Float64(Cᵩ), SVector{4,Float64}(Rᵩ), SVector{4,Float64}(Iᵩ)
)

@inline function unit_tangent(s::T, t::T, m::ThreeDimensionalStandingWaveFlagellum) where {T<:Number}
    θ = get_angle(s, t, m.ω, m.C_θ, m.R_θ, m.I_θ)
    ϕ = get_angle(s, t, m.ω, m.Cᵩ, m.Rᵩ, m.Iᵩ)
    sθ, cθ = sincos(θ); sϕ, cϕ = sincos(ϕ)
    SVector(cθ*cϕ, cθ*sϕ, sθ)
end

@inline function unit_tangent_and_dt(s::T, t::T, m::ThreeDimensionalStandingWaveFlagellum) where {T<:Number}
    θ, θdot = get_angle_and_rate(s, t, m.ω, m.C_θ, m.R_θ, m.I_θ)
    ϕ, ϕdot = get_angle_and_rate(s, t, m.ω, m.Cᵩ, m.Rᵩ, m.Iᵩ)
    sθ, cθ = sincos(θ); sϕ, cϕ = sincos(ϕ)
    τ  = SVector(cθ*cϕ, cθ*sϕ, sθ)
    τdot = θdot * SVector(-sθ*cϕ, -sθ*sϕ, cθ) +         # ∂t̂/∂θ
         ϕdot * SVector(-cθ*sϕ,  cθ*cϕ, zero(T))      # ∂t̂/∂ϕ
    (τ, τdot)
end
