struct Flagellum{M <: FlagellumModel, T <: Number} <: FluidBoundary
    model::M
    N::Int
    Q::Int
    config::Configuration{T}
    ϵ::T
end

function Flagellum(; 
    model=PlanarFlagellum(1., 0., 0.3, 0.15, 2π, -2π, 2π),
    N=23, 
    location=SVector(0., 0., 0.),
    orientation=I3,
    Q=127, 
    ϵ=0.01
)
    f = Flagellum(
        model,
        N,
        Q,
        Configuration(location, orientation, zeros(3,N), zeros(3,N), zeros(3,Q), zeros(Int, Q)),
        ϵ
    )
    update_boundary!(f, 0.)
    f.config.nearest .= nearest_neighbour(f.config.force_pts, f.config.quad_pts)
    f
end


function update_boundary!(fl::Flagellum, t::T) where {T <: Number}
    @unpack force_pts, quad_pts, velocity = fl.config
    fl.model(force_pts, velocity, t)
    fl.model(quad_pts, t)
end


 





# mutable struct PlanarFlagellum{T <: Number} <: Flagellum
#     state::FlagellumState{T}

#     # Waveform parameters
#     C::T    # intrinsic curvature
#     R₀::T   # amplitude at s=0
#     R₁::T   # amplitude of envelope 
#     k::T    # wavelength of envelope
#     ϕ::T    # wavelength of phase
#     ω::T    # frequency

#     # Preallocated integrand buffers
#     cosθ::Vector{T}
#     sinθ::Vector{T}
#     ωθ₁sinθsin::Vector{T}
#     ωθ₁cosθsin::Vector{T}
# end







# function PlanarFlagellum(
#     N::Int, Q::Int, ϵ::T, C::T, R₀::T, R₁::T, k::T, L::T, ϕ::T, ω::T;
#     N_int::Int=150, 
#     location=SVector{3,T}(zero(T),zero(T),zero(T)),
#     orientation=SMatrix{3,3,T}(I),
# ) where {T<:Number}

#     pf = PlanarFlagellum(
#         FlagellumState{T}(
#             N, Q, L,
#             collect(LinRange{T}(zero(T), one(T), N_int)),
#             collect(LinRange{T}(zero(T), one(T), N)),
#             collect(LinRange{T}(zero(T), one(T), Q)),
#             Configuration(location, orientation, zeros(T,3,N), zeros(T,3,Q), zeros(T,3,N)),
#             zeros(Int, Q),
#             ϵ
#         ),
#         C, R₀, R₁, k, ϕ, ω,
#         zeros(N_int),
#         zeros(N_int),
#         zeros(N_int),   
#         zeros(N_int)
#     )
#     update!(pf, zero(T))
#     @unpack force_pts, quad_pts = pf.state.config
#     pf.state.nearest .= nearest_neighbour(force_pts, quad_pts)
#     pf
# end


# "Update the reference flagellum configuration to time t"
# function update!(pf::PlanarFlagellum, t::Float64)
#     @unpack L, s_int, s_force, s_quad = pf.state
#     @unpack C, R₀, R₁, k, ϕ, ω, cosθ, sinθ, ωθ₁sin = pf

#     @inbounds for i in eachindex(s_int)
#         s = s_int[i]
#         θ₁ = (R₀ + R₁*sin(k*s))
#         θ = C*s + θ₁*cos(ω*t + ϕ*s)
#         cosθ[i] = cos(θ)
#         sinθ[i] = sin(θ)
#         ωθ₁sin = ω * θ₁*sin(ω*t + ϕ*s)
#         ωθ₁sinθsin[i] = ωθ₁sin * sinθ[i]
#         ωθ₁cosθsin[i] = ωθ₁sin * cosθ[i]
#     end



#     spl_x = Spline1D(L*s_force, cosθ)
#     spl_y = Spline1D(L*s_force, sinθ)
#     spl_xdot = Spline1D(L*s_force,  ωθ₁sin .* sinθ)
#     spl_ydot = Spline1D(L*s_force, -ωθ₁sin .* cosθ)

#     @unpack force_pts, quad_pts, velocity = pf.config

#     cumulative_trapz!(force_pts[1,:], L*s_force, cosθ)
#     cumulative_trapz!()
#     ds = s_force[1]
#     @inbounds for i in eachindex(pf.s_force)
#         force_pts[1,i] = integrate(spl_x, 0., L*s_force[i])
#         force_pts[2,i] = integrate(spl_y, 0., L*s_force[i])
        
#         velocity[1,i]  = integrate(spl_xdot, 0., L*s_force[i])
#         velocity[2,i]  = integrate(spl_ydot, 0., L*s_force[i])
#     end

#     @inbounds for i in eachindex(pf.s_force)
#         force_pts[1,i] = integrate(spl_x, 0., L*s_force[i])
#         force_pts[2,i] = integrate(spl_y, 0., L*s_force[i])
        
#         velocity[1,i]  = integrate(spl_xdot, 0., L*s_force[i])
#         velocity[2,i]  = integrate(spl_ydot, 0., L*s_force[i])
#     end

#     @inbounds for i in eachindex(pf.s_quad)
#         quad_pts[1,i] = integrate(spl_x, 0., pf.L*pf.s_quad[i])
#         quad_pts[2,i] = integrate(spl_y, 0., pf.L*pf.s_quad[i])
#     end
# end

