abstract type Problem end

mutable struct SwimmingProblem{T <: Number} <: Problem
    swimmer::Swimmer
    config::Configuration

    lin_prob::LinearProblem
    force_vals::Vector{T} # keep track of the values of the forces at the force points
    wall::Bool
end

function SwimmingProblem(S::Swimmer, x0::SVector{3,T}=SVector(0., 0., 0.), B::SMatrix{3,3,T}=I3) where {T <: Number} 
    N = length(S.config.force_pts)
    sp = SwimmingProblem(
        S,
        Configuration(
            SVector(zero(T), zero(T), zero(T)),
            I3,
            zeros(T, size(S.config.force_pts)),
            zeros(T, size(S.config.velocity)),
            zeros(T, size(S.config.quad_pts)),
            zeros(Int, size(S.config.quad_pts,2))
        ),
        LinearProblem(zeros(T, N + 6, N + 6), zeros(T, N + 6)),
        zeros(T, N + 6),
        false
    )
    move!(sp, x0, B, zero(T))
    sp
end

function move!(sp::SwimmingProblem, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::T) where {T <: Number}
    move!(sp.swimmer, x0, B, t)

    @unpack force_pts, quad_pts, velocity = sp.swimmer.config
    @views begin
        sp.config.force_pts .= x0 .+ B * sp.swimmer.config.force_pts
        sp.config.velocity  .= B * sp.swimmer.config.velocity
        sp.config.quad_pts  .= x0 .+ B * sp.swimmer.config.quad_pts
    end
end

function move!(sp::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3, T}, t::T) where {T <: Number}
    B = hcat(b1, b2, cross(b1, b2))
    move!(sp, x0, B, t)
end

function solve_problem!(sp::SwimmingProblem; μ=1.)
    swimming_matrix!(
        sp.lin_prob.A,
        sp.swimmer.config.location,
        sp.config.force_pts,
        sp.config.quad_pts,
        sp.swimmer.config.nearest,
        sp.swimmer.ϵ
    )
    # @unpack lin_prob, config, force_vals = sp
    # @unpack nearest, ϵ = sp.swimmer
    # x0 = sp.swimmer.config.location

    # swimming_matrix!(lin_prob.A, x0, config.force_pts, config.quad_pts, nearest, ϵ; μ=μ)
    sp.lin_prob.b[1:end-6] .= reshape(sp.config.velocity, :)
    sp.force_vals .= solve(sp.lin_prob, MKLLUFactorization())
end

mutable struct DynamicSwimmingProblem <: Problem
    swimming_problem::SwimmingProblem
    ode_prob::ODEProblem
    sol::Union{Nothing, ODESolution}
end

function DynamicSwimmingProblem(
    S::Swimmer, 
    x0::SVector{3,T}=SVector(0., 0., 0.), 
    B::SMatrix{3,3,T}=I3, 
    t_final::T=20.,
    saveat::T=0.05,    
) where {T <: Number} 
    sp = SwimmingProblem(S, x0, B)

    x0_0 = zeros(3)
    b1_0 = [1., 0., 0.]
    b2_0 = [0., 1., 0.]

    function rhs!(dX, X, p, t)
        X[4:6] .= X[4:6] ./ norm(X[4:6]) # prevent drift in basis normalisation
        X[7:9] .= X[7:9] ./ norm(X[7:9])
        move!(sp, SVector{3}(X[1:3]), SVector{3}(X[4:6]), SVector{3}(X[7:9]), t)
        solve_problem!(sp)

        @views begin
            dX[1:3] .= sp.force_vals[end-5:end-3]
            dX[4:6] .= cross(sp.force_vals[end-2:end], X[4:6])
            dX[7:9] .= cross(sp.force_vals[end-2:end], X[7:9])
        end
    end

    @show rhs!(zeros(9), [x0_0; b1_0; b2_0], nothing, 0.0)

    DynamicSwimmingProblem(
        sp,
        ODEProblem(rhs!, [x0_0; b1_0; b2_0], (0.0, t_final), saveat=saveat),
        nothing
    )
end


function solve_problem!(dsp::DynamicSwimmingProblem; method=Tsit5())
    dsp.sol = solve(dsp.ode_prob, method)
end

