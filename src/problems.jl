abstract type Problem end

mutable struct SwimmingProblem{T <: Number} <: Problem
    swimmer::Swimmer
    points::Discretisation  # IS THIS NECESSARY?
    eps::T   # regularisation parameter
    mu::T    # viscosity

    lin_prob::LinearProblem
    force_vals::Union{Nothing, Vector{T}} #  keep track of the values of the forces at the force points
    wall::Bool
end

function SwimmingProblem(
    S::Swimmer;
    x0::SVector{3,T}=SVector(0., 0., 0.), 
    B::SMatrix{3,3,T}=I3,
    eps=0.01,
    mu=1.
) where {T <: Number} 

    N = length(S.points.force_pts)
    points = NearestDiscretisation(
        zeros(T, 3, S.points.N),
        zeros(T, 3, S.points.Q),
        S.points.nearest
    )
    sp = SwimmingProblem(
        S, points, eps, mu,
        LinearProblem(zeros(T, N + 6, N + 6), zeros(T, N + 6)),
        nothing,
        false
    )
    move_boundary!(sp, x0, B, zero(T))
    sp
end

function check_solved!(prob)
    if isnothing(prob.force_vals)
        solve_problem!(prob)
    end
end

function move_boundary!(sp::SwimmingProblem, x0::SVector{3,T}, B::SMatrix{3,3,T}, t::T) where {T <: Number}
    move_boundary!(sp.swimmer, x0, B, t)

    @unpack force_pts, quad_pts, velocity = sp.swimmer.points
    @views begin
        sp.points.force_pts .= x0 .+ B * sp.swimmer.points.force_pts
        sp.points.velocity  .= B * sp.swimmer.points.velocity
        sp.points.quad_pts  .= x0 .+ B * sp.swimmer.points.quad_pts
    end
end

function move_boundary!(sp::SwimmingProblem, x0::SVector{3,T}, b1::SVector{3,T}, b2::SVector{3, T}, t::T) where {T <: Number}
    B = hcat(b1, b2, cross(b1, b2))
    move_boundary!(sp, x0, B, t)
end

function solve_problem!(sp::SwimmingProblem)
    swimming_matrix!(
        sp.lin_prob.A,
        sp.swimmer.points.location,
        sp.points.force_pts,
        sp.points.quad_pts,
        sp.swimmer.points.nearest,
        sp.eps,
        μ=sp.mu
    )

    @views sp.lin_prob.b[1:end-6] .= reshape(sp.points.velocity, :)
    sp.force_vals = solve(sp.lin_prob, MKLLUFactorization())
end

mutable struct DynamicSwimmingProblem <: Problem
    swimming_problem::SwimmingProblem
    ode_prob::ODEProblem
    sol::Union{Nothing, ODESolution}
end

function DynamicSwimmingProblem(
    S::Swimmer;
    x0::SVector{3,T}=SVector(0., 0., 0.), 
    B::SMatrix{3,3,T}=I3, 
    t_final::T=20.,
    saveat::T=0.05,
    eps=0.01,
    mu=1.
) where {T <: Number} 
    sp = SwimmingProblem(S; x0=x0, B=B, eps=eps, mu=mu)

    x0_0 = zeros(3)
    b1_0 = [1., 0., 0.]
    b2_0 = [0., 1., 0.]

    function rhs!(dX, X, p, t)
        # X[4:6] .= X[4:6] ./ norm(X[4:6]) # prevent drift in basis normalisation
        # X[7:9] .= X[7:9] ./ norm(X[7:9])
        move_boundary!(sp, SVector{3}(X[1:3]), SVector{3}(X[4:6]), SVector{3}(X[7:9]), t)
        solve_problem!(sp)

        @views begin
            dX[1:3] .= sp.force_vals[end-5:end-3]
            dX[4:6] .= cross(sp.force_vals[end-2:end], X[4:6])
            dX[7:9] .= cross(sp.force_vals[end-2:end], X[7:9])
        end
    end

    DynamicSwimmingProblem(
        sp,
        ODEProblem(rhs!, [x0_0; b1_0; b2_0], (0.0, t_final), saveat=saveat),
        nothing
    )
end


function solve_problem!(dsp::DynamicSwimmingProblem; method=Tsit5())
    dsp.sol = solve(dsp.ode_prob, method)
end

