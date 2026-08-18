module MicroSwimmers

using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using LazyArrays
# using Dierckx
using DiffEqBase, OrdinaryDiffEq
using LinearSolve
using SciMLOperators
using LsqFit
using Statistics
# using Meshing
using FastGaussQuadrature
using DSP
using ForwardDiff
using Polyester

function __init__()
    # MKL's OpenMP runtime spin-waits (default 200ms) after a parallel region
    # (e.g. MKLLUFactorization's getrf inside solve!) instead of sleeping.
    # SwimmingTrajectoryProblem alternates threaded assembly (Polyester) and
    # solve! every timestep, so those spinning MKL threads can steal cores
    # from the next assembly step. Must be set before the first factorisation
    # runs, so this belongs in __init__, not a call site. Respects an
    # explicit user override set before `using MicroSwimmers`.
    haskey(ENV, "KMP_BLOCKTIME") || (ENV["KMP_BLOCKTIME"] = "0")
end

include("maths.jl")
include("geometry.jl")
include("frame.jl")
include("discretisations.jl")
include("numerics.jl")
include("boundary.jl")
include("cell_body_models.jl")
include("implicit_body.jl")
include("flagellum_models.jl")
include("flagellum_accessories.jl")
include("microswimmer.jl")
include("trajectories.jl")
include("problems.jl")
include("forces_and_torques.jl")
include("fluid.jl")
include("exports.jl")
include("rigid_body_problem.jl")

end # module