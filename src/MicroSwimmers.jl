module MicroSwimmers

using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using LazyArrays
# using Dierckx
using DiffEqBase, OrdinaryDiffEq
using LinearSolve
using LsqFit
using Statistics
# using Meshing
using FastGaussQuadrature
using DSP
using ForwardDiff

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
include("problem2.jl")
include("forces_and_torques.jl")
include("fluid.jl")
include("exports.jl")

end # module