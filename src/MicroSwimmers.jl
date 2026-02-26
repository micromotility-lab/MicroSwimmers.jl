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

include("maths.jl")
include("geometry.jl")
include("discretisations.jl")
include("numerics.jl")
include("boundary.jl")
include("cell_body_models.jl")
include("cell_body_constructors.jl")
include("flagellum_models.jl")
include("flagellum_constructors.jl")
include("flagellate.jl")
include("trajectories.jl")
include("problems.jl")
include("forces_and_torques.jl")
include("fluid.jl")
include("exports.jl")

end # module