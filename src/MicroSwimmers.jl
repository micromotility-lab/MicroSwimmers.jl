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
# using Meshing
using FastGaussQuadrature

include("maths.jl")
include("geometry.jl")
include("numerics.jl")
include("boundary.jl")
include("nearest.jl")
include("cell_body.jl")
include("flagellum_models.jl")
include("swimmers.jl")
include("flagella.jl")
include("problems.jl")
include("trajectories.jl")
include("forces_and_torques.jl")
include("fluid.jl")
include("exports.jl")

end # module