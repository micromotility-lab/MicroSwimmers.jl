module MicroSwimmers

using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using LazyArrays
# using Dierckx
using DiffEqBase, OrdinaryDiffEq, LinearSolve
# using Meshing

include("maths.jl")
include("geometry.jl")
include("numerics.jl")
include("boundary.jl")
include("nearest.jl")
include("cell_body.jl")
include("flagellum_models.jl")
include("swimmers.jl")
include("problems.jl")
include("exports.jl")

end # module