module MicroSwimmers

using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using LazyArrays
# using Dierckx
using GeometryBasics: Mesh, Point3, Vec3, Vec3f, Mat3, coordinates, faces
using DiffEqBase, OrdinaryDiffEq, LinearSolve
# using Meshing

include("maths.jl")
include("geometry.jl")
include("numerics.jl")
include("boundary.jl")
include("nearest.jl")
include("cell_body.jl")
include("flagellum_models.jl")
include("flagellum.jl")
include("swimmers.jl")
include("problems.jl")
include("exports.jl")

end # module