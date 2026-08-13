using MicroSwimmers
using Test
using Random
using LinearAlgebra
using StaticArrays
using Aqua

@testset "MicroSwimmers.jl" begin
    include("maths_tests.jl")
    include("frame_tests.jl")
    include("geometry_tests.jl")
    include("numerics_tests.jl")
    include("discretisations_tests.jl")
    include("problems_stokes_law_tests.jl")
    include("problems_invariants_tests.jl")
    include("trajectories_tests.jl")
    include("regression_tests.jl")
    include("quality_tests.jl")
end
