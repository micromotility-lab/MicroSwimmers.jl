@testset "discretisations.jl" begin

    @testset "nearest_neighbour (Vector{SVector} method)" begin
        force_pts = [SVector(0.0, 0.0, 0.0), SVector(10.0, 0.0, 0.0), SVector(0.0, 10.0, 0.0)]
        quad_pts  = [SVector(0.1, 0.0, 0.0), SVector(9.9, 0.1, 0.0), SVector(0.0, 9.8, 0.0), SVector(5.0, 5.0, 0.0)]

        nearest = MicroSwimmers.nearest_neighbour(force_pts, quad_pts)
        @test nearest == [1, 2, 3, 1]     # last point equidistant-ish from 1 and 3; force pt 1 is actually closest

        # exact coincidence: nearest of a force point to itself is itself
        nearest2 = MicroSwimmers.nearest_neighbour(force_pts, force_pts)
        @test nearest2 == [1, 2, 3]
    end

    @testset "nearest_neighbour! mutates disc.nearest" begin
        disc = NearestDiscretisation(2, 3)
        disc.force_pts .= [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)]
        disc.quad_pts  .= [SVector(0.1, 0.0, 0.0), SVector(0.4, 0.0, 0.0), SVector(0.9, 0.0, 0.0)]

        MicroSwimmers.nearest_neighbour!(disc)
        @test disc.nearest == [1, 1, 2]
    end

    @testset "tie-breaking picks the first index" begin
        force_pts = [SVector(-1.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)]
        quad_pts  = [SVector(0.0, 0.0, 0.0)]        # exactly equidistant

        nearest = MicroSwimmers.nearest_neighbour(force_pts, quad_pts)
        @test nearest == [1]
    end

end
