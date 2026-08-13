@testset "geometry.jl" begin

    @testset "fibonacci_ellipsoid" begin
        for (a, b, c, N) in [(1.0, 1.0, 1.0, 50), (2.0, 1.0, 0.5, 77)]
            pts = fibonacci_ellipsoid(a, b, c, N)
            @test length(pts) == N

            # every point lies exactly on the ellipsoid surface
            for p in pts
                @test (p[1]/a)^2 + (p[2]/b)^2 + (p[3]/c)^2 ≈ 1.0 atol=1e-10
            end
        end

        # sphere case: every point is equidistant from the centre
        a = 1.7
        pts = fibonacci_ellipsoid(a, a, a, 60)
        r = norm.(pts)
        @test all(x -> isapprox(x, a; atol=1e-10), r)

        # roughly centred coverage (Fibonacci lattice is symmetric about the origin)
        centroid = sum(pts) / length(pts)
        @test norm(centroid) < 0.1 * a
    end

    @testset "is_inside_ellipsoid" begin
        center = [0.0, 0.0, 0.0]
        radii  = [1.0, 2.0, 3.0]

        @test is_inside_ellipsoid([0.0, 0.0, 0.0], center, radii)        # centre
        @test is_inside_ellipsoid([0.5, 0.5, 0.5], center, radii)        # strictly inside
        @test !is_inside_ellipsoid([2.0, 2.0, 2.0], center, radii)       # strictly outside
        @test is_inside_ellipsoid([1.0, 0.0, 0.0], center, radii)        # on the surface (within default tol)

        # translated ellipsoid
        c2 = [5.0, 0.0, 0.0]
        @test is_inside_ellipsoid([5.0, 0.0, 0.0], c2, radii)
        @test !is_inside_ellipsoid([0.0, 0.0, 0.0], c2, radii)
    end

end
