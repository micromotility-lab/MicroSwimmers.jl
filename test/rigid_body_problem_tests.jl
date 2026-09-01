@testset "rigid_body_problem_tests.jl test set" begin

    function ellipsoid_swimmer(; a = 1.0, b = 1.0, c = 1.0, N = 213, Q = 879, orientation = rotation_matrix([0.0, 0.0, 1.0], 0))
        model = EllipsoidBody(a, b, c)
        f = Part(model, N, Q)
        MicroSwimmer([f],orientation=orientation)
    end
    
    #test that body shape parameter matches ellipsoid aspect ratio
    @testset "body_shape_parameter matches ellipsoid aspect ratio" begin
        test_a = [0.5,1.0,2.0,5.0]
        for a in test_a
            ms = ellipsoid_swimmer(a = a, b = 1.0, c = 1.0)
            B = body_shape_parameter(ms; eps = 1e-3)
            Breth = (a^2-1)/(a^2+1)
            @test abs(Breth - B) < 1e-2
        end
    end
    #test that body shape parameter is NaN for asymmetric body
    @testset "body_shape_parameter is NaN for asymmetric body" begin
        ms = ellipsoid_swimmer(a = 1.0, b = 2.0, c = 3.0)
        B = body_shape_parameter(ms; eps = 1e-3)
        @test isnan(B)
    end
    #test that body shape parameter is invariant under rotation of the body
    @testset "body_shape_parameter is invariant under rotation of the body" begin
        test_a = [0.5,2.0]
        test_phi = [ π/4, π/2, 7π/8]
        for a in test_a
            for phi in test_phi
                ms = ellipsoid_swimmer(a = a, b = 1.0, c = 1.0, orientation = rotation_matrix([0.0, 0.0, 1.0], phi))
                B = body_shape_parameter(ms; eps = 1e-3)
                Breth = (a^2-1)/(a^2+1)
                @test abs(Breth - B) < 1e-2
            end
        end
    end
end
