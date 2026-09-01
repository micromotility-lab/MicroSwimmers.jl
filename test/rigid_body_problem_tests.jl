@testset "rigid_body_problem_tests.jl test set" begin

    function ellipsoid_swimmer(; a = 1.0, b = 1.0, c = 1.0, N = 213, Q = 879, orientation = rotation_matrix([0.0, 0.0, 1.0], 0))
        model = EllipsoidBody(a, b, c)
        # eps = 1e-5: this is a Bretherton/Jeffery shape test, so the regularisation is kept
        # far below the node spacing to approach the singular-Stokeslet limit.
        f = Part(model, N, Q; eps = 1e-5)
        MicroSwimmer([f],orientation=orientation)
    end
    
    #test that body shape parameter matches ellipsoid aspect ratio
    @testset "body_shape_parameter matches ellipsoid aspect ratio" begin
        test_a = [0.5,1.0,2.0,5.0]
        for a in test_a
            ms = ellipsoid_swimmer(a = a, b = 1.0, c = 1.0)
            B = body_shape_parameter(ms)
            Breth = (a^2-1)/(a^2+1)
            @test abs(Breth - B) < 1e-2
        end
    end
    #test that body shape parameter is NaN for asymmetric body
    @testset "body_shape_parameter is NaN for asymmetric body" begin
        ms = ellipsoid_swimmer(a = 1.0, b = 2.0, c = 3.0)
        B = body_shape_parameter(ms)
        @test isnan(B)
    end
    #test that body shape parameter is invariant under rotation of the body
    @testset "body_shape_parameter is invariant under rotation of the body" begin
        test_a = [0.5,2.0]
        test_phi = [ π/4, π/2, 7π/8]
        for a in test_a
            for phi in test_phi
                ms = ellipsoid_swimmer(a = a, b = 1.0, c = 1.0, orientation = rotation_matrix([0.0, 0.0, 1.0], phi))
                B = body_shape_parameter(ms)
                Breth = (a^2-1)/(a^2+1)
                @test abs(Breth - B) < 1e-2
            end
        end
    end

    #test that rigid body problem orientation_ode returns expected values for ellipsoid in shear flow
    @testset "RigidBodyProblem returns expected values for ellipsoid in shear flow" begin
        test_a = [0.5,1.0,2.0]
        test_theta = [pi/4, 3*pi/4] #need off axis initial orientation to rotate in shear flow
        for a in test_a
            for theta_0 in test_theta
                initial_orientation = @SMatrix [
                    sin(theta_0)    -cos(theta_0)   0;
                    0      0     1;
                    cos(theta_0)   sin(theta_0)  0
                ]
                ms = ellipsoid_swimmer(a = a, b = 1.0, c = 1.0, orientation = initial_orientation)
                Omega_shear = @SVector [0.0, 0.0, -0.5]
                strain_shear = @SMatrix [0.0 0.5 0.0; 0.5 0.0 0.0; 0.0 0.0 0.0]
                sol, x, b1, b2 = RigidBodyProblem(ms; Omega_shear = Omega_shear, strain_shear = strain_shear)
                #check that Jeffery conserved quantity K is conserved
                K = tan(theta_0)
                theta_num = acos.(clamp.(getindex.(b1, 3), -1.0, 1.0))
                phi_num   = atan.(getindex.(b1, 2), getindex.(b1, 1))
                orbit_const = tan.(theta_num) .* sqrt.(cos.(phi_num).^2 .+ a^2 .* sin.(phi_num).^2)
                @test maximum(abs,orbit_const .- K) < 1e-2
            end
        end
    end
    end

