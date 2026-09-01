@testset "regression tests (from MicroSwimmersPlots/examples/getting_started.jl)" begin

    # These pin the solver's output for realistic, known-good configurations
    # taken directly from the package's example script. Unlike the analytic
    # Stokes'-law tests, there's no closed-form reference here — the point is
    # to catch any change (refactor, threading, sign convention, multi-part
    # assembly bug) that silently alters the numerical result for a case the
    # maintainer has already visually verified via MicroSwimmersPlots.
    #
    # Reference values below were captured by running this exact code against
    # the current `main` branch and are expected to remain stable under
    # non-behavioural changes (performance work, refactors, comments).
    #
    # eps = 0.1 was the problem-level default when these values were captured. It is passed
    # per part now that eps lives on the discretisation, so the references stay comparable.

    @testset "isolated flagellum" begin
        N, Q = 23, 127
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        f = Part(model, N, Q; eps=0.1)
        ms = MicroSwimmer([f])

        prob = SwimmingProblem(ms)
        solve_problem!(prob)

        U = get_U(prob)
        Ω = get_Ω(prob)
        F, T = total_force_and_torque(prob)

        @test norm(F) < 1e-8
        @test norm(T) < 1e-8

        U_ref = SVector(0.001792467908360207, -0.2237906168660196, 0.0)
        Ω_ref = SVector(0.0, 0.0, -0.4903399840842083)
        @test U ≈ U_ref atol=1e-6
        @test Ω ≈ Ω_ref atol=1e-6
    end

    @testset "sperm-like swimmer (body + flagellum, multi-part assembly)" begin
        a = b = c = 1.0
        N_body, Q_body = 213, 917
        N, Q = 23, 127

        body = Part(EllipsoidBody(a, b, c), N_body, Q_body; eps=0.1)
        f = Part(
            PlanarFlagellum(10.0, 0.0, 0.6, 0.5, π/2, 2π, 2π, 0.0),
            N, Q,
            eps=0.1,
            location=[1.0, 0.0, 0.0],
            orientation=rotation_matrix([1.0, 0.0, 0.0], 0.0)
        )
        ms = MicroSwimmer([body, f])

        sprob = SwimmingProblem(ms)
        solve_problem!(sprob)

        U = get_U(sprob)
        Ω = get_Ω(sprob)
        F, T = total_force_and_torque(sprob)

        @test norm(F) < 1e-6
        @test norm(T) < 1e-6

        U_ref = SVector(-0.24650680886215123, -1.458308536376894, -5.4033768937554515e-5)
        Ω_ref = SVector(2.933933320197317e-5, -8.475110965006673e-6, -0.5748740208904173)
        @test U ≈ U_ref atol=1e-4
        @test Ω ≈ Ω_ref atol=1e-4
    end

end
