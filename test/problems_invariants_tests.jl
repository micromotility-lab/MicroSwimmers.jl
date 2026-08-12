@testset "problems.jl — physical invariants" begin

    function small_flagellum_swimmer(; N=10, Q=51)
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        f = Part(model, N, Q)
        MicroSwimmer([f])
    end

    @testset "SwimmingProblem is force- and torque-free" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms; eps=0.1)
        solve_problem!(prob)

        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "force/torque-free holds at several phases of the beat" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms; eps=0.1)
        for t in (0.0, 0.13, 0.37, 0.81)
            update_boundary!(prob, t)
            solve_problem!(prob)
            F, T = total_force_and_torque(prob)
            @test norm(F) < 1e-8
            @test norm(T) < 1e-8
        end
    end

    @testset "stresslet_tensor is symmetric and traceless" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms; eps=0.1)
        solve_problem!(prob)
        S = stresslet_tensor(prob)
        @test S ≈ S' atol=1e-10
        @test abs(tr(S)) < 1e-10
    end

    @testset "total_power is non-negative for a dragged sphere" begin
        a = 1.0
        body = Part(EllipsoidBody(a, a, a), 30, 130)
        ms = MicroSwimmer([body])
        prob = ResistanceProblem(ms; eps=0.06)
        add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
        solve_problem!(prob)
        @test MicroSwimmers.total_power(prob) > 0
    end

end
