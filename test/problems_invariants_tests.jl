@testset "problems.jl — physical invariants" begin

    function small_flagellum_swimmer(; N=10, Q=51)
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        f = Part(model, N, Q; eps=0.1)
        MicroSwimmer([f])
    end

    @testset "SwimmingProblem is force- and torque-free" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms)
        solve_problem!(prob)

        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "force/torque-free holds at several phases of the beat" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms)
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
        prob = SwimmingProblem(ms)
        solve_problem!(prob)
        S = stresslet_tensor(prob)
        @test S ≈ S' atol=1e-10
        @test abs(tr(S)) < 1e-10
    end

    @testset "total_power is non-negative for a dragged sphere" begin
        a = 1.0
        body = Part(EllipsoidBody(a, a, a), 30, 130; eps=0.06)
        ms = MicroSwimmer([body])
        prob = ResistanceProblem(ms)
        add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
        solve_problem!(prob)
        @test MicroSwimmers.total_power(prob) > 0
    end

    @testset "total_power equals F·U for a rigidly translating body" begin
        # A rigidly translating body has u = U at every node, so the rate of working is
        #   P = Σ_q ⟨f_{nearest[q]}, U⟩ = ⟨Σ_q f_{nearest[q]}, U⟩ = ⟨total_force, U⟩
        # exactly — the same quadrature-point sum that total_force uses. Summing over
        # force points instead drops the patch multiplicity and breaks this by ~Q/N,
        # which is how the pre-7b3326f version of total_power went unnoticed.
        a    = 1.0
        body = Part(EllipsoidBody(a, a, a), 40, 170; eps=0.05)
        ms   = MicroSwimmer([body])
        prob = ResistanceProblem(ms)
        U    = SVector(0.3, -0.2, 0.5)
        add_rigid_body_motion!(body, U, zero(SVector{3,Float64}))
        solve_problem!(prob)

        F, _ = total_force_and_torque(prob)
        @test total_power(prob) ≈ dot(F, U) rtol=1e-10

        # the identity is sharp: the force-point sum differs by the patch multiplicity,
        # so this test genuinely discriminates between the two conventions
        forces = get_forces(prob)
        @test !isapprox(sum(dot(f, U) for f in forces), dot(F, U); rtol=0.5)
    end

end
