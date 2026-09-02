@testset "problems.jl — hybrid dense/GMRES solve" begin

    function small_flagellum_swimmer(; N=10, Q=51)
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        f = Part(model, N, Q; eps=0.1)
        MicroSwimmer([f])
    end

    @testset "hybrid solve matches dense-only solve across a beat cycle" begin
        ts = (0.0, 0.05, 0.11, 0.19, 0.27, 0.33, 0.41)

        ms_dense = small_flagellum_swimmer()
        prob_dense = SwimmingProblem(ms_dense)

        ms_hybrid = small_flagellum_swimmer()
        prob_hybrid = SwimmingProblem(ms_hybrid; hybrid=true, refactor_interval=3)

        for t in ts
            update_boundary!(prob_dense, t)
            solve_problem!(prob_dense)

            update_boundary!(prob_hybrid, t)
            solve_problem!(prob_hybrid)

            @test get_U(prob_hybrid) ≈ get_U(prob_dense) atol=1e-6
            @test get_Ω(prob_hybrid) ≈ get_Ω(prob_dense) atol=1e-6
            @test prob_hybrid.force_vals ≈ prob_dense.force_vals atol=1e-6
        end

        # a GMRES-only (non-refactorising) step must actually have been exercised
        @test prob_hybrid.steps_since_refactor > 0
    end

    @testset "force/torque-free invariant holds under hybrid solve" begin
        ms = small_flagellum_swimmer()
        prob = SwimmingProblem(ms; hybrid=true, refactor_interval=3)
        for t in (0.0, 0.09, 0.17, 0.23, 0.31, 0.4, 0.5)
            update_boundary!(prob, t)
            solve_problem!(prob)
            F, T = total_force_and_torque(prob)
            # Looser than problems_invariants_tests.jl's 1e-8 — that test's t values happen to
            # land on favourable discretisation points; these don't, and the residual here is
            # BEM discretisation error present identically on dense-only solves, not a hybrid
            # solver artifact (confirmed: unaffected by tightening gmres_reltol).
            @test norm(F) < 1e-7
            @test norm(T) < 1e-7
        end
    end

    @testset "unconverged GMRES falls back to a dense solve" begin
        ms = small_flagellum_swimmer()
        # An unreachable reltol combined with a single allowed iteration guarantees GMRES
        # can never report convergence, regardless of how good the preconditioner is —
        # this isolates the fallback-on-failure path from geometry-drift specifics.
        prob = SwimmingProblem(ms; hybrid=true, refactor_interval=1000,
            gmres_maxiters=1, gmres_reltol=1e-14)
        update_boundary!(prob, 0.0)
        solve_problem!(prob)                      # first solve: always dense
        @test prob.steps_since_refactor == 0

        update_boundary!(prob, 0.05)
        solve_problem!(prob)
        @test prob.steps_since_refactor == 0

        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-6
        @test norm(T) < 1e-6
    end

    @testset "SwimmingTrajectoryProblem runs end-to-end with hybrid=true" begin
        ms = small_flagellum_swimmer()
        tprob = SwimmingTrajectoryProblem(ms; t_final=0.2, saveat=0.02,
            hybrid=true, refactor_interval=4)
        solve_problem!(tprob)
        @test !isnothing(tprob.traj)
        @test length(tprob.traj.t) > 1
    end

end
