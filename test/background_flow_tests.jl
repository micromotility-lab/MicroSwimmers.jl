@testset "background_flow.jl — imposed background flow" begin

    # A rigid sphere is the cleanest check: in any linear flow a free sphere translates with
    # the local flow velocity and spins at half the vorticity (Faxén), whatever its resolution.
    sphere_swimmer(; a=1.0, N=60, Q=250, eps=0.04) =
        MicroSwimmer([Part(EllipsoidBody(a, a, a), N, Q; eps=eps)])

    function place!(ms, x0)
        ms.frame = Frame(SVector{3}(x0), SMatrix{3,3}(1.0I))
        ms
    end

    @testset "constructors" begin
        f = LinearFlow([0 1 0; 0 0 0; 0 0 0], [1, 2, 3])
        @test f isa LinearFlow{Float64}
        @test f(SVector(1.0, 2.0, 3.0)) == SVector(3.0, 2.0, 3.0)
        @test ShearFlow(2.0)(SVector(0.0, 1.5, 0.0)) == SVector(3.0, 0.0, 0.0)
        @test UniformFlow([1.0, 0.0, 0.0])(SVector(5.0, 5.0, 5.0)) == SVector(1.0, 0.0, 0.0)
        @test LinearFlow(U0=[0.0, 1.0, 0.0]).G == zero(SMatrix{3,3,Float64,9})
        f.U0 = SVector(0.0, 0.0, 0.0)   # mutable, so a slider or optimiser can drive it
        @test f.U0 == zero(SVector{3,Float64})
    end

    @testset "NoFlow leaves the solve bit-identical" begin
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        p1 = SwimmingProblem(MicroSwimmer([Part(model, 10, 51; eps=0.1)]))
        p2 = SwimmingProblem(MicroSwimmer([Part(model, 10, 51; eps=0.1)]); flow=NoFlow())
        update_boundary!(p1, 0.13); update_boundary!(p2, 0.13)
        @test solve_problem!(p1) == solve_problem!(p2)
    end

    @testset "free sphere is carried by a uniform flow" begin
        U0 = SVector(0.3, -1.0, 0.5)
        prob = SwimmingProblem(sphere_swimmer(); flow=UniformFlow(U0))
        solve_problem!(prob)
        @test get_U(prob) ≈ U0 atol=1e-8
        @test norm(get_Ω(prob)) < 1e-8
        @test maximum(norm, get_forces(prob)) < 1e-8
    end

    @testset "free sphere in simple shear: U = G x0, Ω = ω/2" begin
        γ  = 2.0
        x0 = SVector(0.4, 1.5, -0.2)
        prob = SwimmingProblem(place!(sphere_swimmer(), x0); flow=ShearFlow(γ))
        solve_problem!(prob)
        @test get_U(prob) ≈ SVector(γ * x0[2], 0.0, 0.0) rtol=1e-2
        @test get_Ω(prob) ≈ SVector(0.0, 0.0, -γ / 2) rtol=1e-2
    end

    @testset "free sphere in pure strain neither moves nor spins" begin
        prob = SwimmingProblem(sphere_swimmer(); flow=LinearFlow([0 1 0; 1 0 0; 0 0 0]))
        solve_problem!(prob)
        @test norm(get_U(prob)) < 1e-2
        @test norm(get_Ω(prob)) < 1e-2
        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "held sphere in a uniform flow: Galilean invariance and Stokes drag" begin
        a, μ = 1.0, 1.0
        U0 = SVector(1.0, 0.0, 0.0)

        ms_flow = sphere_swimmer(a=a)
        p_flow  = ResistanceProblem(ms_flow; mu=μ, flow=UniformFlow(U0))
        solve_problem!(p_flow)

        # at rest in a flow U0 is the same problem as moving at -U0 through still fluid
        ms_move = sphere_swimmer(a=a)
        add_rigid_body_motion!(ms_move.parts[1], -U0, zero(SVector{3,Float64}))
        p_move  = ResistanceProblem(ms_move; mu=μ)
        solve_problem!(p_move)

        @test p_flow.force_vals ≈ p_move.force_vals rtol=1e-12
        F, _ = total_force_and_torque(p_flow)
        @test abs(norm(F) - 6π * μ * a) / (6π * μ * a) < 0.1
    end

    @testset "hybrid GMRES agrees with dense under a background flow" begin
        model = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
        flow  = LinearFlow([0 0.7 0; 0 0 0; 0 0 0], [0.1, 0.0, 0.2])
        pd = SwimmingProblem(MicroSwimmer([Part(model, 10, 51; eps=0.1)]); flow=flow)
        ph = SwimmingProblem(MicroSwimmer([Part(model, 10, 51; eps=0.1)]); flow=flow,
                             hybrid=true, refactor_interval=3)
        for t in (0.0, 0.05, 0.11, 0.19, 0.27)
            update_boundary!(pd, t); solve_problem!(pd)
            update_boundary!(ph, t); solve_problem!(ph)
            @test ph.force_vals ≈ pd.force_vals atol=1e-6
        end
        @test ph.steps_since_refactor > 0
    end

    @testset "sphere trajectory in shear is advected and rotated" begin
        γ  = 1.0
        x0 = SVector(0.0, 2.0, 0.0)
        ms = place!(sphere_swimmer(N=40, Q=170, eps=0.05), x0)
        traj = SwimmingTrajectoryProblem(ms; t_final=0.5, saveat=0.25, flow=ShearFlow(γ))
        solve_problem!(traj)
        tr_ = traj.traj
        t   = tr_.t[end]
        @test tr_.x[end] ≈ x0 + SVector(γ * x0[2] * t, 0.0, 0.0) rtol=1e-2
        θ = -γ / 2 * t
        @test tr_.b1[end] ≈ SVector(cos(θ), sin(θ), 0.0) atol=1e-2
    end

    @testset "FluidVelocity: disturbance by default, total on request" begin
        flow = ShearFlow(1.0)
        prob = SwimmingProblem(place!(sphere_swimmer(N=40, Q=170, eps=0.05), (0.0, 1.0, 0.0));
                               flow=flow)
        solve_problem!(prob)
        x = SVector(3.0, 2.5, 0.5)
        @test FluidVelocity(prob; total=true)(x) ≈ FluidVelocity(prob)(x) + flow(x)
        # no-slip: the total field matches the rigid-body motion on the surface
        y  = prob.disc.force_pts[1]
        ub = get_U(prob) + cross(get_Ω(prob), y - prob.microswimmer.frame.location)
        @test FluidVelocity(prob; total=true)(y) ≈ ub atol=5e-2
    end

    @testset "warnings" begin
        @test_logs (:warn, r"vanish on z = 0") ResistanceProblem(sphere_swimmer(N=20, Q=90);
            wall=true, flow=UniformFlow([1.0, 0.0, 0.0]))
        @test_logs (:warn, r"incompressible") SwimmingProblem(sphere_swimmer(N=20, Q=90);
            flow=LinearFlow([1 0 0; 0 0 0; 0 0 0]))
        @test_logs ResistanceProblem(place!(sphere_swimmer(N=20, Q=90), (0, 0, 2.0));
            wall=true, flow=LinearFlow([0 0 1; 0 0 0; 0 0 0]))
    end
end
