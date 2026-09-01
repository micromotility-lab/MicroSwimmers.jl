@testset "Nystrom (collocation) discretisation" begin

    # There is no make_discretisation returning a NystromDiscretisation, so a Nystrom Part
    # has to be built through the default inner constructor and populated by hand.
    function nystrom_sphere(a, N)
        disc = NystromDiscretisation(N)
        disc.force_pts .= fibonacci_ellipsoid(a, a, a, N)
        disc.velocity  .= Ref(zero(SVector{3,Float64}))
        Part(EllipsoidBody(a, a, a), disc,
             Frame(zero(SVector{3,Float64}), SMatrix{3,3,Float64,9}(I)))
    end

    # Tangential squirmer: u_s = B1 * sin(θ) ê_θ, written without the 1/sin(θ) so it stays
    # regular at the poles. Swims at U = (2/3) B1 ẑ.
    ez = SVector(0.0, 0.0, 1.0)
    squirm(p, a, B1) = B1 * (dot(p/a, ez) * (p/a) - ez)

    @testset "non-hybrid constructor" begin
        # Regression guard: this constructor used to reference an undefined `n` and a
        # typo'd `iVector`, so *any* non-hybrid Nystrom SwimmingProblem threw UndefVarError.
        part = nystrom_sphere(1.0, 60)
        prob = SwimmingProblem(MicroSwimmer([part]); eps=0.3)

        @test prob.disc isa NystromDiscretisation
        @test nf(prob.disc) == 60
        @test nq(prob.disc) == 60
        @test prob.pl_box === nothing         # nothing to preallocate when not hybrid
        @test prob.gmres_cache === nothing
    end

    @testset "solves, and stays force- and torque-free" begin
        a, B1 = 1.0, 1.0
        part = nystrom_sphere(a, 120)
        part.disc.velocity .= [squirm(p, a, B1) for p in part.disc.force_pts]
        prob = SwimmingProblem(MicroSwimmer([part]); eps=0.25)
        solve_problem!(prob)

        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8

        U = get_U(prob)
        @test abs(U[1]) < 1e-3 && abs(U[2]) < 1e-3    # swims along the symmetry axis
        @test U[3] > 0
    end

    @testset "converges to the analytic squirmer speed U = 2/3 B1" begin
        a, B1 = 1.0, 1.0
        analytic = 2/3 * B1
        errs = Float64[]
        for (N, eps) in [(120, 0.25), (250, 0.18), (500, 0.12)]
            part = nystrom_sphere(a, N)
            part.disc.velocity .= [squirm(p, a, B1) for p in part.disc.force_pts]
            prob = SwimmingProblem(MicroSwimmer([part]); eps=eps)
            solve_problem!(prob)
            push!(errs, abs(get_U(prob)[3] - analytic) / analytic)
        end
        @test issorted(errs, rev=true)     # error decreases monotonically with resolution
        @test errs[end] < 0.2              # collocation is lower order than the NEAREST scheme
    end

    @testset "hybrid path agrees with the dense solve" begin
        a, B1 = 1.0, 1.0
        function solve_sphere(hybrid)
            part = nystrom_sphere(a, 120)
            part.disc.velocity .= [squirm(p, a, B1) for p in part.disc.force_pts]
            prob = SwimmingProblem(MicroSwimmer([part]); eps=0.25, hybrid=hybrid)
            solve_problem!(prob)
            get_U(prob), get_Ω(prob)
        end
        Ud, Ωd = solve_sphere(false)
        Uh, Ωh = solve_sphere(true)
        @test Uh ≈ Ud atol=1e-10
        @test Ωh ≈ Ωd atol=1e-10
    end

    @testset "mul_swimming! matches assemble_swimming! on a Nystrom disc" begin
        a, B1 = 1.0, 1.0
        part = nystrom_sphere(a, 80)
        part.disc.velocity .= [squirm(p, a, B1) for p in part.disc.force_pts]
        prob = SwimmingProblem(MicroSwimmer([part]); eps=0.3)
        solve_problem!(prob)          # gathers the part into prob.disc

        disc = prob.disc
        x0   = prob.microswimmer.frame.location
        N3   = 3*nf(disc)
        A    = zeros(N3 + 6, N3 + 6)
        assemble_swimming!(A, x0, disc, prob.kernel; μ=prob.mu)

        Random.seed!(31)
        y = zeros(N3 + 6)
        for _ in 1:3
            x = randn(N3 + 6)
            mul_swimming!(y, x, x0, disc, prob.kernel; μ=prob.mu)
            @test y ≈ A * x atol=1e-10
        end
    end

    @testset "post-processing works on a Nystrom disc" begin
        # total_torque and shear_tensor used to @unpack quad_pts and nearest, neither of
        # which NystromDiscretisation has, so every one of these threw.
        a, B1 = 1.0, 1.0
        part = nystrom_sphere(a, 120)
        part.disc.velocity .= [squirm(p, a, B1) for p in part.disc.force_pts]
        prob = SwimmingProblem(MicroSwimmer([part]); eps=0.25)
        solve_problem!(prob)

        forces = get_forces(prob)
        @test length(forces) == 120
        @test total_force(forces, prob.disc) ≈ sum(forces) atol=1e-12

        S = stresslet_tensor(prob)
        @test S ≈ S' atol=1e-10
        @test abs(tr(S)) < 1e-10

        @test MicroSwimmers.total_power(prob) > 0          # squirming dissipates energy
        @test size(shear_tensor(forces, prob.disc)) == (3, 3)

        # a NystromDiscretisation is never weighted
        @test !MicroSwimmers.is_weighted(prob.disc)
        @test MicroSwimmers._fluid_quad_wts(prob.disc) === nothing
    end

end
