@testset "per-part regularisation parameter" begin

    # eps lives on the discretisation, one value per region, where a region is a contiguous
    # block of quadrature points: `quad_eps[j]` regularises `quad_part_ranges[j]`. The blob is
    # centred on the source (quadrature) point, so an interaction between two parts is
    # regularised at the *source's* eps, not the field point's.

    flagellum_model() = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)

    @testset "3-arg kernel agrees with the 2-arg form" begin
        xi = SVector(0.3, -0.2, 0.7)
        Xj = SVector(-0.1, 0.4, 0.2)
        for eps in (0.01, 0.1, 0.5)
            @test RegStokeslet(eps)(xi, Xj) ≈ RegStokeslet(eps)(xi, Xj, eps)
            @test MicroSwimmers.RegBlakelet(eps)(xi, Xj) ≈ MicroSwimmers.RegBlakelet(eps)(xi, Xj, eps)
            # and the third argument, not the field, is what takes effect
            @test RegStokeslet(0.9)(xi, Xj, eps) ≈ RegStokeslet(eps)(xi, Xj)
        end
    end

    @testset "set_eps! / get_eps" begin
        p = Part(EllipsoidBody(1.0, 1.0, 1.0), 20, 90; eps=0.07)
        @test get_eps(p) == 0.07
        @test p.disc.quad_eps == [0.07]
        set_eps!(p, 0.03)
        @test get_eps(p) == 0.03

        # a vaned flagellum is the one model with more than one region inside a single Part
        vaned = Part(PlanarVanedFlagellum(flagellum_model(), 0.2, 0.8, 0.15), 31, 151)
        @test length(vaned.disc.quad_part_ranges) == 2
        set_eps!(vaned, (0.05, 0.02))
        @test vaned.disc.quad_eps == [0.05, 0.02]
        @test get_eps(vaned) == [0.05, 0.02]        # vector, since the regions disagree
        @test_throws ArgumentError set_eps!(vaned, (0.05, 0.02, 0.01))
    end

    @testset "a region partition out of step with quad_eps is rejected" begin
        # zip would truncate to the shorter of the two and silently drop whole blocks of
        # quadrature points from the assembly, which is very hard to spot in the answer.
        body = Part(EllipsoidBody(1.0, 1.0, 1.0), 20, 90; eps=0.1)
        prob = ResistanceProblem(MicroSwimmer([body]))
        gather!(prob)
        N = nf(prob.disc)
        push!(prob.disc.quad_part_ranges, 1:0)      # a region with no eps of its own
        @test_throws ArgumentError assemble!(zeros(3N, 3N), prob.disc, prob.kernel; μ=prob.mu)
    end

    @testset "uniform per-region eps reproduces the scalar path exactly" begin
        # Guards the `ranges === nothing` fallback in eps_regions: partitioning the quadrature
        # must not change a single bit when every region carries the same value.
        body = Part(EllipsoidBody(1.0, 1.0, 1.0), 40, 170; eps=0.08)
        flag = Part(flagellum_model(), 15, 71; eps=0.08, location=[1.0, 0.0, 0.0])
        prob = ResistanceProblem(MicroSwimmer([body, flag]))
        gather!(prob)

        N, Q = nf(prob.disc), nq(prob.disc)
        A_regions = zeros(3N, 3N)
        A_scalar  = zeros(3N, 3N)
        assemble!(A_regions, prob.disc, prob.kernel; μ=prob.mu)
        assemble!(A_scalar, prob.disc.force_pts, prob.disc.quad_pts, prob.disc.nearest,
                  prob.disc.quad_wts, RegStokeslet(0.08); μ=prob.mu)
        @test A_regions == A_scalar
    end

    @testset "heterogeneous eps uses the source region's value" begin
        body = Part(EllipsoidBody(1.0, 1.0, 1.0), 25, 110; eps=0.2)
        flag = Part(flagellum_model(), 12, 61; eps=0.03, location=[1.0, 0.0, 0.0])
        prob = ResistanceProblem(MicroSwimmer([body, flag]))
        gather!(prob)

        N = nf(prob.disc)
        A = zeros(3N, 3N)
        assemble!(A, prob.disc, prob.kernel; μ=prob.mu)

        # hand-build the same matrix, picking eps by which region the *quadrature* point is in
        expected = zeros(3N, 3N)
        for (j, rng) in enumerate(prob.disc.quad_part_ranges)
            ε = prob.disc.quad_eps[j]
            for q in rng, m in 1:N
                S   = RegStokeslet(ε)(prob.disc.force_pts[m], prob.disc.quad_pts[q])
                col = 3*prob.disc.nearest[q] - 2
                expected[3m-2:3m, col:col+2] .+= S
            end
        end
        expected ./= 8π * prob.mu
        @test A ≈ expected atol=1e-14

        # and it genuinely differs from using either eps everywhere
        for ε in (0.2, 0.03)
            A_uniform = zeros(3N, 3N)
            assemble!(A_uniform, prob.disc.force_pts, prob.disc.quad_pts, prob.disc.nearest,
                      prob.disc.quad_wts, RegStokeslet(ε); μ=prob.mu)
            @test !isapprox(A, A_uniform; rtol=1e-6)
        end
    end

    @testset "mul_swimming! matches assemble_swimming! under heterogeneous eps" begin
        # The hybrid solver preconditions GMRES with a dense LU, so the matrix-free matvec has
        # to see exactly the same per-region eps the dense assembly does.
        body = Part(EllipsoidBody(1.0, 1.0, 1.0), 20, 90; eps=0.2)
        flag = Part(flagellum_model(), 10, 51; eps=0.03, location=[1.0, 0.0, 0.0])
        prob = SwimmingProblem(MicroSwimmer([body, flag]))
        gather!(prob)

        N  = nf(prob.disc)
        x0 = zero(SVector{3,Float64})
        A  = zeros(3N+6, 3N+6)
        assemble_swimming!(A, x0, prob.disc, prob.kernel; μ=prob.mu)

        rng = Random.MersenneTwister(42)
        x   = randn(rng, 3N+6)
        y   = zeros(3N+6)
        mul_swimming!(y, x, x0, prob.disc, prob.kernel; μ=prob.mu)
        @test y ≈ A * x rtol=1e-10
    end

    @testset "a vaned flagellum carries two regions into the problem" begin
        vaned = Part(PlanarVanedFlagellum(flagellum_model(), 0.2, 0.8, 0.15), 31, 151)
        set_eps!(vaned, (0.05, 0.02))
        prob = SwimmingProblem(MicroSwimmer([vaned]))

        # one Part, but two regions at the problem level — this is what the per-Part range
        # partition could not represent
        @test length(prob.disc.quad_part_ranges) == 2
        @test prob.disc.quad_eps == [0.05, 0.02]
        @test last(last(prob.disc.quad_part_ranges)) == nq(prob.disc)

        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "problem-level eps warns and is ignored" begin
        body = Part(EllipsoidBody(1.0, 1.0, 1.0), 20, 90; eps=0.06)
        ms   = MicroSwimmer([body])

        # the warning carries maxlog=1, so ask for it on a fresh value each time
        @test_logs (:warn,) match_mode=:any ResistanceProblem(ms; eps=0.99)

        quiet   = ResistanceProblem(ms)
        ignored = ResistanceProblem(ms; eps=0.99)
        @test quiet.disc.quad_eps == ignored.disc.quad_eps == [0.06]
        @test quiet.kernel.eps == ignored.kernel.eps == 0.06

        add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
        solve_problem!(quiet); solve_problem!(ignored)
        @test quiet.force_vals == ignored.force_vals
    end

    @testset "Part(model) sizes itself from hf, eps and hq" begin
        # Force points use the open rule hf = L/(N+1), quadrature the closed one hq = L/(Q-1).
        L    = 10.0
        flag = Part(PlanarFlagellum(L, 0.0, 0.6, 0.5, π/2, 2π, 2π, 0.0))
        @test L / (nf(flag.disc) + 1) ≈ DEFAULT_HF     rtol=0.1
        @test L / (nq(flag.disc) - 1) ≈ 2*DEFAULT_EPS  rtol=0.1
        @test get_eps(flag) == DEFAULT_EPS

        # surfaces spread N points over the area, so h = sqrt(area/N)
        m    = EllipsoidBody(5.0, 4.0, 4.0)
        body = Part(m)
        @test sqrt(surface_area(m) / nf(body.disc)) ≈ DEFAULT_HF     rtol=0.1
        @test sqrt(surface_area(m) / nq(body.disc)) ≈ 2*DEFAULT_EPS  rtol=0.1

        # raising eps past hf/4 makes the patch-resolution constraint bind instead, keeping
        # Q at roughly 4N rather than letting it collapse to N
        coarse = Part(m; eps=0.25)
        @test nq(coarse.disc) / nf(coarse.disc) ≈ 4 rtol=0.1

        # explicit hf/hq override the defaults
        fine = Part(m; hf=1.0, hq=0.5)
        @test sqrt(surface_area(m) / nf(fine.disc)) ≈ 1.0 rtol=0.1
        @test sqrt(surface_area(m) / nq(fine.disc)) ≈ 0.5 rtol=0.1

        # a model with no arclength or surface_area says so rather than guessing
        @test_throws ArgumentError npoints_for_spacing(
            CylindricalGroovedBody(1.0, 1.0, 1.0, 0.5, 0.5, 0.2), 0.5, false)
    end

    @testset "ellipsoid_area" begin
        @test ellipsoid_area(1.0, 1.0, 1.0) ≈ 4π rtol=1e-12
        @test ellipsoid_area(2.0, 2.0, 2.0) ≈ 16π rtol=1e-12
        # Thomsen's approximation is within ~1.1% of the exact prolate area
        a, c = 1.0, 2.0
        e    = sqrt(1 - a^2/c^2)
        exact = 2π*a^2 * (1 + (c/(a*e)) * asin(e))
        @test ellipsoid_area(a, a, c) ≈ exact rtol=0.011
    end

end
