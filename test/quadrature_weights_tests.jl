@testset "quadrature weights" begin

    # Exact prolate spheroid results (equatorial semi-axis b, polar semi-axis a > b).
    # Used as an independent reference for both the surface area the ray-march weights
    # should integrate to, and the drag the weighted BEM should reproduce.
    function prolate_area(b, a)
        e = sqrt(1 - b^2/a^2)
        2π*b^2*(1 + (a/(b*e))*asin(e))
    end
    function prolate_drag_axial(b, a; mu=1.0)
        e = sqrt(1 - b^2/a^2)
        16π*mu*a*e^3 / (-2e + (1+e^2)*log((1+e)/(1-e)))
    end
    function prolate_drag_transverse(b, a; mu=1.0)
        e = sqrt(1 - b^2/a^2)
        32π*mu*a*e^3 / (2e + (3e^2-1)*log((1+e)/(1-e)))
    end

    raymarch(m, N) = begin
        pts, wts = SVector{3,Float64}[], Float64[]
        MicroSwimmers.raymarch_cloud!(pts, wts, x -> MicroSwimmers.implicit(m, x),
                                      MicroSwimmers.bounding_radius(m), N;
                                      seed=MicroSwimmers.seed(m))
        pts, wts
    end

    @testset "raymarch weights integrate to the surface area" begin
        # A sphere ray-marched from its centre hits at t = a with n parallel to d, so every
        # weight is exactly a^2 * 4pi/N — the area identity holds to machine precision.
        for a in (1.0, 2.5)
            _, wts = raymarch(ImplicitEllipsoid(a, a, a), 2000)
            @test sum(wts) ≈ 4π*a^2 rtol=1e-10
            @test all(w -> w ≈ first(wts), wts)      # uniform for a sphere
        end

        # A prolate spheroid is the real test: t and |n·d| both vary along the surface.
        b, a = 1.0, 2.0
        _, wts = raymarch(ImplicitEllipsoid(b, b, a), 4000)
        @test sum(wts) ≈ prolate_area(b, a) rtol=1e-5
        # guard against a silently-constant weight sneaking in
        @test maximum(wts) / minimum(wts) > 3
    end

    @testset "total_quad_weight" begin
        a = 1.0
        bw = Part(ImplicitEllipsoid(a, a, a), 200, 900; weighted=true)
        bu = Part(ImplicitEllipsoid(a, a, a), 200, 900)

        @test bw.disc.quad_wts !== nothing
        @test bu.disc.quad_wts === nothing
        @test total_quad_weight(bw.disc) ≈ 4π*a^2 rtol=1e-10
        # an unweighted discretisation has unit weights, so this counts quadrature points
        @test total_quad_weight(bu.disc) == nq(bu.disc)
    end

    @testset "weighted=true is rejected for models with no area elements" begin
        @test_throws ArgumentError Part(EllipsoidBody(1.0, 1.0, 1.0), 20, 90; weighted=true)
        @test supports_quadrature_weights(ImplicitEllipsoid(1.0, 1.0, 1.0))
        @test !supports_quadrature_weights(EllipsoidBody(1.0, 1.0, 1.0))
    end

    # ---------------------------------------------------------------- assembly

    force_pts = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)]
    quad_pts  = [SVector(0.1, 0.0, 0.0), SVector(0.2, 0.1, 0.0), SVector(0.9, 0.0, 0.1)]
    nearest   = [1, 1, 2]
    wts       = [0.3, 1.7, 0.8]
    kernel    = RegStokeslet(0.2)
    μ         = 1.3

    @testset "absent weights are the absorbed path" begin
        A5 = zeros(6, 6); assemble!(A5, force_pts, quad_pts, nearest, kernel; μ=μ)
        A6 = zeros(6, 6); assemble!(A6, force_pts, quad_pts, nearest, nothing, kernel; μ=μ)
        @test A5 == A6                                  # bit-identical, not just approx

        Aones = zeros(6, 6)
        assemble!(Aones, force_pts, quad_pts, nearest, ones(3), kernel; μ=μ)
        @test Aones ≈ A5 atol=1e-14
    end

    @testset "weighted assemble!" begin
        A = zeros(6, 6)
        assemble!(A, force_pts, quad_pts, nearest, wts, kernel; μ=μ)

        expected = zeros(6, 6)
        expected[1:3, 1:3] .= wts[1] .* Matrix(kernel(force_pts[1], quad_pts[1])) .+
                              wts[2] .* Matrix(kernel(force_pts[1], quad_pts[2]))
        expected[1:3, 4:6] .= wts[3] .* Matrix(kernel(force_pts[1], quad_pts[3]))
        expected[4:6, 1:3] .= wts[1] .* Matrix(kernel(force_pts[2], quad_pts[1])) .+
                              wts[2] .* Matrix(kernel(force_pts[2], quad_pts[2]))
        expected[4:6, 4:6] .= wts[3] .* Matrix(kernel(force_pts[2], quad_pts[3]))
        expected ./= (8π*μ)

        @test A ≈ expected atol=1e-12
    end

    @testset "weighted assemble_swimming! constraint rows" begin
        fp = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)]
        qp = [SVector(0.05, 0.0, 0.0), SVector(0.9, 0.0, 0.0), SVector(0.0, 0.9, 0.1), SVector(0.1, 0.1, 0.0)]
        nn = [1, 2, 3, 1]
        ww = [0.4, 1.2, 0.7, 2.1]
        x0 = SVector(0.2, 0.2, 0.0)
        N3 = 3*length(fp)

        A = zeros(N3 + 6, N3 + 6)
        assemble_swimming!(A, x0, fp, qp, nn, ww, kernel; μ=μ)

        # BEM block matches the weighted assemble!
        Atop = zeros(N3, N3)
        assemble!(Atop, fp, qp, nn, ww, kernel; μ=μ)
        @test A[1:N3, 1:N3] ≈ Atop atol=1e-12

        # rigid-body columns are unaffected by quadrature weights
        for (m, xm) in enumerate(fp)
            row = 3m - 2
            @test A[row:row+2, N3+1:N3+3] ≈ -Matrix(I, 3, 3) atol=1e-12
            @test A[row:row+2, N3+4:N3+6] ≈ Matrix(MicroSwimmers.skew_symmetric_static(xm - x0)) atol=1e-12
        end

        # force-free / torque-free rows must carry the same weights as the BEM block
        Ftest = zeros(3, N3); Ttest = zeros(3, N3)
        for (q, yq) in enumerate(qp)
            n = nn[q]
            Ftest[:, 3n-2:3n] .+= ww[q] .* Matrix(I, 3, 3)
            Ttest[:, 3n-2:3n] .+= ww[q] .* Matrix(MicroSwimmers.skew_symmetric_static(yq - x0))
        end
        @test A[N3+1:N3+3, 1:N3] ≈ Ftest atol=1e-12
        @test A[N3+4:N3+6, 1:N3] ≈ Ttest atol=1e-12
    end

    @testset "weighted mul_swimming! matches weighted assemble_swimming!" begin
        fp = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)]
        qp = [SVector(0.05, 0.0, 0.0), SVector(0.9, 0.0, 0.0), SVector(0.0, 0.9, 0.1), SVector(0.4, 0.4, 0.0)]
        nn = [1, 2, 3, 1]
        ww = [0.4, 1.2, 0.7, 2.1]
        x0 = SVector(0.2, 0.2, 0.0)
        N3 = 3*length(fp)

        A = zeros(N3 + 6, N3 + 6)
        assemble_swimming!(A, x0, fp, qp, nn, ww, kernel; μ=μ)

        Random.seed!(23)
        y = zeros(N3 + 6)
        for _ in 1:5
            x = randn(N3 + 6)
            mul_swimming!(y, x, x0, fp, qp, nn, ww, kernel; μ=μ)
            @test y ≈ A * x atol=1e-10
        end
    end

    # ---------------------------------------------------------------- end to end

    @testset "Stokes' law for a weighted implicit sphere" begin
        a = 1.0
        errs = Float64[]
        for (N, Q, eps) in [(100, 450, 0.08), (200, 900, 0.05)]
            body = Part(ImplicitEllipsoid(a, a, a), N, Q; weighted=true, eps=eps)
            prob = ResistanceProblem(MicroSwimmer([body]))
            add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
            solve_problem!(prob)
            F, _ = total_force_and_torque(prob)
            push!(errs, abs(norm(F) - 6π*a) / (6π*a))
        end
        @test errs[end] < errs[1]
        @test errs[end] < 0.05
    end

    @testset "prolate spheroid drag with genuinely non-uniform weights" begin
        # Unlike a sphere, this body's weights vary by a factor of ~4 across the surface,
        # so a wrong weight convention cannot cancel out of the answer.
        b, a = 1.0, 2.0
        body = Part(ImplicitEllipsoid(b, b, a), 300, 1400; weighted=true, eps=0.04)
        @test maximum(body.disc.quad_wts) / minimum(body.disc.quad_wts) > 3

        ms = MicroSwimmer([body])
        prob = ResistanceProblem(ms)

        add_rigid_body_motion!(body, SVector(0.0, 0.0, 1.0), zero(SVector{3,Float64}))
        solve_problem!(prob)
        F, _ = total_force_and_torque(prob)
        @test norm(F) ≈ prolate_drag_axial(b, a) rtol=0.05

        add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
        solve_problem!(prob)
        F, _ = total_force_and_torque(prob)
        @test norm(F) ≈ prolate_drag_transverse(b, a) rtol=0.05
    end

    @testset "mixed weighted/unweighted swimmer stays force- and torque-free" begin
        # The sharpest test in the file: it can only pass if the weighted constraint rows in
        # assemble_swimming! and the weighted total_force/total_torque agree on the
        # convention, across a swimmer whose two parts use *different* conventions.
        body = Part(ImplicitEllipsoid(0.5, 0.5, 0.5), 120, 500; weighted=true, eps=0.1)
        flag = Part(PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0), 20, 101;
                    eps=0.1, location=SVector(0.0, 0.0, -0.5))
        prob = SwimmingProblem(MicroSwimmer([body, flag]))

        for t in (0.0, 0.37)
            update_boundary!(prob, t)
            solve_problem!(prob)
            F, T = total_force_and_torque(prob)
            @test norm(F) < 1e-8
            @test norm(T) < 1e-8
        end

        # the global weight vector interleaves the two parts correctly: the body range
        # integrates to its surface area, the flagellum range keeps unit weights
        @test total_quad_weight(prob.disc, 1) ≈ 4π*0.5^2 rtol=1e-10
        @test total_quad_weight(prob.disc, 2) == nq(prob.disc, 2)
    end

    @testset "unit weights reproduce the unweighted solve exactly" begin
        a = 1.0
        function drag(weighted)
            body = Part(ImplicitEllipsoid(a, a, a), 150, 700; weighted=weighted, eps=0.06)
            weighted && fill!(body.disc.quad_wts, 1.0)
            prob = ResistanceProblem(MicroSwimmer([body]))
            add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
            solve_problem!(prob)
            first(total_force_and_torque(prob))
        end
        @test drag(true) ≈ drag(false) atol=1e-12
    end

    @testset "total_power identity holds on the weighted path" begin
        body = Part(ImplicitEllipsoid(1.0, 1.0, 1.0), 200, 900; weighted=true, eps=0.05)
        prob = ResistanceProblem(MicroSwimmer([body]))
        U = SVector(0.3, -0.2, 0.5)
        add_rigid_body_motion!(body, U, zero(SVector{3,Float64}))
        solve_problem!(prob)

        F, _ = total_force_and_torque(prob)
        @test total_power(prob) ≈ dot(F, U) rtol=1e-10
    end

    @testset "FluidVelocity carries the weights through" begin
        # Evaluating the field at a force point reproduces that node's boundary velocity,
        # since that is exactly the row of the system that was solved. It only holds if the
        # field applies the same quadrature weights the solve did.
        body = Part(ImplicitEllipsoid(1.0, 1.0, 1.0), 120, 500; weighted=true, eps=0.08)
        prob = ResistanceProblem(MicroSwimmer([body]))
        add_rigid_body_motion!(body, SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
        solve_problem!(prob)

        fv = FluidVelocity(prob)
        @test fv.quad_wts !== nothing
        for m in (1, 37, 120)
            @test fv(prob.disc.force_pts[m]) ≈ prob.disc.velocity[m] atol=1e-9
        end
    end

end
