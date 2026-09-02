@testset "problems.jl — Stokes' law analytic limit" begin

    # A sphere of radius `a` translating/rotating in quiescent Stokes flow has
    # known drag: F = -6*pi*mu*a*U (translation), T = -8*pi*mu*a^3*Ω (rotation).
    # The regularised-stokeslet BEM solution converges to this as N,Q grow and
    # eps shrinks; we check the numerical drag approaches the analytic value
    # as resolution increases, rather than pinning a single fixed tolerance.

    function sphere_problem(a, N, Q; eps)
        body = Part(EllipsoidBody(a, a, a), N, Q; eps=eps)
        ms = MicroSwimmer([body])
        prob = ResistanceProblem(ms)
        prob, body
    end

    @testset "translational drag" begin
        a  = 1.0
        μ  = 1.0
        Uvec = SVector(1.0, 0.0, 0.0)
        analytic = 6π * μ * a

        errs = Float64[]
        for (N, Q, eps) in [(20, 90, 0.08), (40, 170, 0.05), (80, 330, 0.03)]
            prob, body = sphere_problem(a, N, Q; eps=eps)
            add_rigid_body_motion!(body, Uvec, zero(SVector{3,Float64}))
            solve_problem!(prob)
            F, _ = total_force_and_torque(prob)
            push!(errs, abs(norm(F) - analytic) / analytic)
        end

        @test errs[end] < errs[1]              # error decreases with resolution
        @test errs[end] < 0.1                  # within 10% at the finest resolution tested
    end

    @testset "rotational drag" begin
        a  = 1.0
        μ  = 1.0
        Ωvec = SVector(0.0, 0.0, 1.0)
        analytic = 8π * μ * a^3

        errs = Float64[]
        for (N, Q, eps) in [(20, 90, 0.08), (40, 170, 0.05), (80, 330, 0.03)]
            prob, body = sphere_problem(a, N, Q; eps=eps)
            add_rigid_body_motion!(body, zero(SVector{3,Float64}), Ωvec)
            solve_problem!(prob)
            _, T = total_force_and_torque(prob)
            push!(errs, abs(norm(T) - analytic) / analytic)
        end

        @test errs[end] < errs[1]
        @test errs[end] < 0.1
    end

    @testset "grand_resistance_matrix on a sphere" begin
        a = 1.0
        body = Part(EllipsoidBody(a, a, a), 60, 250; eps=0.04)
        ms = MicroSwimmer([body])
        R = grand_resistance_matrix(ms)

        trans_drag = 6π * a
        rot_drag   = 8π * a^3

        # translation-translation block ~ 6*pi*mu*a * I, rotation-rotation ~ 8*pi*mu*a^3 * I
        @test R[1:3, 1:3] ≈ trans_drag * I atol=0.1*trans_drag
        @test R[4:6, 4:6] ≈ rot_drag * I atol=0.1*rot_drag

        # translation-rotation coupling vanishes for a sphere
        @test norm(R[1:3, 4:6]) < 0.05*trans_drag
        @test norm(R[4:6, 1:3]) < 0.05*rot_drag
    end

end
