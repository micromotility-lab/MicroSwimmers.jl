@testset "maths.jl" begin

    @testset "rotation_matrix" begin
        Random.seed!(1)
        for _ in 1:20
            axis  = normalize(randn(3))
            angle = 2π * rand() - π
            R = rotation_matrix(axis, angle)

            @test R' * R ≈ I atol=1e-12
            @test det(R) ≈ 1.0 atol=1e-12
            @test R * axis ≈ axis atol=1e-10          # fixed point on the rotation axis
        end

        @test rotation_matrix([0.0, 0.0, 1.0], 0.0) ≈ I atol=1e-14

        # composing two rotations about the same axis adds the angles
        axis = [0.0, 0.0, 1.0]
        R1 = rotation_matrix(axis, 0.3)
        R2 = rotation_matrix(axis, 0.7)
        @test R1 * R2 ≈ rotation_matrix(axis, 1.0) atol=1e-12

        # 90 degree rotation about z maps x -> y
        R90 = rotation_matrix([0.0, 0.0, 1.0], π/2)
        @test R90 * [1.0, 0.0, 0.0] ≈ [0.0, 1.0, 0.0] atol=1e-12
    end

    @testset "skew_symmetric_static" begin
        Random.seed!(2)
        for _ in 1:20
            x = SVector{3}(randn(3))
            y = SVector{3}(randn(3))
            K = MicroSwimmers.skew_symmetric_static(x)
            @test K * y ≈ cross(x, y) atol=1e-12
            @test K ≈ -K' atol=1e-12
        end
    end

    @testset "smooth_max" begin
        for _ in 1:10
            a, b = randn(), randn()
            for k in (1.0, 10.0, 100.0)
                @test MicroSwimmers.smooth_max(a, b, k) >= max(a, b) - 1e-12
            end
            # converges to max(a,b) as k grows
            @test abs(MicroSwimmers.smooth_max(a, b, 1e6) - max(a, b)) < 1e-4
        end
    end

    @testset "cumulative_trapz!" begin
        x = collect(0.0:0.1:1.0)

        # constant integrand: exact
        y = fill(2.0, length(x))
        out = similar(x)
        MicroSwimmers.cumulative_trapz!(out, x, y)
        @test out ≈ 2.0 .* x atol=1e-12
        @test out[1] == 0.0

        # linear integrand y=x: exact for trapezoidal rule, ∫_0^s t dt = s^2/2
        y = x
        MicroSwimmers.cumulative_trapz!(out, x, y)
        @test out ≈ x .^ 2 ./ 2 atol=1e-12
    end

    @testset "Helix" begin
        # a helix with zero radius and zero forward velocity is a fixed point
        p = (0.0, 0.0, 0.0, 0.0, 1.0, π/2, 0.0, 0.0, 0.0)
        ts = collect(0.0:0.25:2.0)
        pts = MicroSwimmers.helix(ts, p)
        @test all(pt -> pt ≈ [0.0, 0.0, 0.0], pts)

        # pure translation along z-axis (θ=0) with zero radius
        p2 = (0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
        pts2 = MicroSwimmers.helix(ts, p2)
        @test all(i -> pts2[i] ≈ [0.0, 0.0, ts[i]], eachindex(ts))

        h = MicroSwimmers.Helix(0.0, 0.0, 0.0, 1.0, 2π, π/2, 0.0, 0.5, 0.0)
        @test radius(h) == 0.5
        @test axis_velocity(h) == 1.0
        @test axis_angular_velocity(h) == 2π
        @test pitch(h) ≈ 1.0 atol=1e-12          # P = 2π v/ω = 2π*1/(2π) = 1
        @test chirality_sign(h) == 1
        @test h(ts) ≈ MicroSwimmers.helix(ts, Tuple(h))
    end

end
