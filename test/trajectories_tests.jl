@testset "trajectories.jl" begin

    function linear_trajectory(; n=11, v=SVector(1.0, 2.0, -1.0))
        t = collect(range(0.0, 2.0, n))
        x = [v * ti for ti in t]
        b1 = [SVector(1.0, 0.0, 0.0) for _ in 1:n]
        b2 = [SVector(0.0, 1.0, 0.0) for _ in 1:n]
        Trajectory(t, x, b1, b2, false)
    end

    @testset "average_swimming_velocity" begin
        v = SVector(1.0, 2.0, -1.0)
        traj = linear_trajectory(; v=v)
        @test average_swimming_velocity(traj) ≈ v atol=1e-12
    end

    @testset "translate_trajectory / centred_trajectory" begin
        traj = linear_trajectory()
        shift = SVector(1.0, 1.0, 1.0)
        traj2 = translate_trajectory(traj, shift)
        @test traj2.x ≈ traj.x .- Ref(shift) atol=1e-12
        @test traj2.b1 == traj.b1 && traj2.b2 == traj.b2

        c = centred_trajectory(traj)
        @test norm(sum(c.x)) < 1e-10       # mean position at the origin
    end

    @testset "running_mean" begin
        traj = linear_trajectory()
        rm = running_mean(traj, 1)
        # for a perfectly linear trajectory a running mean leaves interior points unchanged
        @test rm.x[2:end-1] ≈ traj.x[2:end-1] atol=1e-10
    end

    @testset "helix(ts, p) matches Helix struct evaluation" begin
        p = (0.1, -0.2, 0.3, 0.5, 1.0, π/3, π/4, 0.7, 0.2)
        ts = collect(0.0:0.2:2.0)
        h = Helix(p...)
        @test h(ts) ≈ helix(ts, p) atol=1e-12
    end

    @testset "fit_helix recovers exact synthetic helix parameters" begin
        p_true = (0.0, 0.0, 0.0, 0.3, 2π, π/2.2, 0.3, 0.6, 0.1)
        ts = collect(range(0.0, 1.0, 60))
        pts = helix(ts, p_true)
        b1 = [SVector(1.0, 0.0, 0.0) for _ in ts]
        b2 = [SVector(0.0, 1.0, 0.0) for _ in ts]
        traj = Trajectory(ts, pts, b1, b2, false)

        h = fit_helix(traj; smooth=false)

        @test radius(h) ≈ p_true[8] atol=1e-3
        @test axis_velocity(h) ≈ p_true[4] atol=1e-3
        @test axis_angular_velocity(h) ≈ p_true[5] atol=1e-3
    end

    @testset "continue_periodic_trajectory" begin
        n = 9
        t = collect(range(0.0, 1.0, n))
        # a closed periodic loop in b1/b2 (rotate by pi/4 per period) with net translation
        x = [SVector(ti, 0.0, 0.0) for ti in t]
        b1 = [SVector(cos(π/2*ti), sin(π/2*ti), 0.0) for ti in t]
        b2 = [SVector(-sin(π/2*ti), cos(π/2*ti), 0.0) for ti in t]
        traj = Trajectory(t, x, b1, b2, true)

        ext = continue_periodic_trajectory(traj, 3)
        @test length(ext.t) == 3*(n-1) + 1
        @test ext.t[1:n] ≈ traj.t atol=1e-12
        @test ext.x[1:n] ≈ traj.x atol=1e-12

        # continue_periodic_trajectory! extends in place and agrees with the pure version
        traj2 = deepcopy(traj)
        continue_periodic_trajectory!(traj2, 3)
        @test traj2.t ≈ ext.t atol=1e-12
        @test traj2.x ≈ ext.x atol=1e-12
    end

end
