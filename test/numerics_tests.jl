@testset "numerics.jl" begin

    @testset "RegStokeslet" begin
        eps = 0.1
        k = RegStokeslet(eps)

        Random.seed!(4)
        for _ in 1:10
            xi = SVector{3}(randn(3))
            Xj = SVector{3}(randn(3))
            S = k(xi, Xj)
            @test S ≈ S' atol=1e-12                       # symmetric kernel

            # self-interaction closed form: R=0 => diag = 2eps^2/eps^3 = 2/eps
            Sself = k(xi, xi)
            @test Sself ≈ (2/eps) * I atol=1e-10
        end

        # converges to the singular (unregularised) Oseen tensor as eps -> 0, R != 0
        xi = SVector(1.0, 0.0, 0.0)
        Xj = SVector(0.0, 0.0, 0.0)
        R  = xi - Xj
        Sbig = MMatrix{3,3,Float64}(undef)
        MicroSwimmers.stokeslet!(Sbig, R)
        for eps_small in (1e-2, 1e-4, 1e-6)
            Sreg = RegStokeslet(eps_small)(xi, Xj)
            @test Sreg ≈ Matrix(Sbig) atol=10*eps_small^2
        end

        # regularised_stokeslet! (in-place legacy form) agrees with the callable kernel
        Sip = MMatrix{3,3,Float64}(undef)
        regularised_stokeslet!(Sip, R; eps=eps)
        @test Matrix(Sip) ≈ Matrix(k(xi, Xj)) atol=1e-12
    end

    @testset "RegBlakelet wall limit" begin
        eps = 0.05
        xi = SVector(0.3, -0.2, 1.0)

        # as the source height above the wall (z=0) grows, the wall correction
        # should vanish and RegBlakelet should approach plain RegStokeslet
        errs = Float64[]
        for h in (10.0, 100.0, 1000.0)
            Xj = SVector(0.1, 0.4, h)
            Sblake = MicroSwimmers.RegBlakelet(eps)(xi, Xj)
            Sfree  = RegStokeslet(eps)(xi, Xj)
            push!(errs, norm(Sblake - Sfree))
        end
        @test issorted(errs, rev=true)     # error strictly decreasing as h grows
        @test errs[end] < 1e-2
    end

    @testset "assemble!" begin
        eps = 0.2
        μ   = 1.3
        kernel = RegStokeslet(eps)

        force_pts = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)]
        quad_pts  = [SVector(0.1, 0.0, 0.0), SVector(0.2, 0.1, 0.0), SVector(0.9, 0.0, 0.1)]
        nearest   = [1, 1, 2]     # first two quad points nearest force pt 1, third nearest force pt 2

        A = zeros(6, 6)
        assemble!(A, force_pts, quad_pts, nearest, kernel; μ=μ)

        expected = zeros(6, 6)
        expected[1:3, 1:3] .= Matrix(kernel(force_pts[1], quad_pts[1])) .+ Matrix(kernel(force_pts[1], quad_pts[2]))
        expected[1:3, 4:6] .= Matrix(kernel(force_pts[1], quad_pts[3]))
        expected[4:6, 1:3] .= Matrix(kernel(force_pts[2], quad_pts[1])) .+ Matrix(kernel(force_pts[2], quad_pts[2]))
        expected[4:6, 4:6] .= Matrix(kernel(force_pts[2], quad_pts[3]))
        expected ./= (8π*μ)

        @test A ≈ expected atol=1e-12
    end

    @testset "assemble_swimming! rigid-body structure" begin
        eps = 0.2
        μ   = 1.0
        kernel = RegStokeslet(eps)

        force_pts = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)]
        quad_pts  = [SVector(0.05, 0.0, 0.0), SVector(0.9, 0.0, 0.0), SVector(0.0, 0.9, 0.1)]
        nearest   = [1, 2, 3]
        x0 = SVector(0.2, 0.2, 0.0)
        N  = length(force_pts)
        N3 = 3N

        A = zeros(N3 + 6, N3 + 6)
        assemble_swimming!(A, x0, force_pts, quad_pts, nearest, kernel; μ=μ)

        # top-left block matches assemble!
        Atop = zeros(N3, N3)
        assemble!(Atop, force_pts, quad_pts, nearest, kernel; μ=μ)
        @test A[1:N3, 1:N3] ≈ Atop atol=1e-12

        # rigid-body columns: U columns = -I, Omega columns = skew(xm - x0)
        for (m, xm) in enumerate(force_pts)
            row = 3m - 2
            @test A[row:row+2, N3+1:N3+3] ≈ -Matrix(I, 3, 3) atol=1e-12
            @test A[row:row+2, N3+4:N3+6] ≈ Matrix(MicroSwimmers.skew_symmetric_static(xm - x0)) atol=1e-12
        end

        # force-free / torque-free rows
        Ftest = zeros(3, N3)
        Ttest = zeros(3, N3)
        for (q, yq) in enumerate(quad_pts)
            n = nearest[q]
            Ftest[:, 3n-2:3n] .+= Matrix(I, 3, 3)
            Ttest[:, 3n-2:3n] .+= Matrix(MicroSwimmers.skew_symmetric_static(yq - x0))
        end
        @test A[N3+1:N3+3, 1:N3] ≈ Ftest atol=1e-12
        @test A[N3+4:N3+6, 1:N3] ≈ Ttest atol=1e-12
    end

end
