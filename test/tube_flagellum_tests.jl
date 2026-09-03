@testset "tube flagella" begin

    # A LineTubeFlagellum keeps the force points on the centreline and moves only the
    # quadrature points onto a ring of radius `a`. The linear system is therefore unchanged
    # in size and in its collocation points — the accuracy comes for free in solve cost.

    beating()  = PlanarFlagellum(1.0, 0.0, 0.3, 0.15, 2π, 2π, 2π, 0.0)
    straight(L) = PlanarFlagellum(L, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)   # tangent ≡ x̂
    curved(L)   = PlanarFlagellum(L, 0.0, 0.6, 0.3, 2π, 2π, 0.0, 0.0)     # bent, but static

    # translating drag coefficients, along the filament and across it
    function drag(part)
        R = grand_resistance_matrix(MicroSwimmer([part]))
        (abs(R[1,1]), abs(R[2,2]))
    end

    @testset "the force half is exactly the untubed model" begin
        # This is the defining property of a LineTube: only the quadrature moves. If this
        # ever stops holding bitwise, the tube has started changing the unknowns as well as
        # the integration, and it is no longer the cheap variant it claims to be.
        m, N, Q = beating(), 21, 101
        plain = Part(m, N, Q; eps=0.05)
        tube  = Part(LineTubeFlagellum(m, 0.02; Q_cs=12), N, Q; eps=0.002)

        @test nf(tube.disc) == nf(plain.disc) == N
        @test tube.disc.force_pts == plain.disc.force_pts
        @test tube.disc.velocity  == plain.disc.velocity
        @test nq(tube.disc) == Q * 12
    end

    @testset "march_stations agrees with integrate_centreline! bitwise" begin
        # march_stations duplicates integrate_centreline!'s arithmetic so a ring can have the
        # tangent as well as the position. This test is the only thing stopping the two from
        # drifting apart.
        m = beating()
        for (N, inc) in ((21, false), (101, true), (7, true))
            expected = Vector{SVector{3,Float64}}(undef, N)
            MicroSwimmers.integrate_centreline!(expected, m, 0.37; include_endpoints=inc)

            got = Vector{SVector{3,Float64}}(undef, N)
            MicroSwimmers.march_stations(m, N, 0.37; include_endpoints=inc) do i, s, x, τ
                got[i] = x
            end
            @test got == expected
        end
    end

    @testset "ring geometry" begin
        m, Q, Q_cs, a = beating(), 101, 12, 0.02
        tube = Part(LineTubeFlagellum(m, a; Q_cs=Q_cs), 21, Q; eps=0.002)

        # the quadrature stations are the untubed model's closed-rule nodes
        centres = Vector{SVector{3,Float64}}(undef, Q)
        m(centres, 0.0; include_endpoints=true)

        # every ring point sits at distance exactly `a` from its own station
        d = [norm(tube.disc.quad_pts[(i-1)*Q_cs + j] - centres[i]) for i in 1:Q, j in 1:Q_cs]
        @test all(≈(a; atol=1e-12), d)

        # a ring is centred on its station and lies in the plane normal to the tangent
        for i in 1:Q
            ring = @view tube.disc.quad_pts[(i-1)*Q_cs+1 : i*Q_cs]
            @test sum(ring) / Q_cs ≈ centres[i] atol=1e-12
            τ = unit_tangent((i-1)/(Q-1), 0.0, m)
            @test all(p -> abs(dot(p - centres[i], τ)) < 1e-12, ring)
        end
    end

    @testset "one region, and it covers the whole quadrature" begin
        tube = Part(LineTubeFlagellum(beating(), 0.02; Q_cs=12), 21, 101; eps=0.002)
        @test length(tube.disc.force_part_ranges) == 1
        @test length(tube.disc.quad_part_ranges) == 1
        @test get_eps(tube) == 0.002

        prob = SwimmingProblem(MicroSwimmer([tube]))
        @test length(prob.disc.quad_part_ranges) == 1
        @test last(last(prob.disc.quad_part_ranges)) == nq(prob.disc)

        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "nearest maps each ring to its own station" begin
        # A ring point is offset from the axis perpendicular to τ while force points are
        # separated along τ, so its distance to force point f is ~sqrt(d² + a²) with d the
        # axial separation — minimised by the smallest d whatever the radius. Worth pinning,
        # because it is what lets the radius be varied without touching the partition.
        Q, Q_cs = 61, 8
        tube = Part(LineTubeFlagellum(straight(10.0), 0.1; Q_cs=Q_cs), 21, Q; eps=0.01)
        for i in 1:Q
            ring = tube.disc.nearest[(i-1)*Q_cs+1 : i*Q_cs]
            @test all(==(first(ring)), ring)      # a ring is never split across force points
        end
        @test issorted(tube.disc.nearest)         # and the map is monotone in station index
    end

    @testset "a → 0 collapses onto the untubed model" begin
        m = beating()
        ref = SwimmingProblem(MicroSwimmer([Part(m, 21, 101; eps=0.02)]))
        solve_problem!(ref)

        errs = Float64[]
        for a in (1e-2, 1e-3, 1e-4)
            p = SwimmingProblem(MicroSwimmer([
                Part(LineTubeFlagellum(m, a; Q_cs=12), 21, 101; eps=0.02)]))
            solve_problem!(p)
            push!(errs, norm(get_U(p) - get_U(ref)))
        end
        @test issorted(errs, rev=true)     # shrinking the radius shrinks the difference
        @test last(errs) < 1e-6            # and it goes to zero, not to some offset

        # and it does so at O(a²), not O(a): the terms linear in the ring offset cancel
        # around the ring by symmetry, leaving the curvature of the kernel as the leading
        # correction. A decade in `a` should buy two decades of agreement.
        @test all(r -> 50 < r < 200, errs[1:end-1] ./ errs[2:end])
    end

    @testset "the ring quadrature converges spectrally in Q_cs" begin
        # The collocation point is at distance exactly `a` from every point of its own ring,
        # so the ring integrand is smooth and periodic and the equally-weighted trapezoid rule
        # over it converges spectrally. This is why Q_cs is a dozen and not 2πa/eps.
        #
        # Measured on a CURVED filament deliberately: on a straight one every ring is
        # congruent and the symmetry makes even Q_cs=4 exact, which would flatter the rule.
        m = curved(10.0)
        ref = drag(Part(LineTubeFlagellum(m, 0.1; Q_cs=48), 61, 301; eps=0.01))
        err(q) = abs(drag(Part(LineTubeFlagellum(m, 0.1; Q_cs=q), 61, 301; eps=0.01))[1]
                     - ref[1]) / ref[1]
        @test err(4) < 1e-4
        @test err(6) < 1e-7
        @test err(8) < 1e-10
        @test err(DEFAULT_Q_CS) < 1e-12
    end

    @testset "the answer is converged in eps once eps ≲ a/10" begin
        # eps and the radius are independent knobs by design, but the tube is only doing its
        # job when the ring, not the blob, sets the transverse length scale.
        m = straight(10.0)
        c(e) = drag(Part(LineTubeFlagellum(m, 0.1; Q_cs=12), 61, 305; eps=e))[1]
        c10, c100 = c(0.01), c(0.001)
        @test abs(c10 - c100) / c100 < 5e-3
    end

    @testset "a tube of radius a is not a blob of radius a" begin
        # The finding the whole wrapper exists to establish. `eps` on a flagellum has always
        # doubled as the filament radius; an actual filament of radius a is measurably
        # stiffer than that — 9% more translating drag at a=0.02 rising to 16% at a=0.2 —
        # and the gap widens as the filament fattens. So the two are NOT interchangeable,
        # and a study that varies eps is not a study that varies the radius.
        #
        # N=41 gives hf=0.238, which keeps hf > a for every radius below. Sampling this at a
        # finer mesh would put the fatter radii on the wrong side of the LineTube's
        # resolution limit and the comparison would be measuring that instead.
        m, N, Q = straight(10.0), 41, 205
        gaps = Float64[]
        for a in (0.02, 0.05, 0.1, 0.2)
            tube = drag(Part(LineTubeFlagellum(m, a; Q_cs=12), N, Q; eps=a/10))[1]
            blob = drag(Part(m, N, Q; eps=a))[1]
            @test tube > blob                       # the tube is always the stiffer of the two
            push!(gaps, (tube - blob) / blob)
        end
        @test all(>(0.05), gaps)                    # by at least 5% everywhere
        @test issorted(gaps)                        # and increasingly so with radius
    end

    @testset "drag anisotropy trends towards the slender-body limit" begin
        # Convention-free check: whatever constant sits inside the log, slender-body theory
        # has C_perp/C_par -> 2 from below as L/a grows. Testing the trend rather than an
        # absolute value keeps this independent of the half-length/full-length ambiguity that
        # makes published prefactors disagree.
        m, N, Q = straight(10.0), 41, 205        # hf = 0.238 > every a below
        ratios = Float64[]
        for a in (0.2, 0.1, 0.05, 0.02)          # L/a from 50 to 500
            par, perp = drag(Part(LineTubeFlagellum(m, a; Q_cs=12), N, Q; eps=a/10))
            push!(ratios, perp / par)
        end
        @test all(r -> 1 < r < 2, ratios)
        @test issorted(ratios)
    end

    @testset "refining a LineTube past hf ~ a is warned about" begin
        # Collocating on the axis means a ring is `a` from its own force point and
        # sqrt(hf² + a²) from its neighbours'. Once hf < a those are nearly equal, the
        # nearest-neighbour partition stops meaning anything, and the solution degrades under
        # refinement instead of converging. The warning is the only thing standing between a
        # user and a mesh-refinement study that silently gets worse.
        m = straight(10.0)                          # hf = 10/(N+1)
        @test_logs (:warn, r"below its tube radius") match_mode=:any Part(
            LineTubeFlagellum(m, 0.1; Q_cs=8), 241, 601; eps=0.01)
        # comfortably resolved: no warning
        @test_logs Part(LineTubeFlagellum(m, 0.1; Q_cs=8), 41, 201; eps=0.01)
    end

    @testset "eps above radius/2 is warned about" begin
        @test_logs (:warn, r"blob is wider than the") match_mode=:any Part(
            LineTubeFlagellum(beating(), 0.02; Q_cs=8), 21, 101; eps=0.05)
    end

    @testset "sizing and forwarding" begin
        L = 10.0
        m = LineTubeFlagellum(straight(L), 0.05; Q_cs=12)
        @test arclength(m) == L
        @test MicroSwimmers.radius(m) == 0.05
        @test surface_area(m) ≈ 2π * 0.05 * L

        # Part(model) sizes the *stations* by the inherited arclength rule; the ring count
        # comes from the model. Same hf and hq as the untubed model, times Q_cs.
        p     = Part(m; eps=0.005)
        plain = Part(straight(L); eps=0.005)
        @test nf(p.disc) == nf(plain.disc)
        @test nq(p.disc) == nq(plain.disc) * 12
        @test L / (nf(p.disc) + 1) ≈ DEFAULT_HF rtol=0.1
    end

    # -----------------------------------------------------------------------------------
    #  SurfaceTubeFlagellum — the reference case: force points on the surface too
    # -----------------------------------------------------------------------------------

    @testset "surface tube: cloud shape and partition" begin
        m = beating()
        N_cs, Q_cs, a = 8, 16, 0.02
        st = Part(SurfaceTubeFlagellum(m, a; N_cs=N_cs, Q_cs=Q_cs), 21, 101; eps=0.01)

        @test nf(st.disc) == 21 * N_cs
        @test nq(st.disc) == 101 * Q_cs
        @test length(st.disc.quad_part_ranges) == 1

        # every force point is at distance `a` from its own station on the centreline
        centres = Vector{SVector{3,Float64}}(undef, 21)
        m(centres, 0.0; include_endpoints=false)
        d = [norm(st.disc.force_pts[(i-1)*N_cs + j] - centres[i]) for i in 1:21, j in 1:N_cs]
        @test all(≈(a; atol=1e-12), d)

        prob = SwimmingProblem(MicroSwimmer([st]))
        solve_problem!(prob)
        F, T = total_force_and_torque(prob)
        @test norm(F) < 1e-8
        @test norm(T) < 1e-8
    end

    @testset "surface tube: the frame is rotation-minimising" begin
        # A per-station Gram-Schmidt frame jumps when the tangent crosses |τ_x| = 0.9, which
        # would put a seam of displaced unknowns down the tube. Parallel transport does not.
        # Check it directly: consecutive rings should be nearly aligned, so the azimuthal
        # ray from each station to its own j-th point turns only as fast as the centreline
        # does — a seam would show as one huge jump against a background of small ones.
        m = ThreeDimensionalFlagellum(10.0, 1.0, 0.6, 2.0, 6.0, 0.0,
                                            1.0, 0.6, 2.0, 6.0, 0.0, 0.0, π/2)
        N, N_cs, a = 81, 8, 0.05
        st = Part(SurfaceTubeFlagellum(m, a; N_cs=N_cs, Q_cs=16), N, 201; eps=0.02)

        centres = Vector{SVector{3,Float64}}(undef, N)
        m(centres, 0.0; include_endpoints=false)
        spokes = [normalize(st.disc.force_pts[(i-1)*N_cs + 1] - centres[i]) for i in 1:N]
        turns  = [acos(clamp(dot(spokes[i], spokes[i+1]), -1, 1)) for i in 1:N-1]
        # no single station rotates the frame by anything like the π-ish jump a
        # Gram-Schmidt branch flip would produce
        @test maximum(turns) < 0.2
        # and the tangent really does cross the branch, so this is a live test
        τs = [unit_tangent((i-1)/(N-1), 0.0, m) for i in 1:N]
        @test minimum(abs(τ[1]) for τ in τs) < 0.9 < maximum(abs(τ[1]) for τ in τs)
    end

    @testset "surface tube: ForwardDiff velocities match finite differences" begin
        # There is no local formula for the surface velocity — the transported normal at a
        # station depends on the one before it, so ∂n/∂t is a recurrence along the curve.
        # The Dual-valued march gets it exactly; check that against a central difference.
        m = beating()
        N_cs, a, t = 8, 0.02, 0.31
        st = Part(SurfaceTubeFlagellum(m, a; N_cs=N_cs, Q_cs=16), 21, 101; eps=0.01)
        update_boundary!(st, t)
        got = copy(st.disc.velocity)

        h = 1e-6
        update_boundary!(st, t + h); fwd = copy(st.disc.force_pts)
        update_boundary!(st, t - h); bwd = copy(st.disc.force_pts)
        fd = (fwd .- bwd) ./ (2h)

        @test maximum(norm.(got .- fd)) < 1e-7
        @test maximum(norm.(got)) > 1e-3        # and the velocities are not all zero
    end

    @testset "surface tube: agrees with the LineTube on drag" begin
        # The two discretise the same physical filament, so they should give the same drag to
        # within their discretisation error. This is the check that the LineTube's cheap
        # trick — collocating on the axis instead of the surface — is not distorting the
        # answer, which is the whole reason the surface variant exists.
        m, a = straight(10.0), 0.1
        line = drag(Part(LineTubeFlagellum(m, a; Q_cs=16), 41, 205; eps=a/10))
        surf = drag(Part(SurfaceTubeFlagellum(m, a; N_cs=8, Q_cs=16), 41, 205; eps=a/4))
        @test isapprox(surf[1], line[1]; rtol=0.15)
        @test isapprox(surf[2], line[2]; rtol=0.15)
        # and both are stiffer than a blob of the same nominal radius
        blob = drag(Part(m, 41, 205; eps=a))
        @test surf[1] > blob[1]
    end

    @testset "surface tube: lopsided cells are warned about" begin
        # Sized by the arclength rule alone a slender tube gets very long thin surface cells,
        # because the axial spacing comes from L and the ring spacing from a.
        @test_logs (:warn, r"aspect ratio") match_mode=:any Part(
            SurfaceTubeFlagellum(straight(10.0), 0.01; N_cs=8, Q_cs=16), 21, 101; eps=0.005)
    end

    @testset "a tube can still be asked for its centreline" begin
        # The wrapper forwards unit_tangent, so the generic FlagellumModel machinery — the
        # (s,t) grid method, the plotted centreline — keeps working on a tube.
        m    = beating()
        tube = LineTubeFlagellum(m, 0.05; Q_cs=12)
        a = Vector{SVector{3,Float64}}(undef, 31)
        b = Vector{SVector{3,Float64}}(undef, 31)
        m(a, 0.25; include_endpoints=true)
        tube(b, 0.25; include_endpoints=true)
        @test a == b
        @test unit_tangent(0.3, 0.25, tube) == unit_tangent(0.3, 0.25, m)
        @test unit_tangent_and_dt(0.3, 0.25, tube) == unit_tangent_and_dt(0.3, 0.25, m)
    end

    @testset "the tube deforms with time like any other flagellum" begin
        tube = Part(LineTubeFlagellum(beating(), 0.02; Q_cs=12), 21, 101; eps=0.002)
        before = copy(tube.disc.quad_pts)
        update_boundary!(tube, 0.3)
        @test tube.disc.quad_pts != before
        # and it still is a tube afterwards
        centres = Vector{SVector{3,Float64}}(undef, 101)
        beating()(centres, 0.3; include_endpoints=true)
        @test all(i -> norm(tube.disc.quad_pts[(i-1)*12 + 1] - centres[i]) ≈ 0.02, 1:101)
    end
end
