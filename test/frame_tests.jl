@testset "frame.jl" begin

    Random.seed!(3)

    function random_frame()
        loc = SVector{3}(randn(3))
        axis = normalize(randn(3))
        angle = 2π * rand() - π
        Frame(loc, rotation_matrix(axis, angle))
    end

    @testset "identity frame" begin
        F0 = Frame{Float64}()
        @test F0.location ≈ zeros(3) atol=1e-14
        @test F0.orientation ≈ I atol=1e-14

        X = SVector{3}(1.0, 2.0, 3.0)
        @test F0(X) ≈ X atol=1e-14
    end

    @testset "composition and inverse" begin
        for _ in 1:10
            P = random_frame()
            C = random_frame()
            X = SVector{3}(randn(3))

            # composition matches direct matrix/vector composition
            PC = P * C
            @test PC.location ≈ P.location + P.orientation * C.location atol=1e-10
            @test PC.orientation ≈ P.orientation * C.orientation atol=1e-10

            # action of the composed frame matches applying child then parent
            @test PC(X) ≈ P(C(X)) atol=1e-10

            # F * inv(F) is the identity frame
            Finv = inv(P)
            comp = P * Finv
            @test comp.location ≈ zeros(3) atol=1e-8
            @test comp.orientation ≈ I atol=1e-10

            # applying then un-applying returns the original point
            @test Finv(P(X)) ≈ X atol=1e-8
        end
    end

end
