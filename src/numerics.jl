abstract type Kernel end 

struct RegStokeslet{T} <: Kernel
    eps::T
end

@inline function (k::RegStokeslet)(xi, Xj)
    R     = xi - Xj
    rsqr  = dot(R, R)
    eps2 = k.eps^2
    denom = inv(sqrt(rsqr + eps2)^3)      
    diag  = (rsqr + 2eps2) * denom
    diag * I + (denom * R) * R'
end


function regularised_stokeslet!(S::StaticMatrix{3,3,T}, R::StaticVector{3,T}; eps::T=1e-6) where {T <: Number}
    rsqr = dot(R,R)
    diag = rsqr + 2eps^2
    denom = 1 / sqrt(rsqr + eps^2)^3

    @inbounds for i in 1:3, j in 1:3
        S[i,j] = diag * (i == j) + R[i] * R[j]
        S[i,j] *= denom
    end
end

function stokeslet!(S::StaticMatrix{3,3,T}, R::StaticVector{3,T}) where {T<:Number}
    rsqr = dot(R, R)

    # singular at r = 0
    if rsqr == zero(T)
        @inbounds for i in 1:3, j in 1:3
            S[i,j] = (i == j) ? T(Inf) : zero(T)
        end
        return
    end

    denom = inv(rsqr * sqrt(rsqr))  # 1 / r^3

    @inbounds for i in 1:3, j in 1:3
        S[i,j] = rsqr * (i == j) + R[i] * R[j]
        S[i,j] *= denom
    end
end


function regularised_blakelet!(B, T, x, X; eps=1e-6)
    @assert length(x) == 3
    @assert length(X) == 3

    # Clear
    fill!(B, zero(eltype(B)))

    # -------- real-space regularised stokeslet --------
    R = x .- X
    regularised_stokeslet!(B, R; eps=eps)   # your existing kernel

    # -------- image geometry (reflect source in z=0) --------
    Y = @SVector [X[1], X[2], -X[3]]
    Rimg = x .- Y
    X1, X2, X3 = Rimg
    h = X[3]

    Rsq  = X1^2 + X2^2 + X3^2
    dist = sqrt(Rsq + eps^2)
    iR   = inv(dist)
    iR3  = iR^3
    iR5  = iR^5

    # Convenient matrices
    Δ = @SMatrix [1.0 0.0 0.0;
                  0.0 1.0 0.0;
                  0.0 0.0 -1.0]

    # P = [X1X1 X1X2 -X1X3; X2X1 X2X2 -X2X3; X3X1 X3X2 -X3X3]
    P = @SMatrix [X1*X1  X1*X2  -X1*X3;
                  X2*X1  X2*X2  -X2*X3;
                  X3*X1  X3*X2  -X3*X3]

    # -------- image regularised stokeslet (Smith: just negative) --------
    regularised_stokeslet!(T, Rimg; eps=eps)
    B .-= T

    # -------- higher order terms (Smith) --------

    # Blob term: BT = -2 h^2 * kron(Δ, phi)  with phi = 3 eps^2 iR^5
    phi = 3 * eps^2 * iR5
    B .+= (-2h^2 * phi) .* Δ

    # Potential source dipole:
    # PD = 2 h^2 * ( Δ*iR3 - 3*iR5*P )
    B .+= (2h^2) .* (iR3 .* Δ .- 3*iR5 .* P)

    # Regularised stokes dipole:
    # SD = 2 h * ( A - Δ*(X3*iR3) + C + 3*iR5*X3*P )
    #
    # A: only row 3 nonzero
    val = Rsq + 4eps^2
    a1 = X1 * val * iR5
    a2 = X2 * val * iR5
    a3 = -X3 * val * iR5
    A = @SMatrix [0.0 0.0 0.0;
                  0.0 0.0 0.0;
                  a1  a2  a3]

    # -Δ*(X3*iR3)  (scalar times Δ)
    Dterm = (X3 * iR3) .* Δ

    # C: only column 3 nonzero  (THIS is what your Julia version was missing)
    C = @SMatrix [0.0 0.0 X1*iR3;
                  0.0 0.0 X2*iR3;
                  0.0 0.0 X3*iR3]

    SD = (2h) .* (A .- Dterm .+ C .+ (3*iR5*X3) .* P)
    B .+= SD

    # Rotlet difference term:
    # RD = -(6 h eps^2 iR^5) * ( [0;0;X1 X2 X3] - X3*I )
    # For 3×3: [[-X3,0,0],[0,-X3,0],[X1,X2,0]]
    M = @SMatrix [-X3  0.0 0.0;
                  0.0 -X3 0.0;
                  X1   X2  0.0]
    RD = (-(6h*eps^2*iR5)) .* M
    B .+= RD
end



function resistance_matrix!(
    A::AbstractMatrix{T},
    force_pts::AbstractMatrix{T},
    quad_pts::AbstractMatrix{T},
    nearest::AbstractVector{Int},
    eps::T;
    μ::T=one(T),
    wall::Bool=false
) where {T <: Number}
    fill!(A, zero(T))
    S = MMatrix{3,3,T}(undef)
    S2 = MMatrix{3,3,T}(undef)

    for i in axes(force_pts, 2)
        xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
        for j in axes(quad_pts, 2)
            Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]] 
            if wall
                regularised_blakelet!(S, S2, xi, Xj;  eps=eps)
            else
                R = xi - Xj
                regularised_stokeslet!(S, R; eps=eps)
            end
            # R = xi - Xj
            # regularised_stokeslet!(S, R; eps=eps)
            # stokeslet!(S, R)

            n = nearest[j]
            @inbounds for p in 1:3, q in 1:3
                A[3i-3+p, 3n-3+q] -= S[p,q]
            end
        end
    end
    A ./= (-T(8) * T(π) * μ)
end

function resistance_matrix!(
    A::AbstractMatrix{T},
    force_pts::AbstractVector{T},
    quad_pts::AbstractMatrix{T},
    nearest::AbstractVector{Int},
    eps::T;
    μ::T=one(T),
    wall::Bool=false
) where {T <: Number}   
    resistance_matrix!(A, reshape(force_pts, 3, 1), quad_pts, nearest, eps; μ=μ, wall=wall)
end

function swimming_matrix!(
    A::Matrix{T},
    x0::SVector{3,T},
    force_pts::Matrix{T},
    quad_pts::Matrix{T},
    nearest::Vector{Int},
    eps::T;
    μ::T=one(T),
    wall::Bool=false
) where {T <: Number}
    fill!(A, zero(T))

    S = MMatrix{3,3,T}(undef)
    S2 = MMatrix{3,3,T}(undef)
    diffvec = MVector{3,T}(undef)

    for i in axes(force_pts, 2)
        xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
        for j in axes(quad_pts, 2)
            Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]]

            n = nearest[j]
            if wall
                regularised_blakelet!(S, S2, xi, Xj;  eps=eps)
            else
                R = xi - Xj
                regularised_stokeslet!(S, R; eps=eps)
            end
            @inbounds for p in 1:3, q in 1:3
                A[3i-3+p, 3n-3+q] -= S[p,q]
            end
        end

        @inbounds for p in 1:3
            A[3i-3+p, end-6+p] = -one(T)
            diffvec[p] = force_pts[p,i] - x0[p]
        end

        K = skew_symmetric_static(diffvec)
        @inbounds for p in 1:3, q in 1:3
            A[3i-3+p, end-3+q] = K[p,q]
        end
    end

    nf = length(force_pts)
    @inbounds A[1:nf, 1:nf] ./= (-T(8) * T(π) * μ)

    for j in axes(quad_pts, 2)
        n = nearest[j]

        @inbounds for d in 1:3
            A[end-6+d, 3n-3+d] += one(T)
        end

        rq = SVector{3,T}(quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]) - x0

        Kq = skew_symmetric_static(rq)
        @inbounds for p in 1:3, q in 1:3
            A[end-3+p, 3n-3+q] += Kq[p,q]
        end
    end
end

function swimming_matrix!(
    A::Matrix{T},
    x0::SVector{3,T},
    points::NearestDiscretisation,
    eps::T;
    μ::T=one(T),
    wall::Bool=false
) where {T <: Number}
    @unpack force_pts, quad_pts, nearest = points
    swimming_matrix!(A, x0, force_pts, quad_pts, nearest, eps; μ=μ, wall=wall)
end

assemble!(A, disc::NearestDiscretisation, kernel; μ=one(eltype(A))) =
    assemble!(A, disc.force_pts, disc.quad_pts, disc.nearest, kernel; μ=μ)

assemble!(A, disc::NystromDiscretisation, kernel; μ=one(eltype(A))) =
    assemble!(A, disc.force_pts, disc.force_pts, 1:nf(disc), kernel; μ=μ)

function assemble!(A, force_pts, quad_pts, nearest, kernel; μ=one(eltype(A)))
    fill!(A, zero(eltype(A)))
    for (q, yq) in enumerate(quad_pts)
        col = 3*nearest[q] - 2
        for (m, xm) in enumerate(force_pts)
            S   = kernel(xm, yq)
            row = 3m - 2
            @inbounds for b in 1:3, a in 1:3
                A[row+a-1, col+b-1] += S[a,b]
            end
        end
    end
    A .*= inv(8π*μ)
end

# Dispatch wrappers for assemble_swimming! — same pattern as assemble!
assemble_swimming!(A, x0, disc::NearestDiscretisation, kernel; μ=one(eltype(A))) =
    assemble_swimming!(A, x0, disc.force_pts, disc.quad_pts, disc.nearest, kernel; μ=μ)

assemble_swimming!(A, x0, disc::NystromDiscretisation, kernel; μ=one(eltype(A))) =
    assemble_swimming!(A, x0, disc.force_pts, disc.force_pts, 1:nf(disc), kernel; μ=μ)

function assemble_swimming!(A::AbstractMatrix, x0::SVector{3},
                             force_pts, quad_pts, nearest, kernel; μ=one(eltype(A)))
    T  = eltype(A)
    N  = length(force_pts)
    N3 = 3N
    fill!(A, zero(T))

    # BEM block: A[3m-2:3m, 3nearest[q]-2:3nearest[q]] += kernel(xm, yq)
    for (q, yq) in enumerate(quad_pts)
        col = 3*nearest[q] - 2
        for (m, xm) in enumerate(force_pts)
            S   = kernel(xm, yq)
            row = 3m - 2
            @inbounds for b in 1:3, a in 1:3
                A[row+a-1, col+b-1] += S[a,b]
            end
        end
    end
    A[1:N3, 1:N3] .*= inv(8π * T(μ))

    # Rigid body columns: U (cols N3+1:N3+3) and Ω (cols N3+4:N3+6)
    @inbounds for (m, xm) in enumerate(force_pts)
        row = 3m - 2
        for p in 1:3
            A[row+p-1, N3+p] = -one(T)
        end
        K = skew_symmetric_static(xm - x0)
        for p in 1:3, q in 1:3
            A[row+p-1, N3+3+q] = K[p,q]
        end
    end

    # Force-free / torque-free rows (rows N3+1:N3+6)
    @inbounds for (q, yq) in enumerate(quad_pts)
        n = nearest[q]
        for d in 1:3
            A[N3+d, 3n-3+d] += one(T)
        end
        rq = yq - x0
        Kq = skew_symmetric_static(rq)
        for p in 1:3, r in 1:3
            A[N3+3+p, 3*(n-1)+r] += Kq[p,r]
        end
    end
end