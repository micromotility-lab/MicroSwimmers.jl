abstract type Kernel end

# Whether assembly should thread the outer (force-point) loop for this kernel.
# Measured on real problems: threading wins for RegBlakelet (~4x, compute-bound
# per-call cost amortises the per-thread launch overhead) but loses for the much
# cheaper RegStokeslet (~5x slower — dominated by that overhead). Default false.
threaded_assembly(::Kernel) = false

struct RegStokeslet{T} <: Kernel
    eps::T
end

# Two forms throughout: the 2-argument one regularises at the kernel's own `eps`, the
# 3-argument one takes it per call. Assembly uses the latter so each region of the boundary
# can carry its own regularisation — the blob is centred on the source point Xj, so the eps
# that applies is the one belonging to whichever region Xj came from.
@inline (k::RegStokeslet)(xi, Xj) = k(xi, Xj, k.eps)

@inline function (k::RegStokeslet)(xi, Xj, eps)
    R     = xi - Xj
    rsqr  = dot(R, R)
    eps2  = eps^2
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

struct RegBlakelet{T} <: Kernel
    eps::T
end

threaded_assembly(::RegBlakelet) = true

@inline (k::RegBlakelet)(xi, Xj) = k(xi, Xj, k.eps)

@inline function (k::RegBlakelet)(xi, Xj, eps)
    eps2 = eps^2
    h    = Xj[3]

    stokeslet = RegStokeslet(eps)

    # -------- real-space regularised stokeslet --------
    S = stokeslet(xi, Xj, eps)

    # -------- image geometry (reflect source in z = 0) --------
    Y    = @SVector [Xj[1], Xj[2], -Xj[3]]
    Rimg = xi - Y
    X1, X2, X3 = Rimg

    Rsq  = dot(Rimg, Rimg)
    dist = sqrt(Rsq + eps2)
    iR   = inv(dist)
    iR3  = iR^3
    iR5  = iR^5

    Δ = @SMatrix [one(eps)  zero(eps) zero(eps);
                  zero(eps) one(eps)  zero(eps);
                  zero(eps) zero(eps) -one(eps)]

    P = @SMatrix [X1*X1  X1*X2  -X1*X3;
                  X2*X1  X2*X2  -X2*X3;
                  X3*X1  X3*X2  -X3*X3]

    # -------- image regularised stokeslet (Smith: subtracted) --------
    # Rsq, iR3 already computed above for the higher-order terms — reuse them
    # instead of recomputing dot(Rimg,Rimg) and a second sqrt/inv via stokeslet(xi, Y)
    diag_img = (Rsq + 2eps2) * iR3
    Simg = diag_img * I + (iR3 .* Rimg) * Rimg'

    # -------- blob (finite-eps) term --------
    phi = 3eps2 * iR5
    BT  = (-2h^2 * phi) .* Δ

    # -------- potential source dipole --------
    PD = (2h^2) .* (iR3 .* Δ .- 3iR5 .* P)

    # -------- regularised stokes dipole --------
    val = Rsq + 4eps2
    a1  =  X1 * val * iR5
    a2  =  X2 * val * iR5
    a3  = -X3 * val * iR5
    A = @SMatrix [zero(eps) zero(eps) zero(eps);
                  zero(eps) zero(eps) zero(eps);
                  a1        a2        a3]

    Dterm = (X3 * iR3) .* Δ

    C = @SMatrix [zero(eps) zero(eps) X1*iR3;
                  zero(eps) zero(eps) X2*iR3;
                  zero(eps) zero(eps) X3*iR3]

    SD = (2h) .* (A .- Dterm .+ C .+ (3iR5*X3) .* P)

    # -------- rotlet difference term --------
    M = @SMatrix [-X3       zero(eps) zero(eps);
                  zero(eps) -X3       zero(eps);
                  X1        X2        zero(eps)]
    RD = (-(6h*eps2*iR5)) .* M

    S .- Simg .+ BT .+ PD .+ SD .+ RD
end


# function regularised_blakelet!(B, T, x, X; eps=1e-6)
#     @assert length(x) == 3
#     @assert length(X) == 3

#     # Clear
#     fill!(B, zero(eltype(B)))

#     # -------- real-space regularised stokeslet --------
#     R = x .- X
#     regularised_stokeslet!(B, R; eps=eps)   # your existing kernel

#     # -------- image geometry (reflect source in z=0) --------
#     Y = @SVector [X[1], X[2], -X[3]]
#     Rimg = x .- Y
#     X1, X2, X3 = Rimg
#     h = X[3]

#     Rsq  = X1^2 + X2^2 + X3^2
#     dist = sqrt(Rsq + eps^2)
#     iR   = inv(dist)
#     iR3  = iR^3
#     iR5  = iR^5

#     # Convenient matrices
#     Δ = @SMatrix [1.0 0.0 0.0;
#                   0.0 1.0 0.0;
#                   0.0 0.0 -1.0]

#     # P = [X1X1 X1X2 -X1X3; X2X1 X2X2 -X2X3; X3X1 X3X2 -X3X3]
#     P = @SMatrix [X1*X1  X1*X2  -X1*X3;
#                   X2*X1  X2*X2  -X2*X3;
#                   X3*X1  X3*X2  -X3*X3]

#     # -------- image regularised stokeslet (Smith: just negative) --------
#     # Rsq, iR3 already computed above — avoid a second sqrt/dot inside regularised_stokeslet!
#     diag_img = Rsq + 2*eps^2
#     @inbounds for i in 1:3, j in 1:3
#         T[i,j] = diag_img * (i == j) + Rimg[i] * Rimg[j]
#         T[i,j] *= iR3
#     end
#     B .-= T

#     # -------- higher order terms (Smith) --------

#     # Blob term: BT = -2 h^2 * kron(Δ, phi)  with phi = 3 eps^2 iR^5
#     phi = 3 * eps^2 * iR5
#     B .+= (-2h^2 * phi) .* Δ

#     # Potential source dipole:
#     # PD = 2 h^2 * ( Δ*iR3 - 3*iR5*P )
#     B .+= (2h^2) .* (iR3 .* Δ .- 3*iR5 .* P)

#     # Regularised stokes dipole:
#     # SD = 2 h * ( A - Δ*(X3*iR3) + C + 3*iR5*X3*P )
#     #
#     # A: only row 3 nonzero
#     val = Rsq + 4eps^2
#     a1 = X1 * val * iR5
#     a2 = X2 * val * iR5
#     a3 = -X3 * val * iR5
#     A = @SMatrix [0.0 0.0 0.0;
#                   0.0 0.0 0.0;
#                   a1  a2  a3]

#     # -Δ*(X3*iR3)  (scalar times Δ)
#     Dterm = (X3 * iR3) .* Δ

#     # C: only column 3 nonzero  (THIS is what your Julia version was missing)
#     C = @SMatrix [0.0 0.0 X1*iR3;
#                   0.0 0.0 X2*iR3;
#                   0.0 0.0 X3*iR3]

#     SD = (2h) .* (A .- Dterm .+ C .+ (3*iR5*X3) .* P)
#     B .+= SD

#     # Rotlet difference term:
#     # RD = -(6 h eps^2 iR^5) * ( [0;0;X1 X2 X3] - X3*I )
#     # For 3×3: [[-X3,0,0],[0,-X3,0],[X1,X2,0]]
#     M = @SMatrix [-X3  0.0 0.0;
#                   0.0 -X3 0.0;
#                   X1   X2  0.0]
#     RD = (-(6h*eps^2*iR5)) .* M
#     B .+= RD
# end



# function resistance_matrix!(
#     A::AbstractMatrix{T},
#     force_pts::AbstractMatrix{T},
#     quad_pts::AbstractMatrix{T},
#     nearest::AbstractVector{Int},
#     eps::T;
#     μ::T=one(T),
#     wall::Bool=false
# ) where {T <: Number}
#     fill!(A, zero(T))
#     S = MMatrix{3,3,T}(undef)
#     S2 = MMatrix{3,3,T}(undef)

#     for i in axes(force_pts, 2)
#         xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
#         for j in axes(quad_pts, 2)
#             Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]] 
#             if wall
#                 regularised_blakelet!(S, S2, xi, Xj;  eps=eps)
#             else
#                 R = xi - Xj
#                 regularised_stokeslet!(S, R; eps=eps)
#             end
#             # R = xi - Xj
#             # regularised_stokeslet!(S, R; eps=eps)
#             # stokeslet!(S, R)

#             n = nearest[j]
#             @inbounds for p in 1:3, q in 1:3
#                 A[3i-3+p, 3n-3+q] -= S[p,q]
#             end
#         end
#     end
#     A ./= (-T(8) * T(π) * μ)
# end

# function resistance_matrix!(
#     A::AbstractMatrix{T},
#     force_pts::AbstractVector{T},
#     quad_pts::AbstractMatrix{T},
#     nearest::AbstractVector{Int},
#     eps::T;
#     μ::T=one(T),
#     wall::Bool=false
# ) where {T <: Number}   
#     resistance_matrix!(A, reshape(force_pts, 3, 1), quad_pts, nearest, eps; μ=μ, wall=wall)
# end

# function swimming_matrix!(
#     A::Matrix{T},
#     x0::SVector{3,T},
#     force_pts::Matrix{T},
#     quad_pts::Matrix{T},
#     nearest::Vector{Int},
#     eps::T;
#     μ::T=one(T),
#     wall::Bool=false
# ) where {T <: Number}
#     fill!(A, zero(T))

#     S = MMatrix{3,3,T}(undef)
#     S2 = MMatrix{3,3,T}(undef)
#     diffvec = MVector{3,T}(undef)

#     for i in axes(force_pts, 2)
#         xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
#         for j in axes(quad_pts, 2)
#             Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]]

#             n = nearest[j]
#             if wall
#                 regularised_blakelet!(S, S2, xi, Xj;  eps=eps)
#             else
#                 R = xi - Xj
#                 regularised_stokeslet!(S, R; eps=eps)
#             end
#             @inbounds for p in 1:3, q in 1:3
#                 A[3i-3+p, 3n-3+q] -= S[p,q]
#             end
#         end

#         @inbounds for p in 1:3
#             A[3i-3+p, end-6+p] = -one(T)
#             diffvec[p] = force_pts[p,i] - x0[p]
#         end

#         K = skew_symmetric_static(diffvec)
#         @inbounds for p in 1:3, q in 1:3
#             A[3i-3+p, end-3+q] = K[p,q]
#         end
#     end

#     nf = length(force_pts)
#     @inbounds A[1:nf, 1:nf] ./= (-T(8) * T(π) * μ)

#     for j in axes(quad_pts, 2)
#         n = nearest[j]

#         @inbounds for d in 1:3
#             A[end-6+d, 3n-3+d] += one(T)
#         end

#         rq = SVector{3,T}(quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]) - x0

#         Kq = skew_symmetric_static(rq)
#         @inbounds for p in 1:3, q in 1:3
#             A[end-3+p, 3n-3+q] += Kq[p,q]
#         end
#     end
# end

# function swimming_matrix!(
#     A::Matrix{T},
#     x0::SVector{3,T},
#     points::NearestDiscretisation,
#     eps::T;
#     μ::T=one(T),
#     wall::Bool=false
# ) where {T <: Number}
#     @unpack force_pts, quad_pts, nearest = points
#     swimming_matrix!(A, x0, force_pts, quad_pts, nearest, eps; μ=μ, wall=wall)
# end

# Quadrature weight for point q. `nothing` means the weights are absorbed into the force
# unknowns (the default): `true` promotes to a literal 1.0, which LLVM folds out of the
# inner loops, so the absorbed path costs exactly what it did before weights existed.
@inline quad_weight(::Nothing, ::Int)          = true
@inline quad_weight(w::AbstractVector, q::Int) = @inbounds w[q]

# The regularisation partition: `ranges[j]` is a contiguous block of quadrature points all
# regularised at `epsv[j]`. `nothing` means one implicit region at the kernel's own eps, which
# is what every call site did before eps became per-region.
@inline eps_regions(::Nothing, ::Nothing, kernel::Kernel, Q::Int) = ((1:Q, kernel.eps),)

@inline function eps_regions(ranges::AbstractVector, epsv::AbstractVector, ::Kernel, ::Int)
    # zip truncates to the shorter of the two, which would silently drop whole blocks of
    # quadrature points from the assembly. Once per assemble! call, so it costs nothing.
    length(ranges) == length(epsv) || throw(ArgumentError(
        "discretisation has $(length(ranges)) regions but $(length(epsv)) regularisation " *
        "parameters — call set_eps! after changing the region partition"))
    zip(ranges, epsv)
end

assemble!(A, disc::NearestDiscretisation, kernel; μ=one(eltype(A))) =
    assemble!(A, disc.force_pts, disc.quad_pts, disc.nearest, disc.quad_wts, kernel;
              μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

# NystromDiscretisation collocates force and quadrature points, so there are no separate
# quadrature weights to apply -- it always takes the absorbed path.
assemble!(A, disc::NystromDiscretisation, kernel; μ=one(eltype(A))) =
    assemble!(A, disc.force_pts, disc.force_pts, 1:nf(disc), nothing, kernel;
              μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

# Regions outside, points inside: eps is constant across a region, so it lands in a register
# once per region rather than being reloaded per quadrature point.
@inline function accumulate_row!(A, row, xm, quad_pts, nearest, wts, regions, kernel)
    for (rng, ε) in regions
        for q in rng
            S   = kernel(xm, @inbounds(quad_pts[q]), ε) * quad_weight(wts, q)
            col = 3*(@inbounds nearest[q]) - 2
            @inbounds for b in 1:3, a in 1:3
                A[row+a-1, col+b-1] += S[a,b]
            end
        end
    end
end

# no `wts` argument => quadrature weights absorbed into the force unknowns
assemble!(A, force_pts, quad_pts, nearest, kernel; μ=one(eltype(A)), kwargs...) =
    assemble!(A, force_pts, quad_pts, nearest, nothing, kernel; μ=μ, kwargs...)

function assemble!(A, force_pts, quad_pts, nearest, wts, kernel;
                   μ=one(eltype(A)), ranges=nothing, quad_eps=nothing)
    fill!(A, zero(eltype(A)))
    Nf = length(force_pts)
    regions = eps_regions(ranges, quad_eps, kernel, length(quad_pts))
    # m (force points) is the outer loop: each m owns a disjoint row block of
    # A, so threads never race on the same A[row,col] accumulator when threaded.
    if threaded_assembly(kernel)
        @batch for m in 1:Nf
            accumulate_row!(A, 3m - 2, force_pts[m], quad_pts, nearest, wts, regions, kernel)
        end
    else
        for m in 1:Nf
            accumulate_row!(A, 3m - 2, force_pts[m], quad_pts, nearest, wts, regions, kernel)
        end
    end
    A .*= inv(8π*μ)
end

assemble_swimming!(A, x0, disc::NearestDiscretisation, kernel; μ=one(eltype(A))) =
    assemble_swimming!(A, x0, disc.force_pts, disc.quad_pts, disc.nearest, disc.quad_wts, kernel;
                       μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

assemble_swimming!(A, x0, disc::NystromDiscretisation, kernel; μ=one(eltype(A))) =
    assemble_swimming!(A, x0, disc.force_pts, disc.force_pts, 1:nf(disc), nothing, kernel;
                       μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

assemble_swimming!(A::AbstractMatrix, x0::SVector{3},
                   force_pts, quad_pts, nearest, kernel; μ=one(eltype(A)), kwargs...) =
    assemble_swimming!(A, x0, force_pts, quad_pts, nearest, nothing, kernel; μ=μ, kwargs...)

function assemble_swimming!(A::AbstractMatrix, x0::SVector{3},
                             force_pts, quad_pts, nearest, wts, kernel;
                             μ=one(eltype(A)), ranges=nothing, quad_eps=nothing)
    T  = eltype(A)
    N  = length(force_pts)
    N3 = 3N
    fill!(A, zero(T))

    # BEM block: A[3m-2:3m, 3nearest[q]-2:3nearest[q]] += kernel(xm, yq)
    assemble!(view(A, 1:N3, 1:N3), force_pts, quad_pts, nearest, wts, kernel;
              μ=μ, ranges=ranges, quad_eps=quad_eps)

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
    # These rows must carry the same quadrature weight as the BEM block above. The net
    # force is Σ_q w_q f_{nearest[q]}, so weighting one and not the other would make
    # "force-free" mean something different from what total_force reports.
    @inbounds for (q, yq) in enumerate(quad_pts)
        n = nearest[q]
        w = quad_weight(wts, q)
        for d in 1:3
            A[N3+d, 3n-3+d] += w
        end
        rq = yq - x0
        Kq = skew_symmetric_static(rq)
        for p in 1:3, r in 1:3
            A[N3+3+p, 3*(n-1)+r] += w * Kq[p,r]
        end
    end
end

###########################################################################################
### Matrix-free swimming matvec ###########################################################
###########################################################################################
#
# y = A*x for the same (3N+6)×(3N+6) system that assemble_swimming! builds, without ever
# forming A. Used by the hybrid dense-LU/GMRES solver in problems.jl: once a dense
# factorisation goes stale, GMRES needs a matvec (not a fresh assembly) on every iteration —
# assembling the dense matrix on every iteration would cost the same as just refactorising it.

@inline function accumulate_row_matvec!(y, row, xm, quad_pts, nearest, wts, regions, kernel, x, invfactor)
    T   = eltype(y)
    acc = zero(SVector{3,T})
    for (rng, ε) in regions
        for q in rng
            S   = kernel(xm, @inbounds(quad_pts[q]), ε) * quad_weight(wts, q)
            col = 3*(@inbounds nearest[q]) - 2
            xv  = SVector{3,T}(x[col], x[col+1], x[col+2])
            acc = acc + S * xv
        end
    end
    acc = acc * invfactor
    @inbounds for a in 1:3
        y[row+a-1] += acc[a]
    end
end

mul_swimming!(y, x, x0, disc::NearestDiscretisation, kernel; μ=one(eltype(y))) =
    mul_swimming!(y, x, x0, disc.force_pts, disc.quad_pts, disc.nearest, disc.quad_wts, kernel;
                  μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

mul_swimming!(y, x, x0, disc::NystromDiscretisation, kernel; μ=one(eltype(y))) =
    mul_swimming!(y, x, x0, disc.force_pts, disc.force_pts, 1:nf(disc), nothing, kernel;
                  μ=μ, ranges=disc.quad_part_ranges, quad_eps=disc.quad_eps)

mul_swimming!(y::AbstractVector, x::AbstractVector, x0::SVector{3},
              force_pts, quad_pts, nearest, kernel; μ=one(eltype(y)), kwargs...) =
    mul_swimming!(y, x, x0, force_pts, quad_pts, nearest, nothing, kernel; μ=μ, kwargs...)

function mul_swimming!(y::AbstractVector, x::AbstractVector, x0::SVector{3},
                        force_pts, quad_pts, nearest, wts, kernel;
                        μ=one(eltype(y)), ranges=nothing, quad_eps=nothing)
    T  = eltype(y)
    Nf = length(force_pts)
    N3 = 3Nf
    fill!(y, zero(T))
    invfactor = inv(8π*μ)
    regions   = eps_regions(ranges, quad_eps, kernel, length(quad_pts))

    # BEM block: rows 1:N3, cols 1:N3
    if threaded_assembly(kernel)
        @batch for m in 1:Nf
            accumulate_row_matvec!(y, 3m - 2, force_pts[m], quad_pts, nearest, wts, regions, kernel, x, invfactor)
        end
    else
        for m in 1:Nf
            accumulate_row_matvec!(y, 3m - 2, force_pts[m], quad_pts, nearest, wts, regions, kernel, x, invfactor)
        end
    end

    # Rigid body columns' contribution to rows 1:N3: -I*U + K(xm - x0)*Ω
    U = SVector{3,T}(x[N3+1], x[N3+2], x[N3+3])
    Ω = SVector{3,T}(x[N3+4], x[N3+5], x[N3+6])
    @inbounds for (m, xm) in enumerate(force_pts)
        row = 3m - 2
        K = skew_symmetric_static(xm - x0)
        contrib = K * Ω - U
        for p in 1:3
            y[row+p-1] += contrib[p]
        end
    end

    # Force-free / torque-free rows N3+1:N3+6
    Facc = zero(SVector{3,T})
    Tacc = zero(SVector{3,T})
    @inbounds for (q, yq) in enumerate(quad_pts)
        n   = nearest[q]
        col = 3n - 2
        w   = quad_weight(wts, q)
        xv  = SVector{3,T}(x[col], x[col+1], x[col+2])
        Facc = Facc + w * xv
        Kq   = skew_symmetric_static(yq - x0)
        Tacc = Tacc + w * (Kq * xv)
    end
    @inbounds for d in 1:3
        y[N3+d]   += Facc[d]
        y[N3+3+d] += Tacc[d]
    end
    y
end