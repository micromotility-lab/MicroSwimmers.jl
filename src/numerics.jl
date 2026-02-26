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


# function regularised_blakelet!(B, T, x, X; eps=1e-6)
#     @assert length(x) == 3 "x must be a 3D vector"
#     @assert length(X) == 3 "X must be a 3D vector"

#     @info "" eps

#     # Primary stokeslet
#     R = x .- X
#     regularised_stokeslet!(B, R; eps=eps)
    
#     Delta = Diagonal([1., 1., -1.])  
#     # Image stokeslet
#     Rimg     = x .- [X[1]; X[2]; -X[3]]
#     regularised_stokeslet!(T, Rimg; eps=eps)
#     B .-= T
    
#     # Some local shortcuts
#     h   = X[3]          # height of source point above z=0
#     R2  = sum(Rimg.^2)
#     dist = sqrt(R2 + eps^2)
#     iReps   = 1/dist
#     iReps3  = iReps^3
#     iReps5  = iReps^5
#     X1, X2, X3 = Rimg

#     # 1) Blob term:    -2*h^2 * Delta * (3 eps^2 iReps^5)
#     phi = 3eps^2 * iReps5
#     B .+= -2h^2 .* (phi *Delta)

#     # 2) Potential source dipole
#     #    2*h^2 [Delta*iReps3 - 3 iReps5 * (X1X1, X1X2, -X1X3; ...)]
#     M = Rimg .* Rimg'
#     M[:,3] .*= -1

#     B .+= 2h^2 .*(iReps3 * Delta .- 3iReps5 .* M)

#     # 3) Stokes dipole
#     #    2*h * [ (complicated expression) ]
#     # For brevity, we define a small local helper:

#     val = R2 + 4eps^2
#     T .= [
#         0 0 0;
#         0 0 0;
#         X1*val*iReps5   X2*val*iReps5   -X3*val*iReps5
#     ]
#     T .-= (X3 * iReps3) * Delta
#     T[3,:] .+= [X1*iReps3;  X2*iReps3;  X3*iReps3]
#     T .+= 3.0*iReps5*X3*M
#     T .*= 2h
#     B .+= T

#     # 4) Rotlet difference
#     #    -(6 * h * eps^2 * iReps5)*( [0;0; X1 X2 X3] - diag(X3,X3,X3) )
#     T .= [
#         0 0 0;
#         0 0 0;
#         X1 X2 X3
#     ]
#     T[diagind(T)] .-= X3
#     T *= -6*h*eps^2*iReps5
#     B .+= T
# end

# function regularised_blakelet!(B::StaticMatrix, T::StaticMatrix, x::StaticVector, x₀::StaticVector; eps=1e-6)
#     h = x₀[3]
#     r = x - x₀
#     r_norm² = norm(r)^2 + eps^2
#     r_norm = sqrt(r_norm²)
    
#     # Image configuration
#     x₀_image = SVector(x₀[1], x₀[2], -x₀[3])
#     r_image = x - x₀_image
#     r_image_norm² = norm(r_image)^2 + eps^2
#     r_image_norm = sqrt(r_image_norm²)
    
#     I₃ = SMatrix{3,3}(1.0I)
    
#     # Regularized Stokeslet
#     G_reg = ((r_norm² + 2eps^2) * I₃ + (r * r')) / (r_norm^3)
    
#     # Regularized image Stokeslet
#     G_image_reg = ((r_image_norm² + 2eps^2) * I₃ + (r_image * r_image')) / (r_image_norm^3)
    
#     # Regularized doublet terms
#     G_doublet_reg = -2h * (
#         (r_image_norm² + 2eps^2) * I₃ - 3 * (r_image * r_image')
#     ) / (r_image_norm^5)
    
#     e₃ = SVector(0.0, 0.0, 1.0)
#     G_source_doublet_reg = 2h * (
#         (e₃ * r_image') + (r_image * e₃') - 
#         2 * (e₃ ⋅ r_image) * (r_image * r_image') / r_image_norm²
#     ) / (r_image_norm^3)
    
#     B .= G_reg - G_image_reg + G_doublet_reg + G_source_doublet_reg 
# end

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
            
            # @inbounds for k in 1:3
            #     diffvec[k] = force_pts[k,i] - quad_pts[k,j]
            # end
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
    fill!(A, zero(T))

    S = MMatrix{3,3,T}(undef)
    S2 = MMatrix{3,3,T}(undef)
    diffvec = MVector{3,T}(undef)

    for i in axes(force_pts, 2)
        xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
        for j in axes(quad_pts, 2)
            Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]] 
            
            # @inbounds for k in 1:3
            #     diffvec[k] = force_pts[k,i] - quad_pts[k,j]
            # end
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
    points::NystromDiscretisation,
    eps::T;
    μ::T = one(T),
    wall::Bool = false
) where {T<:Number}

    pts = points.pts              # 3 × N
    N   = size(pts, 2)
    fill!(A, zero(T))

    S  = MMatrix{3,3,T}(undef)
    S2 = MMatrix{3,3,T}(undef)
    diffvec = MVector{3,T}(undef)

    # Collocation at i, integrate over j (same nodes)
    @inbounds for i in 1:N
        xi = @SVector [pts[1,i], pts[2,i], pts[3,i]]

        for j in 1:N
            Xj = @SVector [pts[1,j], pts[2,j], pts[3,j]]

            if wall
                regularised_blakelet!(S, S2, xi, Xj; eps=eps)
            else
                R = xi - Xj
                regularised_stokeslet!(S, R; eps=eps)
            end

            for p in 1:3, q in 1:3
                A[3i-3+p, 3j-3+q] -= S[p,q]
            end
        end

        # Translation columns
        for p in 1:3
            A[3i-3+p, end-6+p] = -one(T)
            diffvec[p] = pts[p,i] - x0[p]
        end

        # Rotation columns
        K = skew_symmetric_static(diffvec)
        for p in 1:3, q in 1:3
            A[3i-3+p, end-3+q] = K[p,q]
        end
    end

    # Scale the kernel block
    @inbounds A[1:3N, 1:3N] ./= (-T(8) * T(π) * μ)

    # Net force / net torque rows
    @inbounds for j in 1:N
        # ∑ f_j
        for d in 1:3
            A[end-6+d, 3j-3+d] += one(T)
        end

        # ∑ r_j × f_j
        rj = SVector{3,T}(pts[1,j], pts[2,j], pts[3,j]) - x0
        Kj = skew_symmetric_static(rj)
        for p in 1:3, q in 1:3
            A[end-3+p, 3j-3+q] += Kj[p,q]
        end
    end

    return A
end