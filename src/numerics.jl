function regularised_stokeslet!(S::StaticMatrix{3,3,T}, R::StaticVector{3,T}; eps::T=1e-6) where {T <: Number}
    rsqr = dot(R,R)
    diag = rsqr + 2eps^2
    denom = 1 / sqrt(rsqr + eps^2)^3

    @inbounds for i in 1:3, j in 1:3
        S[i,j] = diag * (i == j) + R[i] * R[j]
        S[i,j] *= denom
    end
end

function resistance_matrix!(
    A::AbstractMatrix{T},
    force_pts::AbstractMatrix{T},
    quad_pts::AbstractMatrix{T},
    nearest::AbstractVector{Int},
    eps::T;
    μ::T=one(T),
) where {T <: Number}
    fill!(A, zero(T))
    S = MMatrix{3,3,T}(undef)

    Threads.@threads for i in axes(force_pts, 2)
        xi = @SVector [force_pts[1,i], force_pts[2,i], force_pts[3,i]]
        for j in axes(quad_pts, 2)
            Xj = @SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]] 
            R = xi - Xj
            regularised_stokeslet!(S, R; eps=eps)

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
) where {T <: Number}   
    resistance_matrix!(A, reshape(force_pts, 3, 1), quad_pts, nearest, eps; μ=μ)
end

function swimming_matrix!(
    A::Matrix{T}, 
    x0::SVector{3,T}, 
    force_pts::Matrix{T}, 
    quad_pts::Matrix{T}, 
    nearest::Vector{Int}, 
    eps::T; 
    μ::T=one(T)
) where {T <: Number}
    fill!(A, zero(T))

    Threads.@threads for i in axes(force_pts, 2)
        S = MMatrix{3,3,T}(undef)
        diffvec = MVector{3,T}(undef)

        for j in axes(quad_pts, 2)
            @inbounds for k in 1:3
                diffvec[k] = force_pts[k,i] - quad_pts[k,j]
            end
            n = nearest[j]
            regularised_stokeslet!(S, diffvec; eps=eps)
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


