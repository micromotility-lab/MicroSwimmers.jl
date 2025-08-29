const I3 = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]

skew_symmetric_static(x::T) where{T} = @SMatrix [0.0   -x[3]   x[2];
                                                 x[3]   0.0   -x[1];
                                                -x[2]   x[1]   0.0]


        
function regularised_stokeslet2!(S::StaticMatrix{3,3,T}, R::StaticVector{3,T}; eps::T=1e-6) where {T <: Number}
    rsqr = dot(R,R)
    diag = rsqr + 2eps^2
    denom = 1 / sqrt(rsqr + eps^2)^3

    @inbounds for i in 1:3, j in 1:3
        S[i,j] = diag * (i == j) + R[i] * R[j]
        S[i,j] *= denom
    end
end

function regularised_stokeslet!(S::MMatrix{3,3,Float64}, R::Union{SVector{3,Float64}, MVector{3,Float64}}; eps=1e-6)
    rsqr = sum(R .^ 2)
    factor = 1 / (rsqr + eps^2)^(3//2)
    @inbounds for i in 1:3, j in 1:3
        S[i,j] = (rsqr + 2eps^2) * (i == j) + R[i] * R[j]
        S[i,j] *= factor
    end
end

function resistance_matrix!(
    A::Matrix{T},
    force_pts::Matrix{T},
    quad_pts::Matrix{T},
    nearest::Vector{Int},
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


function swimming_matrix2!(
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

         Kq = skew_symmetric_static(@SVector [quad_pts[1,j], quad_pts[2,j], quad_pts[3,j]])
         @inbounds for p in 1:3, q in 1:3
             A[end-3+p, 3n-3+q] += Kq[p,q]
         end
    end
end


function swimming_matrix!(A, x0, force_pts, quad_pts, nearest, eps; μ=1.)
    fill!(A, 0.0)
    # T = zeros(3,3)
    Threads.@threads for i in axes(force_pts, 2)
        S = MMatrix{3,3,Float64}(undef)
        diffvec = MVector{3,Float64}(undef)
        for j in axes(quad_pts, 2)
            @inbounds @simd for k in 1:3
                diffvec[k] = force_pts[k,i] - quad_pts[k,j]
            end
            n = nearest[j]
            regularised_stokeslet!(S, diffvec; eps=eps)
            @views A[3i-2:3i, 3n-2:3n] .-= S
        end
        @views A[3i-2:3i, end-5:end-3] .= -I3
        @inbounds @simd for k in 1:3
            diffvec[k] = force_pts[k,i] - x0[k]
        end
        @views A[3i-2:3i, end-2:end] .= skew_symmetric_static(diffvec)
    end

    @views A[1:length(force_pts),1:length(force_pts)] ./= -8π*μ
    
    for j in axes(quad_pts, 2)
        n = nearest[j]
        @views A[end-5:end-3, 3n-2:3n] .+= I3
        @views A[end-2:end, 3n-2:3n] .+= skew_symmetric_static(quad_pts[:,j])
    end
end
