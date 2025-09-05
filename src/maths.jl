const I3 = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]

skew_symmetric_static(x::T) where{T} = @SMatrix [0.0   -x[3]   x[2];
                                                 x[3]   0.0   -x[1];
                                                -x[2]   x[1]   0.0]

smooth_max(a, b, k) = log(exp(k*a) + exp(k*b)) / k

# Cumulative trapezoid integral: returns I[i] = ∫_{x[1]}^{x[i]} y(s) ds
function cumulative_trapz!(
    out::AbstractVector{T},
    x::AbstractVector{T},
    y::AbstractVector{T}
) where {T}
    @assert length(out) == length(x) == length(y)
    @inbounds begin
        out[1] = zero(T)
        for i in eachindex(x)[2:end]
            dx = x[i] - x[i-1]
            out[i] = out[i-1] + (y[i-1] + y[i]) * (dx/2)
        end
    end
    out
end