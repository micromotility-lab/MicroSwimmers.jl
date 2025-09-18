const I3 = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 1.0]

function rotation_matrix(axis::Vector{T}, angle::T) where T
    axis = normalize(axis)  # Ensure it's a unit vector
    x, y, z = axis
    c = cos(angle)
    s = sin(angle)
    C = 1 - c

    @SMatrix [
        c + x^2*C     x*y*C - z*s   x*z*C + y*s;
        y*x*C + z*s   c + y^2*C     y*z*C - x*s;
        z*x*C - y*s   z*y*C + x*s   c + z^2*C
    ]
end

skew_symmetric_static(x::T) where{T} = @SMatrix [0.0   -x[3]   x[2];
                                                 x[3]   0.0   -x[1];
                                                -x[2]   x[1]   0.0]

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