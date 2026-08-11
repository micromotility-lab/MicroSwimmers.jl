abstract type BackgroundFlow end

struct NoFlow <: BackgroundFlow end
@inline (::NoFlow)(x, t) = zero(x)

struct LinearFlow{T} <: BackgroundFlow
    G::SMatrix{3,3,T}      # velocity gradient: u∞ = U0 + G*x
    U0::SVector{3,T}
end
@inline (f::LinearFlow)(x, t) = f.U0 + f.G * x