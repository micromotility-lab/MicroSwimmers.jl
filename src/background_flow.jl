# An imposed flow u∞(x) that the swimmer sits in. The total fluid velocity is u∞ plus the
# boundary's disturbance, so no-slip at a collocation point reads
#     Σ S f − U − Ω×r = u_def − u∞
# i.e. a background flow only moves the right-hand side. The matrix, its factorisation and the
# hybrid GMRES preconditioner never see it.

abstract type BackgroundFlow end

struct NoFlow <: BackgroundFlow end

# Mutable for the same reason the models are: parameters get dragged by sliders and stepped by
# optimisers. Concrete field types keep evaluation a plain load.
mutable struct LinearFlow{T} <: BackgroundFlow
    G::SMatrix{3,3,T,9}    # velocity gradient, G[i,j] = ∂u_i/∂x_j
    U0::SVector{3,T}       # velocity at the lab origin
end

function LinearFlow(G::AbstractMatrix, U0::AbstractVector=zeros(3))
    T = promote_type(eltype(G), eltype(U0), Float64)
    LinearFlow{T}(SMatrix{3,3,T,9}(G), SVector{3,T}(U0))
end

LinearFlow(; G=zeros(3,3), U0=zeros(3)) = LinearFlow(G, U0)

UniformFlow(U0::AbstractVector) = LinearFlow(zeros(3,3), U0)

# u = γ x_grad e_dir; the default is the textbook u = γ y e_x
function ShearFlow(γ::Number; dir::Int=1, grad::Int=2)
    G = zeros(3,3)
    G[dir, grad] = γ
    LinearFlow(G)
end

@inline (f::LinearFlow)(x) = f.U0 + f.G * x

# The NoFlow method is a no-op rather than subtracting zeros, so a quiescent solve stays
# bit-identical to the code before background flows existed.
subtract_flow!(v, x, ::NoFlow) = v
function subtract_flow!(v, x, f::BackgroundFlow)
    @inbounds for i in eachindex(v, x)
        v[i] -= f(x[i])
    end
    v
end

add_flow(u, x, ::NoFlow) = u
add_flow(u, x, f::BackgroundFlow) = u + f(x)

check_flow(::NoFlow, wall::Bool) = nothing

function check_flow(f::LinearFlow, wall::Bool)
    abs(tr(f.G)) > sqrt(eps(eltype(f.G))) * max(1, norm(f.G)) &&
        @warn "LinearFlow has tr(G) = $(tr(f.G)) ≠ 0, so it is not incompressible."
    # The Blakelet already enforces no-slip on z = 0, so u∞ must vanish there too:
    # U0 + G*(x, y, 0) = 0 for every x, y.
    wall && (any(!iszero, f.U0) || any(!iszero, f.G[:, 1:2])) &&
        @warn "With wall=true the background flow must vanish on z = 0 (U0 = 0 and only " *
              "G[:,3] nonzero), or it contradicts the wall's no-slip condition."
    nothing
end
