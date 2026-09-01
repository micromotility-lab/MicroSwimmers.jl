### Ellipsoids

# Knud Thomsen's approximation. The exact surface area needs incomplete elliptic integrals;
# this is within 1.1% for every axis ratio, which is far tighter than the accuracy of picking
# a point count from a target spacing.
ellipsoid_area(a, b, c; p = 1.6075) = 4π * (((a*b)^p + (a*c)^p + (b*c)^p) / 3)^(1/p)

function fibonacci_ellipsoid(a::T, b::T, c::T, num_points::Int) where {T <: Number}
    points = Vector{SVector{3,T}}(undef, num_points)
    phi = (sqrt(5) + 1) / 2 - 1  # Golden ratio minus 1
    ga = 2π * phi  

    for i in 0:num_points-1
        theta = ga * i  
        z = 1 - (2i + 1) / num_points  # height
        r = sqrt(1 - z^2)  # radius of the circle at height z

        x = a * cos(theta) * r
        y = b * sin(theta) * r
        z = c * z

        points[i+1] = SVector(x, y, z)  
    end
    points
end


# rejection sampled version
function fibonacci_ellipsoid_rejection(a, b, c, num_points)
    # area element ||∂_θ r × ∂_z r|| as a function of (θ, z)
    function area_element(θ, z)
        r = sqrt(1 - z^2)
        sqrt(a^2*b^2*z^2/(1-z^2 + eps()) + c^2*(a^2*sin(θ)^2 + b^2*cos(θ)^2))
    end
    
    # maximum area element for rejection sampling
    M = maximum(area_element(θ, z) 
                for θ in range(0, 2π, 100), z in range(-0.99, 0.99, 100))
    
    points = zeros(3, num_points)
    i = 0
    while i < num_points
        θ = 2π * rand()
        z = 2*rand() - 1
        r = sqrt(1 - z^2)
        if rand() < area_element(θ, z) / M
            i += 1
            points[:, i] = [a*r*cos(θ), b*r*sin(θ), c*z]
        end
    end
    points
end

is_inside_ellipsoid(x, center, radii; orientation=I3, tol=1e-8) = sum((orientation' * (x .- center) ./ radii) .^ 2) <= 1.0 + tol

# Needs NonlinearSolve so I'm not sure whether to include just for this
#
# function ellipsoid_intersection(;x0=[-0.2, 0.05], p=[(0., 0.), (0.2, 0.1), (0., .1), (0.2 ,0.1)])
#     function f!(res, u, p)
#         X, A, Y, B = p
#         res .= [sum(((u .- X) ./ A).^2) - 1.,  sum(((u .- Y) ./ B).^2) - 1.]
#     end
#     prob = NonlinearProblem(f!, x0, p)
#     sol = solve(prob)
#     [sol[1]; 0.; sol[2]]
# end


## Cylinder


function fibonacci_cylinder(a::T, b::T, h::T, num_side::Int, num_cap::Int) where {T <: Number}
    # Golden angle
    phi = (sqrt(5) + 1) / 2 - 1
    ga = 2π * phi

    ### 1. Lateral surface: z ∈ (0, h) ###
    side = Array{T}(undef, 3, num_side)
    for i in 0:num_side - 1
        θ = ga * i
        z = h * (i + 0.5) / num_side          # in (0, h)
        x = a * cos(θ)
        y = b * sin(θ)
        side[:, i+1] = SVector(x, y, z)  # or [x, y, z]
    end

    ### 2. Top cap: z = h ###
    top = Array{T}(undef, 3, num_cap)
    for i in 0:num_cap - 1
        θ = ga * i
        r = sqrt((i + 0.5) / num_cap)
        x = a * r * cos(θ)
        y = b * r * sin(θ)
        z = h
        top[:, i+1] = [x, y, z]
    end

    ### 3. Bottom cap: z = 0 ###
    bottom = Array{T}(undef, 3, num_cap)
    for i in 0:num_cap - 1
        θ = ga * i
        r = sqrt((i + 0.5) / num_cap)
        x = a * r * cos(θ)
        y = b * r * sin(θ)
        z = zero(T)
        bottom[:, i+1] = [x, y, z]
    end

    return side, top, bottom
end


is_inside_cylinder(x, center, dims) = x[3] > center[3] - dims[3] && sum(((x[1:2] .- center[1:2]) ./ dims[1:2]).^2) < 1 

## Raymarching alternative

# --- Fibonacci direction on the unit sphere ---
function fib_dir(i, N)
    ϕ = (√5 + 1)/2 - 1
    z = 1 - (2i + 1)/N
    r = sqrt(max(0.0, 1 - z^2))
    θ = 2π*ϕ*i
    SVector(cos(θ)*r, sin(θ)*r, z)
end

# --- bisection on a bracketed sign change of g(t) ---
function bisect(g, lo, hi; tol=1e-12, max_iter=60)
    glo = g(lo)
    for _ in 1:max_iter
        mid = (lo + hi)/2
        gmid = g(mid)
        (hi - lo) < tol && return mid
        sign(gmid) == sign(glo) ? (lo = mid; glo = gmid) : (hi = mid)
    end
    (lo + hi)/2
end

# Ray-march an implicit surface f(x) = 0 from `seed` along N Fibonacci directions, pushing
# every sign change onto `pts` and its surface quadrature weight onto `wts`.
#
# The weight is the Jacobian of the (solid angle -> surface area) map:
#
#     dA = t^2 * ΔΩ / |n·d|
#
# with ΔΩ = 4π/N the solid angle carried by each ray, t the hit distance and n the unit
# surface normal. `sum(wts)` is therefore the surface area — see `total_quad_weight`.
#
# Two caveats, both real for the grooved bodies in implicit_body.jl:
#
#  * The formula is exact only when the surface is star-shaped about `seed`, i.e. each ray
#    crosses exactly once. A ray that crosses 2-3 times (as happens on
#    ImplicitGroovedEllipsoid, which is not star-shaped about `groove_center`) gives every
#    crossing the full ΔΩ of that ray, over-counting the area.
#  * Grazing hits (|n·d| -> 0) send the weight to infinity. We warn rather than clamp: a
#    clamp would silently bias the quadrature, whereas a warning points at the seed choice,
#    which is the actual fix.
function raymarch_cloud!(pts, wts, f, R, N; seed=zero(SVector{3,Float64}),
                        num_steps=200, max_iter=60, tol=1e-12, grazing_tol=1e-3)
    ΔΩ = 4π/N
    grazing = 0
    for i in 0:N-1
        d  = fib_dir(i, N)
        ts = range(0, R; length=num_steps)
        for j in 1:num_steps-1
            p0, p1 = seed + ts[j]*d, seed + ts[j+1]*d
            if sign(f(p0)) != sign(f(p1))
                t = bisect(τ -> f(seed + τ*d), ts[j], ts[j+1]; tol, max_iter)
                p = seed + t*d
                n = normalize(ForwardDiff.gradient(f, p))
                cosθ = abs(dot(n, d))
                cosθ < grazing_tol && (grazing += 1)
                push!(pts, p)
                push!(wts, t^2 * ΔΩ / cosθ)
            end
        end
    end
    grazing > 0 && @warn "raymarch_cloud!: $grazing of $(length(pts)) hits are near-grazing \
        (|n·d| < $grazing_tol); their quadrature weights are unreliable. Consider moving `seed`."
    pts, wts
end

# points-only convenience method: march the surface and discard the weights
function raymarch_cloud!(pts, f, R, N; kwargs...)
    # grazing_tol=0 suppresses the weight warning: no weights are kept, so it is noise here
    raymarch_cloud!(pts, Float64[], f, R, N; grazing_tol=0.0, kwargs...)
    pts
end

# """
# Repulsion relaxation of a point cloud constrained to an implicit surface {f = 0}.

# Takes an arbitrary (e.g. ray-marched, non-uniformly spaced) cloud and relaxes it
# towards a quasi-uniform distribution by tangential Gaussian repulsion plus Newton
# re-projection onto the level set.

# """

# # ---------------------------------------------------------------------------
# # surface primitives
# # ---------------------------------------------------------------------------

# surface_normal(f, p::SVector{3,Float64}) = normalize(ForwardDiff.gradient(f, p))

# """
#     project_to_surface(f, p; tol=1e-12, max_iter=8)

# Pull `p` back onto {f = 0} by Newton iteration along ∇f.

# This is *not* a closest-point projection — it walks along the gradient direction,
# which is only the shortest route in the limit. That is fine here: the tangential
# step that knocked the point off the surface is small by construction.

# Converges quadratically. For an implicit function whose |∇f| varies a lot over
# the surface (e.g. an elongated ellipsoid) the step length varies with it, but
# convergence is unaffected.
# """
# function project_to_surface(f, p::SVector{3,Float64}; tol=1e-12, max_iter=8)
#     for _ in 1:max_iter
#         v = f(p)
#         abs(v) < tol && return p
#         g = ForwardDiff.gradient(f, p)
#         p -= (v / dot(g, g)) * g
#     end
#     return p
# end

# # ---------------------------------------------------------------------------
# # diagnostics
# # ---------------------------------------------------------------------------

# """
#     spacing_stats(pts)

# Nearest-neighbour spacing statistics. `cv` (coefficient of variation) is the
# headline number: it measures how uniform the point set is, independent of scale.

# Rough guide:
#   cv ≳ 0.5    strongly non-uniform (typical raw ray-marched cloud on a 5:1:1 body)
#   cv ≈ 0.15   partially relaxed
#   cv ≈ 0.05   well relaxed; local packing is near-hexagonal
# """
# function spacing_stats(pts::Vector{SVector{3,Float64}})
#     tree = KDTree(pts)
#     _, dists = knn(tree, pts, 2, true)      # dd[1] is the point itself
#     d = [dd[2] for dd in dists]
#     return (min = minimum(d), mean = mean(d), max = maximum(d), cv = std(d) / mean(d))
# end

# # ---------------------------------------------------------------------------
# # the relaxation itself
# # ---------------------------------------------------------------------------

# """
#     relax_cloud!(pts, f; area, kwargs...)

# Relax `pts` in place towards a quasi-uniform distribution on {f = 0}.

# `area` is the total surface area, used only to set the repulsion length scale
# σ ≈ √(area/N). A rough value is fine — `sum(wts)` from a ray march works well,
# because the local weight errors cancel in the sum even when individual weights
# are poor.

# Keyword arguments
#   iters         maximum number of relaxation sweeps
#   σ_mult        repulsion range as a multiple of target spacing √(area/N)
#   η             step size (displacement = η·σ·F)
#   max_step      hard cap on displacement per sweep, in units of σ
#   cutoff        neighbour search radius, in units of σ (3σ captures >99% of the Gaussian)
#   rebuild_every rebuild the KD-tree every this many sweeps
#   tol           stop when the largest displacement falls below tol·σ

# Uses a Jacobi update (all forces computed from the old positions, then all points
# moved at once). Gauss–Seidel converges in fewer sweeps but is order-dependent and
# not thread-safe.
# """
# function relax_cloud!(pts::Vector{SVector{3,Float64}}, f;
#                       area::Real,
#                       iters         = 300,
#                       σ_mult        = 1.2,
#                       η             = 0.15,
#                       max_step      = 0.4,
#                       cutoff        = 3.0,
#                       rebuild_every = 4,
#                       tol           = 1e-4,
#                       verbose       = true)

#     N  = length(pts)
#     h  = sqrt(area / N)              # target mean spacing
#     σ  = σ_mult * h
#     rc = cutoff * σ
#     inv2σ² = 1 / (2σ^2)

#     disp = Vector{SVector{3,Float64}}(undef, N)
#     tree = KDTree(pts)

#     if verbose
#         s = spacing_stats(pts)
#         @info "relax_cloud! start" N σ=round(σ, sigdigits=4) cv=round(s.cv, digits=4)
#     end

#     for it in 1:iters
#         (it - 1) % rebuild_every == 0 && (tree = KDTree(pts))
#         nbrs = inrange(tree, pts, rc)

#         Threads.@threads for i in 1:N
#             p = pts[i]
#             F = zero(SVector{3,Float64})
#             wsum = 0.0
#             @inbounds for j in nbrs[i]
#                 j == i && continue
#                 d  = p - pts[j]
#                 r2 = dot(d, d)
#                 r2 < 1e-30 && continue
#                 # Gaussian repulsion: bounded as r → 0, so a coincident pair
#                 # cannot blow the explicit step up.
#                 k = exp(-r2 * inv2σ²)
#                 F += d * (k / sqrt(r2)); wsum += k
#             end
#             wsum > 0 && (F /= wsum)

#             # keep only the tangential component — the point slides, never lifts
#             n  = surface_normal(f, p)
#             F -= dot(F, n) * n

#             δ  = (η * σ) * F
#             nδ = norm(δ)
#             nδ > max_step * σ && (δ *= (max_step * σ) / nδ)
#             disp[i] = δ
#         end

#         moved = 0.0
#         @inbounds for i in 1:N
#             pts[i] = project_to_surface(f, pts[i] + disp[i])
#             moved  = max(moved, norm(disp[i]))
#         end

#         if verbose && (it == 1 || it % 25 == 0)
#             s = spacing_stats(pts)
#             @info "relax_cloud!" iter=it max_move_over_σ=round(moved/σ, digits=4) cv=round(s.cv, digits=4)
#         end

#         if moved < tol * σ
#             verbose && @info "relax_cloud! converged" iter=it
#             break
#         end
#     end

#     if verbose
#         s = spacing_stats(pts)
#         @info "relax_cloud! done" cv=round(s.cv, digits=4) min_spacing=round(s.min, sigdigits=4)
#     end
#     return pts
# end

