# ============================================================================
#  Tube flagella — a filament with a real cross-sectional radius
#
#  On a bare FlagellumModel the filament radius is carried by `eps`: the force
#  points sit on the centreline and the regularisation blob stands in for the
#  filament's thickness (see the header of defaults.jl). That conflates two
#  different things. `eps` is an isotropic 3D blob centred on the source point,
#  whereas a filament of radius `a` is a surface — regularised transversally
#  only, and carrying azimuthal structure the blob has no way to represent.
#
#  Two wrappers make the radius explicit geometry instead:
#
#    LineTubeFlagellum     force points stay on the centreline; the quadrature
#                          points move onto a ring of radius `a` around it.
#                          The linear system is UNCHANGED in size and in its
#                          collocation points, so this is accuracy for free.
#
#    SurfaceTubeFlagellum  force points and quadrature points both live on the
#                          tube surface. The reference case, against which the
#                          LineTube and the plain blob get validated.
#
#  Both are still regularised-Stokeslet BEM over a nearest-neighbour point
#  cloud. Nothing here is a panel method: a tube is a different distribution of
#  the same kind of points, not a different discretisation of the integral.
# ============================================================================

abstract type TubeFlagellum{FM <: FlagellumModel} <: FlagellumModel end

# Force points on the centreline, quadrature on a ring of radius `radius`.
mutable struct LineTubeFlagellum{FM <: FlagellumModel, T <: Number} <: TubeFlagellum{FM}
    flagellum::FM
    radius::T
    Q_cs::Int          # quadrature points per cross-section
end

LineTubeFlagellum(m::FlagellumModel, a; Q_cs::Int=DEFAULT_Q_CS) =
    LineTubeFlagellum(m, float(a), Q_cs)

# Force points and quadrature points both on the tube surface.
mutable struct SurfaceTubeFlagellum{FM <: FlagellumModel, T <: Number} <: TubeFlagellum{FM}
    flagellum::FM
    radius::T
    N_cs::Int          # force points per cross-section
    Q_cs::Int          # quadrature points per cross-section
end

SurfaceTubeFlagellum(m::FlagellumModel, a; N_cs::Int=DEFAULT_N_CS, Q_cs::Int=DEFAULT_Q_CS) =
    SurfaceTubeFlagellum(m, float(a), N_cs, Q_cs)

# ---------------------------------------------------------------------------
#  Forwarding: a tube is its centreline plus a radius, so every question about
#  the curve goes to the wrapped model. This is what keeps the generic
#  FlagellumModel machinery — the (s,t)-grid method, the plotted centreline,
#  arclength-based sizing — working on a tube.
# ---------------------------------------------------------------------------
arclength(m::TubeFlagellum)    = arclength(m.flagellum)
radius(m::TubeFlagellum)       = m.radius
# Lateral area only: both tubes are open at s=0 and s=1 (see the file header).
surface_area(m::TubeFlagellum) = 2π * m.radius * arclength(m)

unit_tangent(s, t, m::TubeFlagellum)        = unit_tangent(s, t, m.flagellum)
unit_tangent_and_dt(s, t, m::TubeFlagellum) = unit_tangent_and_dt(s, t, m.flagellum)

# ---------------------------------------------------------------------------
#  Station marching
# ---------------------------------------------------------------------------

# March the centreline once, handing each station (i, s, x, τ) to `f`.
#
# This deliberately mirrors integrate_centreline! rather than calling it: a ring
# needs the tangent at each station as well as the position, and asking the
# model for τ a second time afterwards would double the tangent evaluations,
# which are the whole cost of the geometry update. The duplication is guarded by
# a test asserting the two agree bitwise — that test is what stops them drifting.
@inline function march_stations(f::F, m::FlagellumModel, N::Int, t::T;
                                include_endpoints::Bool) where {F, T <: Number}
    s_prev, ds = get_s0_and_ds(T, N, include_endpoints)
    half_L_ds  = T(0.5) * arclength(m) * ds

    τ_prev = unit_tangent(s_prev, t, m)
    x = if include_endpoints
        zero(SVector{3,T})
    else
        τ0 = unit_tangent(zero(T), t, m)          # integrate the first panel from s=0
        (τ0 + τ_prev) * half_L_ds
    end
    f(1, s_prev, x, τ_prev)

    @inbounds for i in 2:N
        s = s_prev + ds
        τ = unit_tangent(s, t, m)
        x = x + (τ_prev + τ) * half_L_ds
        f(i, s, x, τ)
        s_prev, τ_prev = s, τ
    end
    nothing
end

# ---------------------------------------------------------------------------
#  Cross-section frames
# ---------------------------------------------------------------------------

# Any frame orthogonal to τ will do for a *quadrature* ring: the ring is the same
# point set whatever the azimuthal offset, and the equally-weighted periodic
# trapezoid rule over it is invariant to that offset to spectral accuracy. So the
# LineTube can use this per-station Gram-Schmidt construction, branch and all.
#
# It is NOT rotation-minimising — the |τ_x| > 0.9 branch makes it jump — which is
# exactly why SurfaceTubeFlagellum, whose force points move with the frame, uses
# transport_frame below instead.
@inline function ring_frame(τ::SVector{3,T}) where {T <: Number}
    u = abs(τ[1]) > T(0.9) ? SVector{3,T}(0, 1, 0) : SVector{3,T}(1, 0, 0)
    n = normalize(u - dot(τ, u) * τ)
    (n, cross(τ, n))
end

# Parallel transport (Bishop): rotate the previous normal by the minimal rotation
# carrying τ_prev onto τ. Frenet is not an option here — its normal is undefined
# where the curvature vanishes, which on a beating flagellum happens at every
# inflection, i.e. twice per wavelength per beat.
@inline function transport_frame(n::SVector{3,T}, τ_prev::SVector{3,T},
                                 τ::SVector{3,T}) where {T <: Number}
    v  = cross(τ_prev, τ)
    sn = norm(v)
    n_rot = sn > eps(float(real(T))) ? rotation_matrix(v, atan(sn, dot(τ_prev, τ))) * n : n
    # Re-orthogonalise against τ: the rotation is exact in exact arithmetic, but
    # the error accumulates over a few hundred stations if left alone.
    nb = normalize(n_rot - dot(n_rot, τ) * τ)
    (nb, cross(τ, nb))
end

# ---------------------------------------------------------------------------
#  Geometry fill
# ---------------------------------------------------------------------------

# `Q_cs` points evenly spaced round a circle of radius `a` at each station.
function fill_rings!(pts::AbstractVector{SVector{3,T}}, m::FlagellumModel, a, Q_cs::Int,
                     t::T; include_endpoints::Bool) where {T <: Number}
    N  = length(pts) ÷ Q_cs
    aT = T(a)
    march_stations(m, N, t; include_endpoints=include_endpoints) do i, s, x, τ
        n, b = ring_frame(τ)
        base = (i - 1) * Q_cs
        @inbounds for j in 1:Q_cs
            sα, cα = sincos(2 * T(π) * (j - 1) / Q_cs)
            pts[base + j] = x + aT * (cα * n + sα * b)
        end
    end
    pts
end

# A full tube surface: `n_cs` points round each station, on a parallel-transported frame.
#
# Unlike a LineTube's quadrature ring, the azimuthal offset matters here — these points are
# the unknowns, and a frame that jumps between stations would put a seam in them. So the
# normal is carried along the curve by transport_frame rather than rebuilt per station.
function fill_tube!(pts::AbstractVector{<:SVector{3}}, m::FlagellumModel, a, n_cs::Int,
                    t::T; include_endpoints::Bool) where {T <: Number}
    N = length(pts) ÷ n_cs
    carried = Ref((zero(SVector{3,T}), zero(SVector{3,T})))     # (normal, previous tangent)
    march_stations(m, N, t; include_endpoints=include_endpoints) do i, s, x, τ
        n, b = if i == 1
            ring_frame(τ)                                        # seed: any frame will do
        else
            n_prev, τ_prev = carried[]
            transport_frame(n_prev, τ_prev, τ)
        end
        carried[] = (n, τ)
        base = (i - 1) * n_cs
        @inbounds for j in 1:n_cs
            sα, cα = sincos(2 * T(π) * (j - 1) / n_cs)
            pts[base + j] = x + a * (cα * n + sα * b)
        end
    end
    pts
end

# Positions *and* velocities on the tube surface.
#
# There is no local formula for the surface velocity. A point on the surface moves with
# v_c + ω × offset, but ω is the rate of the transported frame, and the transported normal at
# station i depends on the normal at i-1 — so ∂n/∂t obeys a recurrence along the curve, not a
# pointwise identity in τ and ∂τ/∂t. Rather than hand-differentiate that recurrence, run the
# whole march on a Dual-valued t and read the positions off the values and the velocities off
# the partials. Exact, and about twice the cost of positions alone.
#
# Note this needs nothing from the concrete models: march_stations derives its arclength
# station type from `typeof(t)`, so a Dual t promotes s along with it and the existing
# `unit_tangent(s::T, t::T, m)` methods match unchanged.
struct TubeVelocityTag end

function fill_tube!(pts::AbstractVector{SVector{3,T}}, vels::AbstractVector{SVector{3,T}},
                    m::FlagellumModel, a, n_cs::Int, t::T;
                    include_endpoints::Bool) where {T <: Number}
    td  = ForwardDiff.Dual{TubeVelocityTag}(t, one(T))
    dual = Vector{SVector{3,typeof(td)}}(undef, length(pts))
    fill_tube!(dual, m, a, n_cs, td; include_endpoints=include_endpoints)
    @inbounds for i in eachindex(pts)
        p = dual[i]
        pts[i]  = SVector{3,T}(ForwardDiff.value(p[1]),
                               ForwardDiff.value(p[2]),
                               ForwardDiff.value(p[3]))
        vels[i] = SVector{3,T}(ForwardDiff.partials(p[1], 1),
                               ForwardDiff.partials(p[2], 1),
                               ForwardDiff.partials(p[3], 1))
    end
    pts
end

function (m::SurfaceTubeFlagellum)(disc::NearestDiscretisation, t::T) where {T <: Number}
    f_rng = disc.force_part_ranges[1]
    q_rng = disc.quad_part_ranges[1]
    # Both clouds use the transported frame, so a quadrature point and the force point it is
    # nearest to sit at the same azimuth. Seeding them differently would scramble the
    # nearest-neighbour partition even though each cloud on its own looked fine.
    fill_tube!(view(disc.force_pts, f_rng), view(disc.velocity, f_rng),
               m.flagellum, m.radius, m.N_cs, t; include_endpoints=false)
    fill_tube!(view(disc.quad_pts, q_rng), m.flagellum, m.radius, m.Q_cs, t;
               include_endpoints=true)
    nothing
end

function (m::LineTubeFlagellum)(disc::NearestDiscretisation, t::T) where {T <: Number}
    f_rng = disc.force_part_ranges[1]
    q_rng = disc.quad_part_ranges[1]
    # The force half is exactly the untubed model — same nodes, same velocities.
    # That is the point of the LineTube: only the quadrature moves onto the
    # surface, so the matrix keeps both its size and its collocation points.
    m.flagellum(view(disc.force_pts, f_rng), view(disc.velocity, f_rng), t;
                include_endpoints=false)
    fill_rings!(view(disc.quad_pts, q_rng), m.flagellum, m.radius, m.Q_cs, t;
                include_endpoints=true)
    nothing
end

# ---------------------------------------------------------------------------
#  eps vs the tube radius
# ---------------------------------------------------------------------------

check_tube(::Model, N, eps) = nothing
check_tube(::TubeFlagellum, N, eps) = nothing

# What counts as a sensible eps differs between the two tubes, so the checks do not share a
# method. On a LineTube eps must sit well below the ring radius, because the ring is what is
# meant to be setting the transverse length scale. On a SurfaceTube eps is an ordinary
# surface regularisation and should track the point spacing instead — which for a slender
# filament is *smaller* than the radius, so the LineTube rule would fire on every sensible
# SurfaceTube.

# eps and the tube radius are deliberately independent — exploring the (eps, a) plane is the
# entire point of having a tube. But once eps grows past roughly a/2 the blob is wider than
# the ring, the ring geometry stops setting the transverse length scale, and the tube is
# doing nothing the untubed model was not already doing more cheaply.
#
# Measured: at L=10, a=0.1, N=61 the drag plateaus to within 0.1% once eps/a drops below
# about 0.05, and moves by under 0.5% from eps/a = 0.1 to 0.01. So eps ~ a/10 is already in
# the converged regime; a/2 is where it stops being.
#
# The second check is sharper than an accuracy note. A LineTube collocates on the axis, at
# distance exactly `radius` from every point of its own quadrature ring, while the
# neighbouring rings are sqrt(hf^2 + radius^2) away. Once hf drops below the radius those two
# distances are within a few percent of each other, the nearest-neighbour partition stops
# tracking which ring belongs to which force point, and it breaks: measured at L=10, a=0.1,
# translating drag against N runs 16.83, 16.94, 16.94 while hf/a > 0.8, then 16.62, 16.06,
# 19.64 as hf/a falls to 0.28 — non-monotone and diverging. At a=0.02, where hf/a stays above
# 1.3, the same sweep is smooth and monotone throughout.
#
# So refining a LineTube past hf ~ radius makes it worse, not better. That is the opposite of
# what a mesh-refinement study expects, which is why it earns a warning and not a docstring.
function check_tube(m::LineTubeFlagellum, N, eps)
    ε = minimum(eps)
    ε > m.radius / 2 && @warn "Part(LineTubeFlagellum) has eps=$ε against a tube radius of " *
        "$(m.radius). Above eps ~ radius/2 the regularisation blob is wider than the " *
        "cross-section, so the ring geometry no longer sets the transverse length scale " *
        "and the tube behaves like the untubed model. Try eps ~ radius/10." maxlog=1

    hf = arclength(m) / (N + 1)
    hf < m.radius && @warn "Part(LineTubeFlagellum) has force-point spacing hf=$hf below its " *
        "tube radius $(m.radius) (N=$N). The quadrature rings of neighbouring stations are " *
        "then no nearer their own force point than their neighbours', and the solution " *
        "degrades under further refinement instead of converging. Use N < " *
        "$(floor(Int, arclength(m)/m.radius) - 1), or a SurfaceTubeFlagellum." maxlog=1
    nothing
end

# A SurfaceTube has no such refinement limit — its force points are on the surface, so
# refining just makes the surface cloud finer. What it does have is a shape constraint: the
# cells are hf_axial by hf_ring, and the nearest-neighbour patches only mean anything if
# those are within an order of magnitude of each other. A slender filament sized by the
# arclength rule gets very long thin cells unless N_cs is raised to match.
function check_tube(m::SurfaceTubeFlagellum, N, eps)
    axial = arclength(m) / (N + 1)
    ring  = 2π * m.radius / m.N_cs
    r = max(axial, ring) / min(axial, ring)
    r > 10 && @warn "Part(SurfaceTubeFlagellum) has surface cells of aspect ratio " *
        "$(round(r, digits=1)) (axial spacing $(round(axial, sigdigits=3)) against ring " *
        "spacing $(round(ring, sigdigits=3))). Set N_cs ~ " *
        "$(max(round(Int, 2π*m.radius/axial), 3)) to square them up." maxlog=1
    nothing
end
