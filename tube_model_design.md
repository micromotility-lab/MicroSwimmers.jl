# TubeModel / LineTubeModel — Design Note

## Principle

The model owns "how to place points on this geometry". The discretisation is always generic storage (`NewNearestDiscretisation(N_total, Q_total)`). No tube geometry belongs in the disc type.

`TubeFlagellumNearestDiscretisation` and `LineTubeFlagellumNearestDiscretisation` can be deleted once this is implemented.

## Model wrappers

```julia
struct TubeModel{M <: FlagellumModel, T} <: FlagellumModel
    inner::M
    N_cs::Int   # cross-section force points
    Q_cs::Int   # cross-section quad points
    radius::T
end

struct LineTubeModel{M <: FlagellumModel, T} <: FlagellumModel
    inner::M
    Q_cs::Int   # cross-section quad points only (force pts stay on centreline)
    radius::T
end
```

Neither implements `unit_tangent` — they delegate to `inner` for centreline geometry and override only the disc-glue functor:

```julia
function (m::TubeModel)(disc::NewNearestDiscretisation, t::T) where {T <: Number}
    integrate_tube!(disc.force_pts, disc.velocity, m.inner, m.N_cs, m.radius, t)
    integrate_tube!(disc.quad_pts,                 m.inner, m.Q_cs, m.radius, t)
end

function (m::LineTubeModel)(disc::NewNearestDiscretisation, t::T) where {T <: Number}
    m.inner(disc.force_pts, disc.velocity, t; include_endpoints=false)  # centreline force pts
    integrate_tube!(disc.quad_pts, m.inner, m.Q_cs, m.radius, t)        # tube quad pts
end
```

## Usage

```julia
tube_part      = Part(TubeModel(PlanarFlagellum(...), N_cs, Q_cs, radius), N * N_cs, Q * Q_cs)
line_tube_part = Part(LineTubeModel(PlanarFlagellum(...), Q_cs, radius), N, Q * Q_cs)
```

## implement `integrate_tube!`

At each arclength station: get `t̂` from `unit_tangent`, get `ṫ` from `unit_tangent_and_dt`.
Build `n̂` and `b̂` via **parallel transport** (propagate the frame along the centreline) rather than the Frenet formula, to avoid issues at inflection points where curvature → 0.
Distribute `N_cs` / `Q_cs` points around the circumference in the normal–binormal plane.

## VanedFlagellum

`VanedFlagellum` has the same problem — vane parameters are encoded in the Flagellum struct. It could follow the same wrapper pattern as `VaneModel{M}` in a later cleanup pass.
