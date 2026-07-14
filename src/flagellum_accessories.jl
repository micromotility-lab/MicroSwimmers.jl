abstract type Accessory end

# A vane: a sheet extruded in a particular direction from the parent centreline over span [s_start,s_end]
# to height H.
mutable struct Vane{T} <: Accessory
    direction::SVector{3,T}
    s_start::T # s ∈ [0,1]
    s_end::T
    H::T  # H in μm
end

# nearest index to a target fraction of total length
function nearest_index(s, N, include_endpoints)
    s0, ds = get_s0_and_ds(typeof(s), N, include_endpoints)
    clamp(round(Int, (s - s0) / ds) + 1, 1, N)
end
# total vane points for a given parametrisation
function vane_npoints(v::Vane, N, L, include_endpoints)
    start  = nearest_index(v.s_start, N, include_endpoints)
    stop   = nearest_index(v.s_end,   N, include_endpoints)
    height = Nh(v, N, L)
    height * (stop - start + 1)
end

Nh(v::Vane, N_flagellum, L) = floor(Int, (N_flagellum * v.H) ÷ L)
Nv(v::Vane, N_flagellum) = floor(Int, N_flagellum * (v.s_end - v.s_start))
Nstart(v::Vane, N_flagellum) = floor(Int, N_flagellum*v.s_start)


mutable struct PlanarVanedFlagellum{FM <: FlagellumModel} <: FlagellumModel
    flagellum::FM
    vane::Vane
end

PlanarVanedFlagellum(m::FlagellumModel, s_start, s_end, H) = PlanarVanedFlagellum(m, Vane(-ez, Float64(s_start), Float64(s_end), Float64(H)))

# assumes flagellum points are already stored in the first 1:N_flagellum elements of points
function (v::Vane)(points::Vector{SVector{3,T}}, N_flagellum, L, include_endpoints) where {T <: Number}
    start  = nearest_index(v.s_start, N_flagellum, include_endpoints)
    stop   = nearest_index(v.s_end,   N_flagellum, include_endpoints)
    N_v    = stop - start + 1
    height = Nh(v, N_flagellum, L)
    for i in 1:height
        cstart = N_flagellum + (i-1)*N_v + 1
        cend   = N_flagellum + i*N_v
        z = -i*L / N_flagellum
        for (dst, src) in zip(cstart:cend, start:stop)
            p = points[src]
            points[dst] = SVector{3,T}(p[1], p[2], z)
        end
    end
end

# assumes flagellum points are already stored in the first 1:N_flagellum elements of points
function (v::Vane)(points::Vector{SVector{3,T}}, velocities::AbstractVector{SVector{3,T}}, N_flagellum, L, include_endpoints) where {T <: Number}
    start  = nearest_index(v.s_start, N_flagellum, include_endpoints)
    stop   = nearest_index(v.s_end,   N_flagellum, include_endpoints)
    N_v    = stop - start + 1
    height = Nh(v, N_flagellum, L)
    for i in 1:height
        cstart = N_flagellum + (i-1)*N_v + 1
        cend   = N_flagellum + i*N_v
        z = -i*L / N_flagellum
        for (dst, src) in zip(cstart:cend, start:stop)
            p = points[src]
            points[dst] = SVector{3,T}(p[1], p[2], z)
            v = velocities[src]
            velocities[dst] = SVector{3,T}(v[1], v[2], zero(T))
        end
    end
end

# function (v::Vane)(points::AbstractVector{SVector{3,T}}, velocities::AbstractVector{SVector{3,T}}, N_flagellum, L) where {T <: Number}
#     height = Nh(v, N_flagellum, L)
#     N_v = Nv(v, N_flagellum)
#     start = Nstart(v, N_flagellum)
#     # @info "" height N_v start
#     for i in 1:height
#         cstart = N_flagellum + (i-1)*N_v + 1
#         cend   = N_flagellum + i*N_v

#         z = -i*L / N_flagellum
#         for (dst, src) in zip(cstart:cend, start:start+N_v-1)
#             p = points[src]
#             points[dst]     = SVector{3,T}(p[1], p[2], z)
#             v = velocities[src]
#             velocities[dst] = SVector{3,T}(v[1], v[2], zero(T))
#         end
#     end
# end



# function (m::PlanarVanedFlagellum)(disc, t)
#     nb = length(disc.quad_pts) - sum(nquad, m.accessories)   # base count, recovered
#     m.model(subview(disc, 1:nb), t)                          # reuse base's own fill
#     off = nb
#     for v in m.accessories
#         rng = (off+1):(off+nquad(v))
#         fill_vane!(disc, rng, m, v, t, ex)
#         off += nquad(v)
#     end
#     disc
# end

"""Flagellum with a vane (only extends in the z direction currently)"""
function (m::PlanarVanedFlagellum)(
    points::Vector{SVector{3,T}},
    N_f::Int,
    t::T
) where {T <: Number}
    m.flagellum(points, t)  # Fill flagellum points and velocities
    m.vane(points, N_f, m.flagellum.L)
end

function (m::PlanarVanedFlagellum)(
    points::Vector{SVector{3,T}}, velocities::Vector{SVector{3,T}}, 
    N_f::Int, 
    t::T;
    include_endpoints=false
) where {T <: Number}
    f_points = @view points[1:N_f]
    f_velocities = @view velocities[1:N_f]
    m.flagellum(f_points, f_velocities, t, include_endpoints=include_endpoints) 
    m.vane(f_points, f_velocities, N_f, m.flagellum.L)
end

function (m::PlanarVanedFlagellum)(disc::NearestDiscretisation, t::T) where {T <: Number}
    # force pts
    f_rng = disc.force_part_ranges[1]
    m.flagellum(view(disc.force_pts, f_rng), view(disc.velocity, f_rng), t)
    m.vane(disc.force_pts, disc.velocity, size(f_rng, 1), m.flagellum.L, false)

    # quad pts
    q_rng = disc.quad_part_ranges[1]
    m.flagellum(view(disc.quad_pts, q_rng), t, include_endpoints=false)
    m.vane(disc.quad_pts, size(q_rng, 1), m.flagellum.L, true)
end
