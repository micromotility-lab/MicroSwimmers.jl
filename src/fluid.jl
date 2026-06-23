# Calculate the velocity at an arbitrary point in the fluid

struct FluidVelocity{T <: Number}
    quad_pts::Vector{SVector{3,T}}
    nearest::Vector{Int}          
    eps::T
    mu::T
    force_vals::Vector{T}   # change this to AbstractVector to accept a view      
    A::Matrix{T}                  # 3×N
    wall::Bool
end

function FluidVelocity(prob::InstantaneousProblem)
    check_solved!(prob)
    N = length(prob.disc.force_pts)
    A = zeros(3, N)
    FluidVelocity(
        prob.disc.quad_pts,
        prob.disc.nearest,
        prob.eps,
        prob.mu,
        prob.force_vals[1:N], # leave out U and Ω for SwimmingProblem
        A,
        prob.wall
    )
end

function (fv::FluidVelocity)(x)
    resistance_matrix!(fv.A, x, fv.quad_pts, fv.nearest, fv.eps; μ=fv.mu, wall=fv.wall)
    SVector{3}(fv.A * fv.force_vals)     
end

struct PlanarVelocityField{T <: Number}
    plane::Symbol   # :xy, :xz or :yz
    a_range::AbstractVector{T}
    b_range::AbstractVector{T}
    c::T
    points::Vector{SVector{3,T}}
    velocities::Vector{SVector{3,T}}
end

na(vf::PlanarVelocityField) = length(vf.a_range)
nb(vf::PlanarVelocityField) = length(vf.b_range)

function points3(a_range::AbstractVector{T}, b_range::AbstractVector{T}, c::T, plane::Symbol) where {T <: Number}
    if plane === :xy
        [SVector{3,T}(x, y, c) for y in b_range for x in a_range] # order works well with reshape
    elseif plane === :xz
        [SVector{3,T}(x, c, z) for z in b_range for x in a_range]
    elseif plane === :yz
        [SVector{3,T}(c, y, z) for z in b_range for y in a_range]
    else
        throw(ArgumentError("plane must be :xy, :xz, or :yz"))
    end
end

function PlanarVelocityField(prob::InstantaneousProblem, a_range::AbstractVector{T}, b_range::AbstractVector{T}; c::T=0.0, plane::Symbol=:xy) where {T <: Number}
    points = points3(a_range, b_range, c, plane)
    fv = FluidVelocity(prob)
    PlanarVelocityField(plane, a_range, b_range, c, points, fv.(points))
end

function PlanarVelocityField(fv::FluidVelocity, a_range::AbstractVector{T}, b_range::AbstractVector{T}; c::T=0.0, plane=:xy) where {T <: Number}
    points = points3(a_range, b_range, c, plane)
    PlanarVelocityField(plane, a_range, b_range, c, points, fv.(points))
end


function TimeAveragedPlanarVelocityField(prob::InstantaneousProblem, 
    a_range::AbstractVector{T}, 
    b_range::AbstractVector{T}, 
    pre_transform!::Function=update_boundary!; 
    c::T=0.0, plane::Symbol=:xy, period=1.0, num_t=30
) where {T <: Number}
    new_vf = prob -> PlanarVelocityField(prob, a_range, b_range; c=c, plane=plane)
    vfs = time_collect!(prob, pre_transform!, new_vf, period, num_t; endpoint=false)
    PlanarVelocityField(plane, a_range, b_range, c, vfs[1].points, mean(vf.velocities for vf in vfs))
end

function TimeAveragedVelocityField(prob::InstantaneousProblem, pts::Matrix{T}, pre_transform::Function=update_boundary!; period=1.0, num_t=30) where {T <: Number}
    new_vf = prob -> VelocityField(prob, pts)
    vfs = time_collect!(prob, pre_transform, new_vf, period, num_t; endpoint=false)
    VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs))
end

function velocity_flux(u, z_bot, z_top, y_min, y_max; x=0., N=20)
    # Gauss–Legendre nodes and weights on [-1, 1]
    ys_raw, wys = gausslegendre(N)
    ss_raw, wss = gausslegendre(N)

    # Affine transform to [y_min, y_max] and [0, 1]
    ys = 0.5*(y_max - y_min) * (ys_raw .+ 1) .+ y_min
    wys .= 0.5*(y_max - y_min) * wys

    ss = 0.5 * (ss_raw .+ 1)  # [0,1]
    wss .= 0.5 * wss

    z(y, s) = z_bot(y)*(1 - s) + s*z_top(y)
    sum(w1*w2*u([x, yi, z(yi, si)])[1] for (yi, w1) in zip(ys, wys), (si, w2) in zip(ss, wss))
end

function velocity_flux_polar(u, x, y0, z0, R; Nr=20, Nθ=20)
    rs_raw, wrs = gausslegendre(Nr)
    θs_raw, wθs = gausslegendre(Nθ)

    # Affine transforms
    rs = 0.5 * R * (rs_raw .+ 1)  # r ∈ [0, R]
    wrs .= 0.5 * R * wrs          # Jacobian for r

    θs = π * (θs_raw .+ 1)        # θ ∈ [0, 2π]
    wθs .= π * wθs                # Jacobian for θ

    total_flux = 0.0
    for (r, wr) in zip(rs, wrs), (θ, wθ) in zip(θs, wθs)
        y = y0 + r * cos(θ)
        z = z0 + r * sin(θ)
        vel = u([x, y, z])
        total_flux += vel[1] * r * wr * wθ  # extra r from polar area element
    end

    total_flux
end






