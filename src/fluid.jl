# Calculate the velocity at an arbitrary point in the fluid

struct FluidVelocity{T <: Number}
    quad_pts::Matrix{T}
    nearest::Vector{Int}          
    eps::T
    mu::T
    force_vals::Vector{T}         
    A::Matrix{T}                  # 3×N
end

function FluidVelocity(prob::InstantaneousProblem)
    check_solved!(prob)
    N = length(prob.points.force_pts)
    A = zeros(3, N)
    FluidVelocity(
        prob.points.quad_pts,
        prob.points.nearest,
        prob.eps,
        prob.mu,
        prob.force_vals[1:N], # leave out U and Ω for SwimmingProblem
        A
    )
end

function (fv::FluidVelocity)(x)
    resistance_matrix!(fv.A, x, fv.quad_pts, fv.nearest, fv.eps; μ=fv.mu)
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
### Old stuff


# struct VelocityFunction{T <: Number}
#     prob::InstantaneousProblem
#     A::Matrix{T}
# end

# function VelocityFunction(prob::InstantaneousProblem)
#     check_solved!(prob)    
#     N = length(prob.points.force_pts)
#     A = zeros(3, N)
#     VelocityFunction(prob, A)
# end

# function (vfn::VelocityFunction)(x)
#     resistance_matrix!(
#         vfn.A, x, 
#         vfn.prob.points.quad_pts, 
#         vfn.prob.points.nearest, 
#         vfn.prob.eps;
#         μ=vfn.prob.mu,
#     )

#     vfn.A * vfn.prob.force_vals     
# end

# # Calculate a velocity field on a grid of points

# struct VelocityField{T <: Number}
#     points::Vector{SVector{3, T}}
#     velocities::Vector{SVector{3,T}}
# end

# get_velocities(vf::VelocityField) = vf.velocities

# function VelocityField(prob::InstantaneousProblem, points::Matrix{T}) where {T <: Number}
#     check_solved!(prob)
#     static_points = [SVector{3}(pt) for pt in eachcol(points)]
#     u = VelocityFunction(prob)
#     VelocityField(static_points, u.(static_points))
# end


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

    return total_flux
end


# struct AverageVelocityFunction
#     vfns::Vector{VelocityFunction}
# end

# function AverageVelocityFunction(prob::ResistanceProblem; T=1.0, num_t=20)
#     check_solved!(prob)
#     vfns = VelocityFunction[]
#     for t in range(0,T,num_t)[1:end-1]
#         update_boundary!(prob, t)
#         solve_problem!(prob)
#         push!(vfns, VelocityFunction(prob))
#     end
#     AverageVelocityFunction(vfns)
# end

# (ave_vfn::AverageVelocityFunction)(x) = mean(vfn(x) for vfn in ave_vfn.vfns)




# function get_velocity_function(prob::SwimmingProblem; body_frame=false, xy_plane=false)
#     check_solved!(prob)    
#     N = length(prob.points.force_pts)
#     A = zeros(3, N)
#     U = get_U(prob)
#     Ω = get_Ω(prob)
    
#     function vel(x)
#         if xy_plane
#             resistance_matrix!(
#                 A, [x[1], x[2], 0.],
#                 prob.points.quad_pts,
#                 prob.points.nearest,
#                 prob.eps;
#                 μ=prob.mu
#             )
#         else
#             resistance_matrix!(
#                 A, x, 
#                 prob.points.quad_pts, 
#                 prob.points.nearest, 
#                 prob.eps;
#                 μ=prob.mu,
#             )
#         end

#         u = A * prob.force_vals[1:N]     
#         if body_frame
#             u -=  U + cross(Ω, [x[1], x[2], 0.] - prob.microswimmer.points.location)
#         end
        
#         if xy_plane
#             return SVector{2}(u[1:2])
#         else
#             return SVector{3}(u)
#         end
#     end
# end

# function get_velocity_function(prob::Problem)
#     solve_problem!(prob)
#     N = prob.points.N
#     A = zeros(3, 3N)
    
#     function u(x; exclude_fn=nothing)

#         if !isnothing(exclude_fn)
#              if exclude_fn(x) return NaN end
#         end
        
#         resistance_matrix!(
#                 A, x, 
#                 prob.points.quad_pts, 
#                 prob.points.nearest, 
#                 prob.eps;
#                 μ=prob.mu,
#         )
#         SVector{3}(A * prob.force_vals[1:3N])
#     end
# end




## xy plane, MAYBE UNNECESSARY, A LOT OF REDUNDANCY HERE 
# struct VelocityField{T <: Number}
#     points::Union{Vector{SVector{2,T}}, Vector{SVector{3, T}}}
#     velocities::Union{Vector{SVector{2,T}}, Vector{SVector{3,T}}}
# end

# function VelocityField(prob::SwimmingProblem{T}, points::Matrix{T}; body_frame=false) where {T <: Number}
#     check_solved!(prob)
#     planar = size(points, 1) == 2
#     static_points = planar ? [SVector{2}(pt) for pt in eachcol(points)] : [SVector{3}(pt) for pt in eachcol(points)]
#     u = get_velocity_function(prob, body_frame=body_frame, xy_plane=planar)
#     VelocityField(static_points, u.(eachcol(points)))
# end

# function VelocityField(prob::ResistanceProblem{T}, points::Matrix{T}; exclude_fn=nothing) where {T <: Number}
#     static_points = [SVector{3}(pt) for pt in eachcol(points)]
#     u = get_velocity_function(prob)
#     VelocityField(static_points, u.(eachcol(points); exclude_fn=exclude_fn))
# end

# function VelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y_points::AbstractVector{T}, z::T) where {T <: Number}
#     pts = reduce(hcat, vec([[x,y,z] for x in x_points, y in y_points]))
#     VelocityField(prob, pts)
# end

# function VelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y::T, z_points::AbstractVector{T}) where {T <: Number}
#     pts = reduce(hcat, vec([[x,y,z] for x in x_points, z in z_points]))
#     VelocityField(prob, pts)
# end



# function TimeAveragedVelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, pts::Matrix{T}; period=1.0, num_ts=30) where {T <: Number}
#     vfs = VelocityField[]
#     for t in range(0, period, num_ts)[1:end-1]
#         update_boundary!(prob, t)
#         solve_problem!(prob)
#         push!(vfs, VelocityField(prob, pts))
#     end
#     VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs) )
# end

# function TimeAveragedVelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y_points::AbstractVector{T}, z::T, ; period=1.0, num_ts=30) where {T <: Number}
#     vfs = VelocityField[]
#     for t in range(0, period, num_ts)[1:end-1]
#         update_boundary!(prob, t)
#         solve_problem!(prob)
#         push!(vfs, VelocityField(prob, x_points, y_points, z))
#     end
    
#     VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs))
# end

# function TimeAveragedVelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y::T,  z_points::AbstractVector{T}; period=1.0, num_ts=30) where {T <: Number}
#     vfs = VelocityField[]
#     for t in range(0, period, num_ts)[1:end-1]
#         update_boundary!(prob, t)
#         solve_problem!(prob)
#         push!(vfs, VelocityField(prob, x_points, y, z_points))
#     end
#     VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs))
# end

# using StaticArrays, LinearAlgebra

# function rotation_align_to_x(swim_dir) 
#     v = normalize(swim_dir)
#     c = clamp(dot(v, ex), -1.0, 1.0)

#     # # If already along +x
#     # if isapprox(c, one(T); atol=1e-8)
#     #     return I(3)
#     # end

#     # # If along -x, rotate π around any axis ⟂ x (e.g. ẑ)
#     # if isapprox(c, -one(T); atol=1e-8)
#     #     axis = @SVector [zero(T), zero(T), one(T)]
#     #     return rotation_matrix(axis, π)
#     # end

#     axis = normalize(cross(v, ex))          # rotation axis
#     angle = acos(c)                         # rotation angle
#     return rotation_matrix(axis, angle)
# end


# function TimeAveragedDisturbanceField(prob::SwimmingTrajectoryProblem, x_points::AbstractVector{T}, y_points::AbstractVector{T}, z::T; num_ts=30) where {T <: Number}
#     check_solved!(prob)
#     vfs = VelocityField[]
#     traj = prob.traj
#     U, S, V = svd(reduce(hcat, traj.x))
#     displacement = traj.x[end] - traj.x[1]
#     swimming_direction = dot(displacement, U[:,1]) > 0.0 ? U[:,1] : -U[:,1]
#     swimming_angle = atan(swimming_direction[2], swimming_direction[1]) 
#     # R = rotation_matrix([0., 0., 1.], π + swimming_angle)
#     R = rotation_align_to_x(swimming_direction)
#     sprob = prob.swimming_problem

#     for t in eachindex(traj.t)
#         # move_boundary!(sprob, traj.x[t], traj.b1[t], traj.b2[t], traj.t[t])
#         # centroid = mean(eachcol(sprob.points.quad_pts))
#         # R_body = hcat(traj.b1[t], traj.b2[t], cross(traj.b1[t], traj.b2[t]))
#         # move_boundary!(sprob, R*(traj.x[t] - centroid), R*R_body, traj.t[t])
#         move_boundary!(sprob, traj.x[1], traj.b1[t], traj.b2[t], traj.t[t])
#         solve_problem!(sprob)
#         push!(vfs, VelocityField(sprob, x_points, y_points, z))
#     end
#     VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs))
# end

# function TimeAveragedDisturbanceField(prob::SwimmingTrajectoryProblem, x_points::AbstractVector{T}, y::T, z_points::AbstractVector{T}; num_ts=30) where {T <: Number}
#     check_solved!(prob)
#     vfs = VelocityField[]
#     traj = prob.traj
#     U, S, V = svd(reduce(hcat, traj.x))
#     displacement = traj.x[end] - traj.x[1]
#     swimming_direction = dot(displacement, U[:,1]) > 0.0 ? U[:,1] : -U[:,1]
#     swimming_angle = atan(swimming_direction[2], swimming_direction[1]) 
#     # R = rotation_matrix([0., 0., 1.], π + swimming_angle)
#     R = rotation_align_to_x(swimming_direction)
#     sprob = prob.swimming_problem

#     for t in eachindex(traj.t)
#         # move_boundary!(sprob, traj.x[t], traj.b1[t], traj.b2[t], traj.t[t])
#         # centroid = mean(eachcol(sprob.points.quad_pts))
#         # R_body = hcat(traj.b1[t], traj.b2[t], cross(traj.b1[t], traj.b2[t]))
#         # move_boundary!(sprob, R*(traj.x[t] - centroid), R*R_body, traj.t[t])
#         move_boundary!(sprob, traj.x[1], traj.b1[t], traj.b2[t], traj.t[t])
#         solve_problem!(sprob)
#         push!(vfs, VelocityField(sprob, x_points, y, z_points))
#     end
#     VelocityField(vfs[1].points, mean(vf.velocities for vf in vfs))
# end





