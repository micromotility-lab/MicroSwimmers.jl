# DO THIS ALL 3D and then just plot the projection

function get_velocity_function(prob::SwimmingProblem; body_frame=false, xy_plane=false)
    check_solved!(prob)    
    N = length(prob.points.force_pts)
    A = zeros(3, N)
    U = get_U(prob)
    Ω = get_Ω(prob)
    
    function vel(x)
        if xy_plane
            resistance_matrix!(
                A, [x[1], x[2], 0.],
                prob.points.quad_pts,
                prob.points.nearest,
                prob.eps;
                μ=prob.mu
            )
        else
            resistance_matrix!(
                A, x, 
                prob.points.quad_pts, 
                prob.points.nearest, 
                prob.eps;
                μ=prob.mu,
            )
        end

        u = A * prob.force_vals[1:N]     
        if body_frame
            u -=  U + cross(Ω, [x[1], x[2], 0.] - prob.swimmer.points.location)
        end
        
        if xy_plane
            return SVector{2}(u[1:2])
        else
            return SVector{3}(u)
        end
    end
end

function get_velocity_function(prob::Problem)
    solve_problem!(prob)
    N = prob.points.N
    A = zeros(3, 3N)
    
    function u(x)
        resistance_matrix!(
                A, x, 
                prob.points.quad_pts, 
                prob.points.nearest, 
                prob.eps;
                μ=prob.mu,
        )
        SVector{3}(A * prob.force_vals[1:3N])
    end
end



## xy plane, MAYBE UNNECESSARY, A LOT OF REDUNDANCY HERE 
struct VelocityField{T <: Number}
    points::Union{Vector{SVector{2,T}}, Vector{SVector{3, T}}}
    velocities::Union{Vector{SVector{2,T}}, Vector{SVector{3,T}}}
end

function VelocityField(prob::SwimmingProblem{T}, points::Matrix{T}; body_frame=false) where {T <: Number}
    check_solved!(prob)
    planar = size(points, 1) == 2
    static_points = planar ? [SVector{2}(pt) for pt in eachcol(points)] : [SVector{3}(pt) for pt in eachcol(points)]
    u = get_velocity_function(prob, body_frame=body_frame, xy_plane=planar)
    VelocityField(static_points, u.(eachcol(points)))
end

function VelocityField(prob::ResistanceProblem{T}, points::Matrix{T}) where {T <: Number}
    static_points = [SVector{3}(pt) for pt in eachcol(points)]
    u = get_velocity_function(prob)
    VelocityField(static_points, u.(eachcol(points)))
end

function VelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y_points::AbstractVector{T}, z::T) where {T <: Number}
    pts = reduce(hcat, vec([[x,y,z] for x in x_points, y in y_points]))
    VelocityField(prob, pts)
end


function TimeAveragedVelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, pts::Matrix{T}; period=1.0, num_ts=30) where {T <: Number}
    vfs = VelocityField[]
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(prob, t)
        solve_problem!(prob)
        push!(vfs, VelocityField(prob, pts))
    end
    VelocityField(vfs[1].points, sum(vf.velocities for vf in vfs) / num_ts)
end

function TimeAveragedVelocityField(prob::Union{ResistanceProblem, SwimmingProblem}, x_points::AbstractVector{T}, y_points::AbstractVector{T}, z::T, ; period=1.0, num_ts=30) where {T <: Number}
    vfs = VelocityField[]
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(prob, t)
        solve_problem!(prob)
        push!(vfs, VelocityField(prob, x_points, y_points, z))
    end
    VelocityField(vfs[1].points, sum(vf.velocities for vf in vfs) / num_ts)
end





