function get_velocity_function(prob::SwimmingProblem)
    check_solved!(prob)    
    N = length(prob.points.force_pts)
    A = zeros(3, N)
    
    function u(x)
        resistance_matrix!(
            A, 
            reshape(x,3,1), 
            prob.points.quad_pts, 
            prob.points.nearest, 
            prob.eps;
            μ=prob.mu,
        )
        SVector{3}(A * prob.force_vals[1:N])
    end
end
