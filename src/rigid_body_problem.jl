# body shape parameter - Bretherton number based on the terms of the body shape tensor B, which is the contraction of the grand mobility matrix and the shear tensor.
#if the coefficients corresponding to the Bretherton number are not consistent, then the Bretherton number is set to NaN and a warning is printed.
function body_shape_parameter(ms::MicroSwimmer; eps = nothing)
    B = body_mobility_tensor(ms)
    if abs((B[5,1,3])+(B[6,1,2]))< 1e-1
    Breth =  B[6,1,2]- B[5,1,3]
    else
        Breth = NaN
        #println("not symmetric enough")
    end
    Breth
end
#use grand mobility tensor and shear-force/torque tensor to determine mobitly tensor in body fixed frame, 
#can then be used to calculate rigid body motion byor determine number of independent coefficients in mobility problem
function body_mobility_tensor(ms::MicroSwimmer; eps = nothing)   
    move_boundary!(ms,SVector(0.,0.,0.),  SVector(1.0,0.,0.), SVector(0.,1.,0.), 0.0)
    R = grand_resistance_matrix(ms)
    K = R^-1
    G = force_and_torque_shear(ms)
    G_mat = reshape(G, (6,9))
    Bmat = K*G_mat
    B = reshape(Bmat, (6,3,3))
    B
end

function body_shape_oscillation(ms::MicroSwimmer; period=1.0, num_ts=30, eps = nothing)
    Breths = Float64[]
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(ms, t)
        push!(Breths, body_shape_parameter(ms))
    end
    Breths
end


