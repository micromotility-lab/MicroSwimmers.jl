# body shape parameter - Bretherton number based on the terms of the body shape tensor B, which is the contraction of the grand mobility matrix and the shear tensor.
#if the coefficients corresponding to the Bretherton number are not consistent, then the Bretherton number is set to NaN and a warning is printed.

function body_shape_parameter(ms::MicroSwimmer; eps = 0.1)   
    B = body_shape_tensor(ms; eps = eps)
    if abs(abs(B[5,1,3])-abs(B[6,1,2]))< 1e-1
    Breth =  B[6,1,2]- B[5,1,3]
    else
        Breth = NaN
        println("not consistent with Bretherton")
    end
    Breth
end
function body_shape_tensor(ms::MicroSwimmer; eps = 0.1)   
    R = grand_resistance_matrix(ms,eps = eps)
    K = R^-1
    G = force_and_torque_shear(ms,eps = eps)
    G_mat = reshape(G, (6,9))
    Bmat = K*G_mat
    B = reshape(Bmat, (6,3,3))
    B
end

function body_shape_oscillation(ms::MicroSwimmer; period=1.0, num_ts=30, eps = 0.1)
    Breths = Float64[]
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(ms, t)
        Breths = [Breths; body_shape_parameter(ms; eps = eps)]
    end
    Breths
end

function average_body_shape(ms::MicroSwimmer; period=1.0, num_ts=30, eps = 0.1)
    Breths = body_shape_oscillation(ms; period=period, num_ts=num_ts, eps = eps)
    sum(Breths)/num_ts
end
