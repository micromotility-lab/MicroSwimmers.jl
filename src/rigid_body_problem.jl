# body shape parameter - Bretherton number based on the terms of the body shape tensor B, which is the contraction of the grand mobility matrix and the shear tensor.
#if the coefficients corresponding to the Bretherton number are not consistent, then the Bretherton number is set to NaN and a warning is printed.
function body_shape_parameter(ms::MicroSwimmer; eps = nothing)
    B = body_mobility_tensor(ms; eps = eps)
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
    ms = body_fixed(ms)
    R = grand_resistance_matrix(ms,eps = eps)
    K = R^-1
    G = force_and_torque_shear(ms,eps = eps)
    G_mat = reshape(G, (6,9))
    Bmat = K*G_mat
    B = reshape(Bmat, (6,3,3))
    B
end

function body_shape_oscillation(ms::MicroSwimmer; period=1.0, num_ts=30, eps = nothing)
    Breths = Float64[]
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(ms, t)
        push!(Breths, body_shape_parameter(ms; eps = eps))
    end
    Breths
end

# function to ensure that the body frame is fixed in body so resistance matrix is constant and can be used to compute the body shape tensor
function body_fixed(ms::MicroSwimmer)
    parts = [Part(p.model, p.disc, p.frame) for p in ms.parts]
    MicroSwimmer(parts)    
end

#rigid body problem  for  - angular velocity of the body frame in a shear flow with constant shear
#depends on the  body shape tensor B which is constant in the body frame
#tensor contract B_torque with rate of strain of background flow, use body frame basis to simplify calculation 
function orientation_ode(X,B_torque,strain_shear,Omega_shear,t)
    b1 = X[1:3]
    b2 = X[4:6]
    b3 = cross(b1,b2)
    Q = hcat(b1,b2,b3) 
    Omega_b1 =sum((Q*B_torque[1,:,:]*Q').*strain_shear)
    Omega_b2 =sum((Q*B_torque[2,:,:]*Q').*strain_shear)
    Omega_b3 =sum((Q*B_torque[3,:,:]*Q').*strain_shear)
    db1 = cross(Omega_shear,b1) - Omega_b2*b3 + Omega_b3*b2
    db2 = cross(Omega_shear,b2) + Omega_b1*b3 - Omega_b3*b1
    return SVector{6,eltype(X)}(db1..., db2...)
end


function RigidBodyProblem(ms::MicroSwimmer; 
    Omega_shear = SVector(0.5, 0.0, 0.0),
    strain_shear = SMatrix{3,3}(0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0),
    B_torque = body_mobility_tensor(ms)[4:6, :, :],
    t_final=20.0,
    saveat=0.05,
    eps=nothing
    )
    warn_problem_eps(eps)

    T = Float64
    X0 = SVector{6,T}(ms.frame.orientation[:,1]..., ms.frame.orientation[:,2]...)
    rhs = (X, p, t) -> orientation_ode(X, B_torque, strain_shear, Omega_shear, t)
    sol = solve(ODEProblem(rhs, X0, (0.0, t_final)), Tsit5(), saveat=saveat)
    u = sol.u
    x = [ms.frame.location for i in eachindex(u)]
    b1 = [SVector{3}(u[i][1:3]...) for i in eachindex(u)]
    b2 = [SVector{3}(u[i][4:6]...) for i in eachindex(u)]
    return sol, x, b1, b2
end

