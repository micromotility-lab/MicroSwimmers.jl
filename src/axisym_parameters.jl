#want to calculate the parameters of minimal models for the swimmer, which is assumed to be axisymmetric
#determine singularity strengths of force dipole, source dipole, quadrupole and rotlet dipole strengths 

# defining axisymmetric basis force dipole + source dipole + quadrupole + rotlet dipole
axi_basis(d) = (x -> GD(x, d, d), x -> Dsrc(x, d), x -> GQ(x, d, d, d), x -> RD(x, d, d))

# axisymmetric basis for O(r^-4) terms
axi_basis4(d) = (x -> GO(x, d, d, d, d), x -> Qsrc(x, d, d))

#helper functions for least squares regression, use design to create matrix for lsq regression
function design(fns, pts, x0, w=ones(length(pts)))
    A = zeros(3length(pts), length(fns))
    for (p, x) in enumerate(pts), (j, f) in enumerate(fns)
        A[3p-2:3p, j] = w[p] * f(x - x0)
    end
    A
end
#stack 3 vectors for flow field at each point into a single vector, scaled by w[p]
stack3(vs, w=ones(length(vs))) = reduce(vcat, [w[p] * vs[p] for p in eachindex(vs)])

# minimum-norm least squares via truncated SVD (bases above are rank deficient by design)
function lsq(A, b; rtol=1e-10)
    F = svd(A); k = count(>(rtol * F.S[1]), F.S)
    F.V[:, 1:k] * ((F.U[:, 1:k]' * b) ./ F.S[1:k])
end

#least squares regression of O(r-2) and O(r^-3) singularity strengths from exact moments, using axisymmetric basis about d
# also returns fraction of each order that is NOT axisymmetric about d, for checking if the swimmer is approximately axisymmetric
function exact_axi(m::MultipoleMoments, d; reg=true)
    UNIT_SPHERE = fibonacci_ellipsoid(1.0, 1.0, 1.0, 3000)
    b = axi_basis(d)
    u2 = [ff2(m, x) for x in UNIT_SPHERE]; u3 = [ff3(m, x; reg=reg) for x in UNIT_SPHERE]
    A2 = design([b[1]], UNIT_SPHERE, zero( SVector{3,Float64})); A3 = design(collect(b[2:4]), UNIT_SPHERE, zero( SVector{3,Float64}))
    c2 = A2 \ stack3(u2); c3 = A3 \ stack3(u3)
    # fraction of each order that is NOT axisymmetric about d
    na2 = norm(A2 * c2 - stack3(u2)) / norm(stack3(u2)); na3 = norm(A3 * c3 - stack3(u3)) / norm(stack3(u3))
    (α=c2[1], β=c3[1], γ=c3[2], τ=c3[3], nonaxi2=na2, nonaxi3=na3)
end

#returns strengths of the r⁻² and r⁻³ singularities about d
# also returns fraction of each order that is NOT axisymmetric about d, for checking if the swimmer is approximately axisymmetric
function singularity_strengths(ms; t = 0.0, d = ms.frame.orientation[:,1], x0 = ms.frame.location, prob = nothing)
    prob = prob === nothing ? SwimmingProblem(ms) : prob
    update_boundary!(prob, t); solve_problem!(prob)
    c = exact_axi(multi_moments(prob, SVector{3}(x0)), normalize(SVector{3}(d)))
    [c.α, c.β, c.γ, c.τ], c.nonaxi2, c.nonaxi3
end

# want to determine Bretherton parameter for an instaneous swimmer, from body mobility tensor 

function body_mobility_tensor(ms::MicroSwimmer; eps = nothing, t= 0.0)   
    move_boundary!(ms,SVector(0.,0.,0.),  SVector(1.0,0.,0.), SVector(0.,1.,0.), t)
    R = grand_resistance_matrix(ms)
    R = (R+R')/2 
    G = force_and_torque_shear(ms)
    G_mat = reshape(G, (6,9))
    Bmat = R\G_mat
    B = reshape(Bmat, (6,3,3))
    B
end

# helper functions
#calculate a fixed basis from a given axis in order to calculate bretherton parameter
function basis_from_axis(d; normal = SVector(0.0, 0.0, 1.0))
    d = normalize(SVector{3}(d)); n_out = normalize(normal - dot(normal, d)*d)
    hcat(d, cross(n_out, d), n_out)
end
#express body tensor T in new basis V 
rotate_tensor3(T, V) = [sum(V[a,i]*V[b,j]*V[c,k]*T[a,b,c] for a in 1:3, b in 1:3, c in 1:3) for i in 1:3, j in 1:3, k in 1:3]

#generalised Jeffery equation for rotation of rigid body in a shear flow, with constant shear rate and rate of strain, expressed in the lab frame
function orientation_ode(X,B_torque,strain_shear,Omega_shear)
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

#calculate bretherton parameter by lsq fit of generalised jeffery to jeffery over basis of strain tensors
function axisym_shape_parameter(ms::MicroSwimmer; t = 0.0, V = I(3),
        orientation = (SVector(0.8,0.0,0.6), SVector(-0.6,0.0,0.8)), C = nothing)
     STRAIN_BASIS = let s = 1/sqrt(2)
        [SMatrix{3,3}(0, s, 0, s, 0, 0, 0, 0, 0), SMatrix{3,3}(0, 0, s, 0, 0, 0, s, 0, 0),
        SMatrix{3,3}(0, 0, 0, 0, 0, s, 0, s, 0), SMatrix{3,3}(s, 0, 0, 0, -s, 0, 0, 0, 0),
        SMatrix{3,3}(Diagonal([1.0, 1.0, -2.0]/sqrt(6)))]
     end
    C = C === nothing ? [body_mobility_tensor(ms; t = t)[3+i, j, k] for i in 1:3, j in 1:3, k in 1:3] : C   # body frame
    V = Matrix{Float64}(V); det(V) < 0 && (V[:,3] .*= -1)
    B_torque = rotate_tensor3(C, V)                    # rotation tensor in the V basis
    b1, b2 = SVector{3}(orientation[1]), SVector{3}(orientation[2])
    Q = hcat(b1, b2, cross(b1, b2))                    # body (V) axes in the lab
    X0 = SVector{6}(b1..., b2...); e = b1
    strain_of_e(E) = (de = orientation_ode(X0, B_torque, E,  SVector(0.0, 0.0, 0.0)); de[1:3])
    jeff(E) = E*e - dot(e, E*e)*e
    num = 0.0; den = 0.0; tot = 0.0
    for E in STRAIN_BASIS
        lhs = strain_of_e(E); rhs = jeff(E)
        num += dot(rhs, lhs); den += dot(rhs, rhs); tot += dot(lhs, lhs)
    end
    B = num/den
    err = sqrt(max(sum(norm(strain_of_e(E) - B*jeff(E))^2 for E in STRAIN_BASIS), 0.0)/tot)
    return (B = B, err = err)
end

#function that assigns Bretherton number based on body mobility tensor,
function body_shape_parameter(ms::MicroSwimmer; eps = nothing, t= 0.0)
    B = body_mobility_tensor(ms; t = t)
    if abs((B[5,1,3])+(B[6,1,2]))< 1e-1
    Breth =  B[6,1,2]- B[5,1,3]
    else
        Breth = NaN
    end
    Breth
end