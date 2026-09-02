# Every aggregator below sums over quadrature points, matching the force-free/torque-free
# rows that assemble_swimming! builds. quad_weight(nothing, q) is 1, so the unweighted case
# reduces exactly to the plain patch-multiplicity sums these were before weights existed.
total_force(forces, disc::NearestDiscretisation) =
    sum(quad_weight(disc.quad_wts, q) * forces[disc.nearest[q]] for q in eachindex(disc.nearest))

function total_torque(forces, disc::NearestDiscretisation)
    @unpack quad_pts, nearest, quad_wts = disc
    sum(quad_weight(quad_wts, q) * cross(quad_pts[q], forces[nearest[q]]) for q in eachindex(quad_pts))
end

total_force(forces, ::NystromDiscretisation) = sum(forces)

function total_torque(forces::AbstractVector{<:SVector}, disc::NystromDiscretisation)
    sum(cross(disc.force_pts[i], forces[i]) for i in eachindex(disc.force_pts))
end

function total_force_and_torque(prob::InstantaneousProblem)
    check_solved!(prob)
    forces = get_forces(prob)
    total_force(forces, prob.disc), total_torque(forces, prob.disc)
end

function average_stresslet_tensor(prob::InstantaneousProblem; period=1.0, num_ts=30)
    check_solved!(prob)
    Ss = []
    for t in range(0, period, num_ts)[1:end-1]
        update_boundary!(prob, t)
        solve_problem!(prob)
        push!(Ss, stresslet_tensor(prob))
    end
    mean(Ss)
end

# function total_power(prob::InstantaneousProblem)
#     check_solved!(prob)
#     forces = get_forces(prob)
#     vels   = prob.disc.velocity
#     sum(dot(forces[n], vels[n]) for n in eachindex(forces))
# end


function total_power(prob::InstantaneousProblem)
    check_solved!(prob)
    _total_power(get_forces(prob), prob.disc)
end

# Rate of working of the surface tractions, P = sum_q w_q <f_{nearest[q]}, u_{nearest[q]}>.
# The sum runs over quadrature points, like total_force and total_torque: summing over force
# points instead drops the patch multiplicity, which is the bug fixed in 7b3326f (introduced
# by fecbebf during the old-API port, and shipped in the v0.2.0 tag).
function _total_power(forces, disc::NearestDiscretisation)
    @unpack nearest, velocity, quad_wts = disc
    sum(quad_weight(quad_wts, q) * dot(forces[nearest[q]], velocity[nearest[q]])
        for q in eachindex(nearest))
end

_total_power(forces, disc::NystromDiscretisation) =
    sum(dot(forces[i], disc.velocity[i]) for i in eachindex(forces))

function total_energy_dissipated(prob::SwimmingTrajectoryProblem)
    check_solved!(prob)
    traj = prob.traj
    sprob = prob.swimming_problem

    Es = Float64[]

    for (i, t) in enumerate(traj.t)
        move_boundary!(sprob, traj.x[i], traj.b1[i], traj.b2[i], t)
        solve_problem!(sprob)
        push!(Es, total_power(sprob))
    end
    # Trapezoidal integration
    sum(0.5 * (Es[i] + Es[i+1]) * (traj.t[i+1] - traj.t[i]) for i in 1:length(traj.t)-1)
end

# First moment of the surface traction, M = Σ_q w_q f_q p_q'.
# Symmetric traceless part is the stresslet; the skew part carries the torque.
function shear_tensor(forces, disc::NearestDiscretisation)
    @unpack quad_pts, nearest, quad_wts = disc
    sum(quad_weight(quad_wts, q) * (forces[nearest[q]] * quad_pts[q]')
        for q in eachindex(quad_pts))
end

shear_tensor(forces, disc::NystromDiscretisation) =
    sum(forces[i] * disc.force_pts[i]' for i in eachindex(disc.force_pts))

_stresslet_from_shear(M) = 0.5 * (M + M') - (1/3) * tr(M) * I
_stresslet_tensor(forces, disc) = _stresslet_from_shear(shear_tensor(forces, disc))

function stresslet_tensor(prob::InstantaneousProblem)
    check_solved!(prob)
    _stresslet_tensor(get_forces(prob), prob.disc)
end

function force_and_torque_shear(ms::MicroSwimmer; eps=nothing)
    Gamma = zeros(6,3,3)
    prob = ResistanceProblem(ms)
    for (i, n) in enumerate([ex, ey, ez]) 
        add_rigid_body_motion!!(prob.microswimmer, n, zero(SVector{3,Float64}))
        solve_problem!(prob)
        Gamma[i,:,:] = stresslet_tensor(prob)
    
        add_rigid_body_motion!!(prob.microswimmer, zero(SVector{3,Float64}), n)
        solve_problem!(prob)
        Gamma[3+i,:,:] = stresslet_tensor(prob)
        #symmetry and traceless of the rate of strain means that we can make shear tensor symmetric and traceless too
    end
    Gamma
end