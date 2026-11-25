#consider higher order moments of the force distribution by expanding about a point x0
struct MultipoleMoments
    F::SVector{3,Float64}
    M::SMatrix{3,3,Float64,9}
    Q::Array{Float64,3}
    O::Array{Float64,4}
    E::SVector{3,Float64}
    mu::Float64
end

function multi_moments(prob, x0)
    disc = prob.disc; f = get_forces(prob)
    F = zero( SVector{3,Float64}); M = zeros(3, 3); Q = zeros(3, 3, 3); O = zeros(3, 3, 3, 3); E = zero( SVector{3,Float64})
    epsq = zeros(length(disc.quad_pts))
    for (rng, ε) in zip(disc.quad_part_ranges, disc.quad_eps); epsq[rng] .= ε; end
    for q in eachindex(disc.quad_pts)
        w = disc.quad_wts === nothing ? 1.0 : disc.quad_wts[q]
        fq = w * f[disc.nearest[q]]; y = disc.quad_pts[q] - x0
        F += fq; E += epsq[q]^2 * fq
        for j in 1:3, k in 1:3
            M[j, k] += fq[j] * y[k]
            for l in 1:3
                Q[j, k, l] += fq[j] * y[k] * y[l]
                for m in 1:3; O[j, k, l, m] += fq[j] * y[k] * y[l] * y[m]; end
            end
        end
    end
    MultipoleMoments(F, SMatrix{3,3}(M), Q, O, E, prob.mu)
end

#basis of far-field singularities

function GD(x, d, e)   # force dipole                         
    r = norm(x)
    ((dot(d, x) * e - dot(e, x) * d - dot(d, e) * x) / r + 3dot(e, x) * dot(d, x) * x / r^3) / r^2
end

function GQ(x, c, d, e)    #quadrupole                                                             
    r = norm(x); cx, dx, ex = dot(c, x), dot(d, x), dot(e, x)
    g1 = (dot(d, e) * cx + dot(c, e) * dx + dot(c, d) * ex) * x + dx * ex * c + cx * ex * d - cx * dx * e
    (dot(d, e) * c + dot(c, e) * d - dot(c, d) * e - 3g1 / r^2 + 15cx * dx * ex * x / r^4) / r^3
end

function GO(x, b, c, d, e)    # octupole                                                            
    r = norm(x); bx, cx, dx, ex = dot(b, x), dot(c, x), dot(d, x), dot(e, x)
    bc, bd, be, cd, ce, de = dot(b, c), dot(b, d), dot(b, e), dot(c, d), dot(c, e), dot(d, e)
    g2 = (cd * ex + ce * dx + de * cx) * b + (de * bx + be * dx + bd * ex) * c +
         (ce * bx + be * cx + bc * ex) * d - (cd * bx + bd * cx + bc * dx) * e +
         (cd * be + ce * bd + bc * de) * x
    g3 = cx * dx * ex * b + bx * dx * ex * c + bx * cx * ex * d - bx * cx * dx * e +
         (be * cx * dx + bd * cx * ex + bc * dx * ex + cd * bx * ex + ce * bx * dx + de * bx * cx) * x
    3 / r^4 * (g2 / r - 5g3 / r^3 + 35bx * cx * dx * ex * x / r^5)
end
#source dipole
Dsrc(x, e) = (r = norm(x); (-e + 3dot(e, x) * x / r^2) / r^3)

#source quadrupole
function Qsrc(x, d, e)                                                                  
    r = norm(x); dx, ex = dot(d, x), dot(e, x)
    -3 / r^4 * ((dx * e + ex * d + dot(d, e) * x) / r - 5ex * dx * x / r^3)
end
#rotlet dipole                                            
RD(x, d, e) = cross(d, e) / norm(x)^3 + 3dot(d, x) * cross(e, x) / norm(x)^5         

# far-field expansion at each order from the moments, using the basis of singularities
#quadrupole is symmetric & traceless part of second moment, source dipole is the trace part of the second moment, rotlet dipole is the skew part of the second moment
const EB = (ex, ey, ez)

ff2(m::MultipoleMoments, x) = sum(m.M[j, k] * GD(x, EB[k], EB[j]) for j in 1:3, k in 1:3) / (8π * m.mu)
ff3(m::MultipoleMoments, x; reg=true) =
    (0.5sum(m.Q[j, k, l] * GQ(x, EB[l], EB[k], EB[j]) for j in 1:3, k in 1:3, l in 1:3) -
     (reg ? 0.5Dsrc(x, m.E) : zero(SVector{3,Float64}))) / (8π * m.mu)
ff4(m::MultipoleMoments, x) =
    sum(m.O[j, k, l, n] * GO(x, EB[n], EB[l], EB[k], EB[j]) for j in 1:3, k in 1:3, l in 1:3, n in 1:3) / (6 * 8π * m.mu)
farfield(m::MultipoleMoments, x; order=3, reg=true) =
    Gs(x, m.F) / (8π * m.mu) + ff2(m, x) + (order >= 3 ? ff3(m, x; reg=reg) : zero(SVector{3,Float64})) +
    (order >= 4 ? ff4(m, x) : zero(SVector{3,Float64}))



