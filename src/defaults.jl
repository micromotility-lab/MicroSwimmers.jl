# Defaults for the model-only `Part(model)` constructor. Lengths are microns: everything in
# `MicroSwimmersPlots/examples` works in μm, so a 1.0 here is 1 μm.
#
# `eps` plays two different roles depending on what it sits on. On a cell body it is a purely
# numerical regularisation of the singular Stokeslet; on a slender flagellum it is the physical
# radius of the filament. Both are served by the same default because what constrains the
# discretisation is the ratio hq/eps, not eps itself.
#
# hq = 2*eps: the quadrature has to resolve the blob it is integrating. Accuracy is insensitive
# to this ratio until roughly hq ~ 5*eps, so 2 leaves a 2.5x margin while keeping Q/N at ~5 on a
# flagellum and ~25 on a surface — comfortably inside the "use Q > 4N" rule of thumb, and close
# to what the hand-tuned examples already use.
const DEFAULT_HF        = 0.5    # target force-point spacing, μm
const DEFAULT_EPS       = 0.05   # regularisation parameter, μm
const DEFAULT_HQ_FACTOR = 2.0    # hq = DEFAULT_HQ_FACTOR * eps

# Cross-section point counts for the tube flagella (src/flagellum_tubes.jl).
#
# A quadrature ring is a periodic trapezoid rule over a smooth periodic integrand: the
# collocation point sits at distance exactly `radius` from every point of its own ring, so
# there is no near-singularity to resolve and convergence is spectral. This is NOT the
# 2*pi*a/eps that resolving a blob on the ring would cost.
#
# Measured on a curved filament (L=10, a=0.1, N=61), relative error in the translating drag
# against a Q_cs=48 reference: Q_cs=4 -> 8e-6, 6 -> 5e-9, 8 -> 3e-12, 12 -> machine precision.
# 12 is past converged and cheap, so it is the default. Do not read that off a *straight*
# filament: there all the rings are congruent and the symmetry makes Q_cs=4 exact.
const DEFAULT_Q_CS = 12   # quadrature points per cross-section
const DEFAULT_N_CS = 8    # force points per cross-section (SurfaceTubeFlagellum only)
