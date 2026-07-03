# MicroSwimmers.jl

MicroSwimmers.jl is a Julia package for simulating low-Reynolds number swimmers like protists and spermatozoa. The core method is an implementation of the regularised-stokeslet boundary-element method with nearest-neighbour discretisation ([Smith 2009](https://doi.org/10.1016/j.jcp.2017.12.008)). The package also contains tools for constructing generic microswimmers through parameterised body and cilia/flagella models, calculating fluid flows and quantifying outputs such as swimming trajectories and fluxes. Visualisation tools for MicroSwimmers.jl simulations can be found in the separate package [MicroSwimmersPlots.jl](https://github.com/micromotility-lab/MicroSwimmersPlots.jl).

This package is developed by James Cass as a postdoc in the [micromotility lab](https://micromotility.com/) led by Kirsty Wan, in the University of Exeter's Living Systems Institute.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/micromotility-lab/MicroSwimmers.jl")
```

## Quick start

```julia
using MicroSwimmers

# # A spherical cell body (semi-axes a = b = c = 1 μm)
body = EllipsoidBody(1.0, 1.0, 1.0),

# discretise the body with N = 213 force points and Q=917 quadrature points
body_disc = Part(body, 213, 917)

# Define a planar flagellum beating pattern (tangent-angle model, Gallagher et al. 2018):
#   θ(s,t) = Cs + (R₀ + R₁ sin(ks/L)) cos(ωt - ϕs/L)
flagellum = PlanarFlagellum(
    10.0,  # L:  length (μm)
    0.0,   # C:  static curvature
    0.6,   # R₀: amplitude envelope
    0.5,   # R₁: spatial modulation of amplitude
    π/2,   # k:  envelope wavenumber
    2π,    # ϕ:  travelling-wave wavenumber
    2π,    # ω:  angular frequency
    0.0,   # δ:  overall phase
)

# discretise the flagellum, attached at the edge of the body on the x-axis
flagellum_disc = Part(
    flagellum, 23, 127,
    location=[1.0, 0.0, 0.0],
    orientation=rotation_matrix([1.0, 0.0, 0.0], 0.0),
)

# Assemble the swimmer from its parts
ms = MicroSwimmer([body_disc, flagellum_disc])

# Solve the swimming problem for the rigid-body velocity U, angular velocity Ω,
# and the force distribution
prob = SwimmingProblem(ms)
solve_problem!(prob)

U      = get_U(prob)
Ω      = get_Ω(prob)
forces = get_forces(prob)

# An isolated swimmer has zero net-force and net-torque
F, T = total_force_and_torque(prob)
```

This builds a sperm-like swimmer — a cell body with an attached beating flagellum — solves for how it swims, and confirms the net force and torque vanish.

## Related packages

- [MicroSwimmersPlots.jl](https://github.com/micromotility-lab/MicroSwimmersPlots.jl) — visualisation for MicroSwimmers.jl
- [MicroSwimmersExamples.jl](link) — worked examples (coming soon)

## Citation

If you use this package in your research, please cite:

```bibtex
@article{cass2026simulation,
  title={Simulation-driven discovery of morphology-function relationships in microswimmers},
  author={Cass, James F and Wan, Kirsty Y},
  journal={bioRxiv},
  pages={2026--06},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```

## License

This project is licensed under the [MIT License](LICENSE).
