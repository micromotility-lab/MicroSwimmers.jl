# Core structures & problems
export FluidBoundary, update_boundary!, CellBody # , Flagellum, Swimmer, Problem
export Configuration

# Bodies
export SphericalBody, EllipsoidBody, EllipsoidalGroovedBody # , CylindricalGroovedBody, PNASExcavateBody

# Flagella
export Flagellum, PlanarFlagellum, QuasiPlanarFlagellum, TubePlanarFlagellum #, ThreeDimensionalFlagellum, StandingWaveFlagellum, Vane, VanedFlagellum
# export get_pts!, get_pts_and_velocity!
# Swimmers
export UniFlagellate, Flagellate
# export move!, update!

# Problems
export SwimmingProblem, DynamicSwimmingProblem #,  FeedingProblem, MultipleSwimmerProblem
export solve_problem! #, total_body_force_and_torque, get_velocity_function

# Utilities you want public (be selective)
# export nearest_neighbour, fibonacci_ellipsoid, ...
