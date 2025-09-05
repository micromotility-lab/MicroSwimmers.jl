# Core structures & problems
export FluidBoundary, update_boundary!, move_boundary!
export CellBody
export Configuration

# Bodies
export SphericalBody, EllipsoidBody, EllipsoidalGroovedBody # , CylindricalGroovedBody, PNASExcavateBody

# Flagella
export Flagellum, TubeFlagellum, PlanarFlagellum, QuasiPlanarFlagellum, TubePlanarFlagellum #, ThreeDimensionalFlagellum, StandingWaveFlagellum, Vane, VanedFlagellum
# export get_pts!, get_pts_and_velocity!
# Swimmers
export UniFlagellate, Flagellate, Swimmer
# export move!, update!

# Problems
export SwimmingProblem, DynamicSwimmingProblem #,  FeedingProblem, MultipleSwimmerProblem
export solve_problem! #, total_body_force_and_torque, get_velocity_function

# Trajectories
export Trajectory, swimming_velocity, continue_periodic_trajectory!


export total_force, total_torque, total_power, total_force_and_torque, get_velocity_function