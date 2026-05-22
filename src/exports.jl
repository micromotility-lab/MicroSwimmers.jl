# Utilities 
export rotation_matrix, rotation_align_to_x
export hf, hq
# Core structures & problems
export FluidBoundary, update_boundary!, move_boundary!
export CellBody, grand_resistance_matrix
export Configuration

export fibonacci_cylinder
# Bodies
export RigidMotionBody, add_rigid_body_motion!, add_velocity!, add_angular_velocity!, reset_velocity!
export SphericalBody, EllipsoidBody, EllipsoidalGroovedBody, CylindricalGroovedBody, FlatGroovedBody #, PNASExcavateBody
    
# Flagella
export Flagellum, BareFlagellum, TubeFlagellum, LineTubeFlagellum,VanedFlagellum 
export FlagellumModel, PlanarFlagellum, StandingWaveFlagellum, QuasiPlanarFlagellum, ThreeDimensionalFlagellum#, TubePlanarFlagellum #, ThreeDimensionalFlagellum, StandingWaveFlagellum, Vane, VanedFlagellum
export get_vane_pts
# export get_pts!, get_pts_and_velocity!
# Swimmers
export UniFlagellate, Flagellate, MicroSwimmer
export discretisation, NearestDiscretisation, NystromDiscretisation

export Colony
# export move!, update!

export is_inside_ellipsoid

export regularised_blakelet!

# Problems
export SwimmingProblem, SwimmingTrajectoryProblem, get_U, get_Ω, get_force_pts, get_forces, get_velocities, get_quad_pt_velocities #,  FeedingProblem, MultipleSwimmerProblem
export ResistanceProblem, ParticleTrajectoryProblem
export solve_problem! #, total_body_force_and_torque, get_velocity_function
export translate_problem!, rotate_problem!

# Trajectories
export Trajectory, average_swimming_velocity, continue_periodic_trajectory!, continue_periodic_trajectory, running_mean, centred_trajectory, translate_trajectory
export Helix, helix, fit_helix, pitch, pitch_angle, curvature, torsion, chirality_sign, axis_angular_velocity, axis_azimuthal_angle, axis_velocity, axis_polar_angle, radius, axis_vector, initial_point, translate_helix
export initial_helix_pars


# Quantification
export mean_std
export total_force, total_torque, total_power, total_force_and_torque, stresslet_tensor, average_stresslet_tensor, disturbance_stresslet_tensor, total_energy_dissipated
export VelocityFunction, FluidVelocity, AverageVelocityFunction, VelocityField, PlanarVelocityField, get_velocity_function, velocity_flux, velocity_flux_polar, TimeAveragedPlanarVelocityField, TimeAveragedDisturbanceField

export spacing

export get_sol!

export check_boundary_conditions, check_body_boundary_conditions
