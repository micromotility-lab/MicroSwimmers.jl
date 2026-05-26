# Utilities 
export rotation_matrix
export hf, hq


export FluidBoundary, update_boundary!, move_boundary!
export CellBody, grand_resistance_matrix

export fibonacci_ellipsoid, fibonacci_cylinder, is_inside_cylinder, is_inside_ellipsoid

# Bodies
export SphericalBody, EllipsoidBody, EllipsoidalGroovedBody, CylindricalGroovedBody, FlatGroovedBody
export add_rigid_body_motion!
    
# Flagella
export Flagellum, BareFlagellum, VanedFlagellum 
export FlagellumModel, PlanarFlagellum, StandingWaveFlagellum
export get_vane_pts

# Swimmers
export UniFlagellate, Flagellate, MicroSwimmer
export discretisation, NearestDiscretisation


# Problems
export SwimmingProblem, SwimmingTrajectoryProblem, get_U, get_Ω, get_force_pts, get_forces, get_velocities, get_quad_pt_velocities
export ResistanceProblem, ParticleTrajectoryProblem
export solve_problem!
export translate_problem!, rotate_problem!

# Trajectories
export Trajectory, average_swimming_velocity, continue_periodic_trajectory, running_mean, centred_trajectory, translate_trajectory
export Helix, helix, fit_helix, pitch, pitch_angle, curvature, torsion, chirality_sign, axis_angular_velocity, axis_azimuthal_angle, axis_velocity, axis_polar_angle, radius, axis_vector, initial_point, translate_helix
export initial_helix_pars


# Quantification
export total_force, total_torque, total_power, total_force_and_torque, stresslet_tensor, average_stresslet_tensor, total_energy_dissipated
export VelocityFunction, FluidVelocity, AverageVelocityFunction, VelocityField, PlanarVelocityField, get_velocity_function, velocity_flux, velocity_flux_polar, TimeAveragedPlanarVelocityField, TimeAveragedDisturbanceField


export get_sol

export check_boundary_conditions, check_body_boundary_conditions
