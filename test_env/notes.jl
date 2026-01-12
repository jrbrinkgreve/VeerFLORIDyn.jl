#=
notes for veer implementations




...
runFLORIS funcs:
        - prepare_rotor_points!
        - handle_single_turbine!
        - setup_computation_buffers!
        - compute_wake_effects!
        - compute_final_wind_shear!
    

        



Essentially, for runFLORIS replacements: 
- prepare_rotor_points!: no change
- handle_single_turbine!: no change
- setup_computation_buffers!: add veer related buffers to Floris struct and Params struct
- compute_wake_effects!: add veer related computations, probably in new file veer
- compute_final_wind_shear!: probably different way of combining wakes, as veer wake model already provides velocities









note to self: next time continue at checking centerline stuff,
    or just implement the wake model and see what variables can be transferred










somehow, the wakes are combined in a single number which is the effective wind speed.
 but, i think there should be a way to combine the velocity profile, also for power calculations
 this might require changes to the data structure


and also start writing code diagrams in the notebook!



wake superposition: directly via wind speeds



=#