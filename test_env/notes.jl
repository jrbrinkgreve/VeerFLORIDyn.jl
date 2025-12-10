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
- compute_final_wind_shear!: no change









note to self: next time continue at checking centerline stuff,
    or just implement the wake model and see what variables can be transferred







=#