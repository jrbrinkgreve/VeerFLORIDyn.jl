#=
notes for veer implementations


#call structure for setUpTmpWFAndRun!

...

for all turbines
    if single turbine
    - runFLORIS()
    


    runFLORIS:
        - prepare_rotor_points!
        - handle_single_turbine!
        - setup_computation_buffers!
        - compute_wake_effects!
        - compute_final_wind_shear!
    









note to self: next time continue at checking centerline stuff,
    or just implement the wake model and see what variables can be transferred







=#