#=
notes for veer implementations




- Check centerline stuff!

note to self: next time continue at checking centerline stuff,
    or just implement the wake model and see what variables can be transferred







- More advanced power calculations

somehow, the wakes are combined in a single number which is the effective wind speed.
 but, i think there should be a way to combine the velocity profile, also for power calculations
 this might require changes to the data structure


and also start writing code diagrams in the notebook!




- TSR vs CP

Following Tamaro et al. (2024a) (Sect. 3.2), the power coefficient is calculated as
CP(λ, θ, γ ) = CP(λ, θ, 0) ηP(λ, θ, γ )       
!!!! for TSR vs CP stuff - from Tamaro, A robust active power control

https://wes.copernicus.org/articles/10/2705/2025/


- Is wake steering working with the new veer model? test this,
    also comparing to the gaussian model







- overhaul veer modelling:
    sim.floris    struct field
    also: Cp vs tsr overhaul

- Report
    check cos^2(...) instead of cos(..)^2 in report
    list of parameters in methods section


- also check sign in code of veer implementation
- convention of yaw in gaussian: power calculation code whether 





- Optimisation crashing with large yaw angles: 
    - Gaussian wake model centerline! not accepting large yaws, will crash


=#




