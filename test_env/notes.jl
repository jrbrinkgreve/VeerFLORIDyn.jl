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






#small wf

turbines:
    - id: 1
      type: DTU 10MW
      x: 1347.6
      y: 919.0
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 2
      type: DTU 10MW
      x: 1640.9
      y: 1662.5
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 3
      type: DTU 10MW
      x: 2248.0
      y: 1001.0
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 4
      type: DTU 10MW
      x: 1934.2
      y: 2406.0
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 5
      type: DTU 10MW
      x: 3149.3
      y: 1083.9
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 6
      type: DTU 10MW
      x: 2227.5
      y: 3149.5
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 7
      type: DTU 10MW
      x: 4218.2
      y: 1564.4
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 8
      type: DTU 10MW
      x: 2520.8
      y: 3893.0
      z: 0
      a: 0.33
      yaw: 30
      ti: 0.06
    - id: 9
      type: DTU 10MW
      x: 4097.0
      y: 2666.3
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 10
      type: DTU 10MW
      x: 4036.4
      y: 3711.2
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06








=#




