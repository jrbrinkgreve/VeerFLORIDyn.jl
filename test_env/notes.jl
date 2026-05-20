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












      #54T



turbines:
    - id: 1 
      type: DTU 10MW 
      x: 2524 
      y: 5303 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 2 
      type: DTU 10MW 
      x: 5217 
      y: 5396 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 3 
      type: DTU 10MW 
      x: 6601 
      y: 5041 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 4 
      type: DTU 10MW 
      x: 7610 
      y: 5097 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 5 
      type: DTU 10MW 
      x: 8620 
      y: 5097 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 6 
      type: DTU 10MW 
      x: 9630 
      y: 5172 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 7 
      type: DTU 10MW 
      x: 10621 
      y: 5284 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 8 
      type: DTU 10MW 
      x: 11631 
      y: 5396 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 9 
      type: DTU 10MW 
      x: 804 
      y: 4704 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 10 
      type: DTU 10MW 
      x: 748 
      y: 3900 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 11 
      type: DTU 10MW 
      x: 729 
      y: 2647 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 12 
      type: DTU 10MW 
      x: 710 
      y: 1862 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 13 
      type: DTU 10MW 
      x: 654 
      y: 1058 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 14 
      type: DTU 10MW 
      x: 1739 
      y: 1189 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 15 
      type: DTU 10MW 
      x: 2767 
      y: 1245 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 16 
      type: DTU 10MW 
      x: 3758 
      y: 1638 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 17 
      type: DTU 10MW 
      x: 2599 
      y: 2049 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 18 
      type: DTU 10MW 
      x: 1645 
      y: 2311 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 19 
      type: DTU 10MW 
      x: 1645 
      y: 3601 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 20 
      type: DTU 10MW 
      x: 1664 
      y: 4630 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 21 
      type: DTU 10MW 
      x: 2599 
      y: 4349 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 22 
      type: DTU 10MW 
      x: 2524 
      y: 3302 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 23 
      type: DTU 10MW 
      x: 3459 
      y: 2984 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 24 
      type: DTU 10MW 
      x: 4338 
      y: 2647 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 25 
      type: DTU 10MW 
      x: 5254 
      y: 2404 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 26 
      type: DTU 10MW 
      x: 6114 
      y: 2086 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 27 
      type: DTU 10MW 
      x: 7274 
      y: 2143 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 28 
      type: DTU 10MW 
      x: 8284 
      y: 2273 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 29 
      type: DTU 10MW 
      x: 9106 
      y: 1843 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 30 
      type: DTU 10MW 
      x: 10172 
      y: 2068 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 31 
      type: DTU 10MW 
      x: 11220 
      y: 2311 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 32 
      type: DTU 10MW 
      x: 12285 
      y: 2573 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 33 
      type: DTU 10MW 
      x: 13277 
      y: 2797 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 34 
      type: DTU 10MW 
      x: 12603 
      y: 3639 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 35 
      type: DTU 10MW 
      x: 12622 
      y: 4592 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 36 
      type: DTU 10MW 
      x: 11668 
      y: 4387 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 37 
      type: DTU 10MW 
      x: 11594 
      y: 3395 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 38 
      type: DTU 10MW 
      x: 10602 
      y: 3152 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 39 
      type: DTU 10MW 
      x: 10602 
      y: 4143 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 40 
      type: DTU 10MW 
      x: 12659 
      y: 5639 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 41 
      type: DTU 10MW 
      x: 9630 
      y: 3395 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 42 
      type: DTU 10MW 
      x: 8676 
      y: 4349 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 43 
      type: DTU 10MW 
      x: 8415 
      y: 3208 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 44 
      type: DTU 10MW 
      x: 7442 
      y: 3171 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 45 
      type: DTU 10MW 
      x: 6395 
      y: 3078 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 46 
      type: DTU 10MW 
      x: 5310 
      y: 3395 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 47 
      type: DTU 10MW 
      x: 6152 
      y: 4069 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 48 
      type: DTU 10MW 
      x: 7199 
      y: 4087 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 49 
      type: DTU 10MW 
      x: 5236 
      y: 4368 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 50 
      type: DTU 10MW 
      x: 4357 
      y: 4667 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 51 
      type: DTU 10MW 
      x: 3440 
      y: 4985 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 52 
      type: DTU 10MW 
      x: 3459 
      y: 3994 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 53 
      type: DTU 10MW 
      x: 4357 
      y: 3676 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06
    - id: 54 
      type: DTU 10MW 
      x: 822 
      y: 5490 
      z: 0
      a: 0.33
      yaw: 0
      ti: 0.06




turbine_groups:
  - name: "all"
    id: 0
    turbines: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54]





================================================================================================================================================

if we just make a con.yaw constructor which fixed the past control input, and perturbs the
  wind.dir data, we can easily make an MPC-type controller which has the established flowfield data,
  predicts the best instruction for the future, and then applied the first step!


in pseudocode:


data = load_past_data
weather_prediction = load_weather_prediction
con.yaw_data = optimise(past_info, prediction)
apply_control(con.yaw_data)

repeat




























=#






