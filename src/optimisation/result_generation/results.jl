


include("helper_functions.jl")
#code for averaged result generation
array = [1e3, 2e3, 3e3, 5e3, 7e3, 1e4, 2e4, 3e4, 5e4]
for set_lambda_l1 in array
    @eval set_lambda_l1 = $set_lambda_l1
    include("../parallel_yaw_optimisation.jl")
end


#=
YAW LIMITS

=======================================================
Inf deg limit


Total optimization time: 322.63 seconds
Maximum yaw misalignment: 28.233 degrees, at time step CartesianIndex(1309, 6)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3985.68 kW
Increase over baseline: 1.77 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 54.79 deg



=======================================================
25deg limit

Total optimization time: 240.63 seconds
Maximum yaw misalignment: 24.991 degrees, at time step CartesianIndex(1264, 6)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3983.84 kW
Increase over baseline: 1.72 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 51.44 deg







=======================================================
20deg limit


Total optimization time: 208.44 seconds
Maximum yaw misalignment: 19.99 degrees, at time step CartesianIndex(1211, 1)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3974.89 kW
Increase over baseline: 1.49 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 44.09 deg








=======================================================
15deg limit

Total optimization time: 178.79 seconds
Maximum yaw misalignment: 14.979 degrees, at time step CartesianIndex(1215, 4)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3959.8 kW
Increase over baseline: 1.11 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 36.92 deg






=======================================================
10deg limit

Total optimization time: 132.33 seconds
Maximum yaw misalignment: 9.998 degrees, at time step CartesianIndex(1037, 1)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3937.83 kW
Increase over baseline: 0.55 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 31.27 deg


























================================================================================================
SPARSITY LIMITS
NOTE: COMPARE BASELINE TO L0==4 VALUE

L0==4 baseline power:  





L0==1

Total optimization time: 160.14 seconds
Maximum yaw misalignment: 23.473 degrees, at time step CartesianIndex(1, 4)

Baseline average power per turbine:  3780.11 kW
Optimized average power per turbine: 3798.91 kW
Increase over baseline: -3.00 %

Baseline average L1 yaw change norm: 0.0 deg
Optimized average L1 yaw change norm: 0.0 deg




L0==2

Total optimization time: 247.3 seconds
Maximum yaw misalignment: 15.446 degrees, at time step CartesianIndex(2855, 1)

Baseline average power per turbine:  3864.92 kW
Optimized average power per turbine: 3910.74 kW
Increase over baseline: -0.16 %

Baseline average L1 yaw change norm: 22.78 deg
Optimized average L1 yaw change norm: 20.45 deg








L0==3

Total optimization time: 249.49 seconds
Maximum yaw misalignment: 19.588 degrees, at time step CartesianIndex(1203, 4)

Baseline average power per turbine:  3906.92 kW
Optimized average power per turbine: 3938.71 kW
Increase over baseline: 0.57 %

Baseline average L1 yaw change norm: 30.37 deg
Optimized average L1 yaw change norm: 28.33 deg




L0==4
-







L0==5

Total optimization time: 237.97 seconds
Maximum yaw misalignment: 24.989 degrees, at time step CartesianIndex(1367, 2)

Baseline average power per turbine:  3919.71 kW
Optimized average power per turbine: 3989.31 kW
Increase over baseline: 1.86 %

Baseline average L1 yaw change norm: 35.44 deg
Optimized average L1 yaw change norm: 55.87 deg






L0==6

Total optimization time: 222.48 seconds
Maximum yaw misalignment: 24.985 degrees, at time step CartesianIndex(1322, 6)

Baseline average power per turbine:  3922.78 kW
Optimized average power per turbine: 3996.65 kW
Increase over baseline: 2.05 %

Baseline average L1 yaw change norm: 36.2 deg
Optimized average L1 yaw change norm: 61.11 deg















================================================================================================
REGULARISATION PARAMETER LAMBDA L1

1e3

Total optimization time: 274.34 seconds
Maximum yaw misalignment: 24.976 degrees, at time step CartesianIndex(1244, 4)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3983.07 kW
Increase over baseline: 1.7 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 49.95 deg








2e3
Total optimization time: 269.2 seconds
Maximum yaw misalignment: 24.904 degrees, at time step CartesianIndex(1259, 6)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3983.23 kW
Increase over baseline: 1.71 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 48.31 deg





3e3

Total optimization time: 270.36 seconds
Maximum yaw misalignment: 24.68 degrees, at time step CartesianIndex(1260, 2)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3980.87 kW
Increase over baseline: 1.65 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 44.69 deg


5e3

Total optimization time: 300.64 seconds
Maximum yaw misalignment: 23.271 degrees, at time step CartesianIndex(1237, 2)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3975.08 kW
Increase over baseline: 1.5 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 39.25 deg



7e3

Total optimization time: 292.26 seconds
Maximum yaw misalignment: 21.586 degrees, at time step CartesianIndex(1231, 2)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3965.61 kW
Increase over baseline: 1.26 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 32.78 deg





1e4

Total optimization time: 318.95 seconds
Maximum yaw misalignment: 20.055 degrees, at time step CartesianIndex(1209, 4)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3949.66 kW
Increase over baseline: 0.85 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 24.77 deg



2e4

Total optimization time: 314.93 seconds
Maximum yaw misalignment: 18.297 degrees, at time step CartesianIndex(1213, 4)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3916.45 kW
Increase over baseline: 0.0 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 14.72 deg






3e4

Total optimization time: 308.0 seconds
Maximum yaw misalignment: 21.257 degrees, at time step CartesianIndex(1, 6)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3886.91 kW
Increase over baseline: -0.75 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 9.83 deg






5e4

Total optimization time: 301.23 seconds
Maximum yaw misalignment: 22.186 degrees, at time step CartesianIndex(1, 4)

Baseline average power per turbine:  3916.4 kW
Optimized average power per turbine: 3848.13 kW
Increase over baseline: -1.74 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 5.14 deg

















================================================================================================
veer gradient:
0.01deg/m

Total optimization time: 277.91 seconds
Maximum yaw misalignment: 24.95 degrees, at time step CartesianIndex(3000, 1)

Baseline average power per turbine:  3914.99 kW
Optimized average power per turbine: 3981.16 kW
Increase over baseline: 1.69 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 51.46 deg






0.03 deg/m

Total optimization time: 286.1 seconds
Maximum yaw misalignment: 24.92 degrees, at time step CartesianIndex(1207, 4)

Baseline average power per turbine:  3904.37 kW
Optimized average power per turbine: 3966.13 kW
Increase over baseline: 1.58 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 51.55 deg




0.05 deg/m
Total optimization time: 294.83 seconds
Maximum yaw misalignment: 24.926 degrees, at time step CartesianIndex(1195, 4)

Baseline average power per turbine:  3884.84 kW
Optimized average power per turbine: 3944.19 kW
Increase over baseline: 1.53 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 50.27 deg




0.07deg/m

Total optimization time: 303.18 seconds
Maximum yaw misalignment: 24.985 degrees, at time step CartesianIndex(1148, 2)

Baseline average power per turbine:  3859.09 kW
Optimized average power per turbine: 3914.75 kW
Increase over baseline: 1.44 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 49.21 deg






0.1 deg/m
Total optimization time: 300.02 seconds
Maximum yaw misalignment: 24.868 degrees, at time step CartesianIndex(1096, 2)

Baseline average power per turbine:  3814.67 kW
Optimized average power per turbine: 3866.63 kW
Increase over baseline: 1.36 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 46.51 deg







0.15 deg/m
Total optimization time: 331.28 seconds
Maximum yaw misalignment: 24.817 degrees, at time step CartesianIndex(1068, 2)

Baseline average power per turbine:  3740.25 kW
Optimized average power per turbine: 3790.86 kW
Increase over baseline: 1.35 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 41.85 deg








0.2 deg/m
Total optimization time: 333.24 seconds
Maximum yaw misalignment: 23.871 degrees, at time step CartesianIndex(1037, 2)

Baseline average power per turbine:  3678.85 kW
Optimized average power per turbine: 3728.18 kW
Increase over baseline: 1.34 %

Baseline average L1 yaw change norm: 34.17 deg
Optimized average L1 yaw change norm: 39.49 deg








=#
nothing
