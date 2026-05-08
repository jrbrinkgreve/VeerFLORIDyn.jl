#https://arc.aiaa.org/doi/epdf/10.2514/6.2018-0755



using Statistics
using Plots
using StatsPlots







run = false

if run
    include("../../../examples/veer_main_mini.jl")

    num_timesteps_to_skip = 125 #skip the first 125 timesteps to avoid startup effects, which is 25 seconds at 5Hz

    powervector = md.PowerGen[1+(num_timesteps_to_skip+1)*wf.nT:end] * 1000.0

    t1_powers = mean(powervector[1:wf.nT:end] )
    println("Turbine 1 average power: $t1_powers kW")

    t2_powers = mean(powervector[2:wf.nT:end] )
    println("Turbine 2 average power: $t2_powers kW")
end

#=
NREL 5 MW case:


noting the results:
no veer, aligned
1162.3150380264303





no veer, -0.5D offset
1510.4336204466147




no veer, +0.5D offset
1510.4336204466151






-0.06349206349206349 deg/m veer, aligned
1263.8728416886713




-0.06349206349206349 deg/m veer, -0.5D offset
1530.8998239500309




-0.06349206349206349 deg/m veer, +0.5D offset
1516.5150035524189



=#




# ── Data 0.062TI──────────────────────────────────────────────────────────────────────
P_no_veer = [1510.4336204466151, 1162.3150380264303,  1510.4336204466147]  # -0.5D, 0D, +0.5D
P_veer    = [1516.5150035524189, 1263.8728416886713, 1530.8998239500309]  # -0.5D, 0D, +0.5D


P_no_veer 
P_veer

P_ref = P_no_veer[2]  # no-veer, aligned (0D)
 
y_no_veer = P_no_veer ./ P_ref
y_veer    = P_veer    ./ P_ref




# ── Plot ──────────────────────────────────────────────────────────────────────
data = hcat(y_no_veer, y_veer)   # nx2 matrix: each row = one group
 
groupedbar(data;
    bar_position = :dodge,
    bar_width    = 0.6,
    color        = [colorant"#4db8b8" colorant"#f0a830"],
    label        = ["no veer" "veer"],
    xticks       = (1:3, ["-0.5D", "0D", "+0.5D"]),
    yticks       = 0:0.25:2.0,
    ylims        = (0, 1.9),
    ylabel       = "Normalized Mean Power",
    legend       = :topleft,
    framestyle   = :box,
    grid         = false,
    size         = (600, 420),
    dpi          = 300,
)
 
savefig("output/double_turbine_validation_NREL_5MW.pdf")





nothing


