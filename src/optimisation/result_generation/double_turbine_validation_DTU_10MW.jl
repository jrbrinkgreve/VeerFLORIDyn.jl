
#this script aims to validate the double turbine interaction:

#paper:
#https://arc.aiaa.org/doi/10.2514/6.2018-0755
using Statistics
using Plots
using StatsPlots


#run this section for result gathering




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
DTU 10 MW case:


noting the results:
no veer, aligned
2324.068099891971





no veer, -0.5D offset
3020.1369503441606






no veer, +0.5D offset
3020.1369503441606








-0.0448 deg/m veer, aligned
2542.5564954195224





-0.0448 deg/m veer, -0.5D offset
3070.3326394805354



-0.0448 deg/m veer, +0.5D offset
3031.649484766626



=#
















#plotting using this data:


plotting = true
if plotting

#=   0.062TI case:
# ── Data ──────────────────────────────────────────────────────────────────────
P_no_veer = [ 3020.1369503441606 , 2324.068099891971, 3020.1369503441606]  # -0.5D, 0D, +0.5D
P_veer    = [3031.649484766626, 2542.5564954195224,  3070.3326394805354]   # -0.5D, 0D, +0.5D
=#



#= 0.03TI case
P_no_veer = [2662.4996563033765, 1329.7786554051588, 2662.4996563033765]
P_veer = [2619.988313792356,1700.4866183828588, 2569.120063139294] 
=#



# 0.043 TI case
P_no_veer = [2835.7873453666125 ,1857.6424942284361,2835.7873453666125]
P_veer = [2812.071133213757, 2151.362451989349, 2856.738985978027 ]




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
 
savefig("output/double_turbine_validation_DTU_10MW.pdf")


end
nothing


