#calculating the veer magitudes for the tests to be conducted:
#validation is based on the paper: 
#https://www.sciencedirect.com/science/article/pii/S0960148126006804?via%3Dihub



using Plots
using LaTeXStrings

# Physical parameters (adjust as needed for your turbine)
p_p = 2.2 # FLORIS exponent
D = 178.4 # Example diameter
R = D / 2

cwp00_veer_gradient = 0.0 / D
vwp15_veer_gradient = 15.0 / D
vwp25_veer_gradient = 25.0 / D
vwp35_veer_gradient = 35.0 / D




# Discrete yaw angles to match the reference plot
yaw_discrete = [-30, -20, -10, 0, 10, 20, 30]

# Define veer gradients
# 0.0 case is excluded as it is the 0% reference baseline
gradients = Dict(
"VWP15" => 15.0 / D,
"VWP25" => 25.0 / D,
"VWP35" => 35.0 / D
)

# Plot configuration with LaTeX labels and specific Y-axis range
p = plot(
xlabel=L"Yaw angle $\gamma$ (deg)",
ylabel=L"$\eta_{\bar{P}}$ (%)",
xticks=(yaw_discrete, string.(yaw_discrete)),
yticks=-12:2:4,
ylims=(-12.5, 4.5),
legend=:bottomright,
size=(600, 350),
leftmargin=5Plots.mm,
rightmargin=5Plots.mm,
bottommargin=7Plots.mm,
grid=false,
legendfontsize=10,
gridfontsize=12,
tickfontsize=12,
guidefontsize = 12,
)

# Add the 0% reference line (dashed grey)
hline!(p, [0], color=:grey, linestyle=:dash, label="")

# Define specific styles to match the reference image
plot_config = Dict(
"VWP15" => (color=:red, marker=:circle),
"VWP25" => (color=:teal, marker=:dtriangle),
"VWP35" => (color=:dodgerblue, marker=:utriangle)
)

# Iterate through each veer case, calculate theoretical factor, and plot
# Use sort to ensure VWP15 is on top in the legend
for key in sort(collect(keys(gradients)), rev=true)
local alpha = gradients[key]
local alpha_rad = deg2rad(alpha)

# Calculate theoretical correction factor from formula
cf = 1 .-  ((p_p * alpha_rad^2 * R^2) / 8) .* (1 .- (p_p - 1) .* tand.(yaw_discrete).^2)

# Convert to percentage change relative to baseline (1.0)
eta_p_percent = (cf .- 1.0) .* 100

config = plot_config[key]

# Plot with discrete markers and specific colors
plot!(p, yaw_discrete, eta_p_percent,
label=key,
color=config[:color],
marker=config[:marker],
markerstrokecolor=config[:color],
markerfacecolor=:white, # Open markers like reference
markersize=8,
linewidth=2,

)
end



display(p)
savefig("output/ma_images.pdf-3.pdf")