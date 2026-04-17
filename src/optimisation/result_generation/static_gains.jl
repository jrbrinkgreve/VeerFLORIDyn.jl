using CSV
using DataFrames
using Dates

# ── Configuration ─────────────────────────────────────────
const INPUT_CSV = "data/2021_9T_Data/WindDir.csv"
const OUTPUT_CSV = "data/2021_9T_Data/WindDir.csv"   # overwrites in-place; change if you want a copy
const LOG_FILE = "static_angle_sweep_gains.log"

# Wind directions to sweep over
wind_directions = collect(0:15:345)  # every 15 degrees, adjust as needed

# ── CSV helpers ───────────────────────────────────────────

function load_csv(path::String)::DataFrame
    CSV.read(path, DataFrame; header=["timestamp", "wind_direction"], skipto=1)
end

function set_wind_direction(df::DataFrame, new_dir::Number)::DataFrame
    df[!, :wind_direction] .= new_dir
    return df
end

function save_csv(df::DataFrame, path::String)
    CSV.write(path, df; header=false)
end


# ── Logging ───────────────────────────────────────────────

function log_result(wind_dir::Number, result::Dict)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    entries = join(["$k=$(v)" for (k, v) in result], ", ")
    line = "$timestamp | wind_dir=$(wind_dir)° | $entries"
    println(line)
    open(LOG_FILE, "a") do f
        println(f, line)
    end
end

# ── Main loop ─────────────────────────────────────────────


#baseline_power_avgs  = Float64[]
#optimized_power_avgs = Float64[]
#energy_increases     = Float64[]




println("Starting sweep over $(length(wind_directions)) wind directions...\n")

for wind_dir in wind_directions
    println("── Wind direction: $(wind_dir)° ──")

    # Edit WindDir.csv
    df = load_csv(INPUT_CSV)
    df = set_wind_direction(df, wind_dir)
    save_csv(df, OUTPUT_CSV)

    # Run simulation
    #result = run_simulation()
    include("../parallel_yaw_optimisation.jl")
   
    # Log result
    log_result(wind_dir, Dict(
        "baseline_power_avg" => baseline_power_avg,
        "optimized_power_avg" => optimized_power_avg,
        "energy_increase_over_baseline" => energy_increase_over_baseline,
        "l1_baseline_yaw_avg" => l1_baseline_yaw_avg,
        "l1_optimized_yaw_avg" => l1_optimized_yaw_avg,
        "feasibility_baseline" => feasibility_baseline,
        "feasibility_optimized" => feasibility_optimized,
        "minimizer" => result.minimizer,
    ))

    push!(baseline_power_avgs,  baseline_power_avg)
    push!(optimized_power_avgs, optimized_power_avg)
    push!(energy_increases,     energy_increase_over_baseline)
end




plot(wind_directions, energy_increases, marker=:o, xlabel="Wind Direction (°)", ylabel="Energy Increase over Baseline (%)", title="Energy Increase vs Wind Direction", legend=false)


