using Dates
using Plots
# Define a struct to hold the parsed data neatly
mutable struct data_record
    timestamp::DateTime
    wind_dir_deg::Float64
    l1_optimized_yaw_avg::Float64
    l1_baseline_yaw_avg::Float64
    energy_increase::Float64
    feasibility_optimized::Float64
    baseline_power_avg::Float64
    feasibility_baseline::Float64
    minimizer::Vector{Float64}
    optimized_power_avg::Float64
end

function parse_wind_log(file_path::String)
    records = data_record[]
    
    open(file_path, "r") do file
        for line in eachline(file)
            isempty(strip(line)) && continue
            
            # Split by the pipe operator
            parts = split(line, "|")
            
            # 1. Parse Timestamp
            ts = DateTime(strip(parts[1]), "yyyy-mm-dd HH:MM:SS")
            
            # 2. Parse Wind Dir (remove the ° symbol)
            dir_str = strip(parts[2])
            wind_dir = parse(Float64, replace(split(dir_str, "=")[2], "°" => ""))
            
            # 3. Parse the remaining KV pairs
            # We use a regex to handle the 'minimizer=[...]' part specifically
            kv_section = strip(parts[3])
            
            # Extract minimizer list first
            m_match = match(r"minimizer=\[(.*?)\]", kv_section)
            minimizer_vals = Float64[]
            if m_match !== nothing
                minimizer_vals = parse.(Float64, split(m_match.captures[1], ","))
            end
            
            # Clean KV section to parse other numeric values
            # (Removing the minimizer part to avoid confusion during splitting)
            clean_kv = replace(kv_section, r"minimizer=\[.*?\]" => "")
            kv_pairs = split(clean_kv, ",")
            
            data_dict = Dict{String, Float64}()
            for pair in kv_pairs
                s_pair = strip(pair)
                if contains(s_pair, "=")
                    k, v = split(s_pair, "=")
                    data_dict[strip(k)] = parse(Float64, v)
                end
            end
            
            # Push to results
            push!(records, data_record(
                ts, wind_dir,
                data_dict["l1_optimized_yaw_avg"],
                data_dict["l1_baseline_yaw_avg"],
                data_dict["energy_increase_over_baseline"],
                data_dict["feasibility_optimized"],
                data_dict["baseline_power_avg"],
                data_dict["feasibility_baseline"],
                minimizer_vals,
                data_dict["optimized_power_avg"]
            ))
        end
    end
    return records
end

# Usage example:
# data = parse_wind_log("your_data_file.txt")
# println("Loaded $(length(data)) rows.")
# println("First row energy increase: $(data[1].energy_increase)%")

data_6 = parse_wind_log("output/6.0_ms_static_angle_sweep_gains.log");
data_8 = parse_wind_log("output/8.2_ms_static_angle_sweep_gains.log");
data_10 = parse_wind_log("output/10.0_ms_static_angle_sweep_gains.log");

using DataFrames

# Convert your array of structs into a DataFrame
df6 = DataFrame(data_6)
df8 = DataFrame(data_8)
df10 = DataFrame(data_10)


wind_directions = collect(0:15:345)  # every 15 degrees, adjust as needed
# First plot: initialize the canvas and labels
plot(wind_directions, df6.energy_increase, 
    marker=:o, 
    label="6.0 m/s", 
    xlabel="Wind Direction (°)", 
    ylabel="Energy Increase over Baseline (%)", 
    title="Energy Increase vs. Wind Direction", 
    legend=:topright)

# Subsequent plots: use plot! to add to the existing canvas
plot!(wind_directions, df8.energy_increase, marker=:o, label="8.2 m/s")
plot!(wind_directions, df10.energy_increase, marker=:o, label="10.0 m/s")

savefig("output/energy_increase_vs_wind_direction.pdf")