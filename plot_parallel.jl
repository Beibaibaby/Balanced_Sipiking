using JLD2
using FilePathsBase
using Plots
include("sim_5.jl") 
# Function to load and average correlation data
function average_correlations(main_dir::String, file_name::String)
    total_corr = nothing
    count = 0

    # Iterate over each subdirectory in the main directory
    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)
            if isfile(file_path)
                @load file_path corr_data
                if total_corr === nothing
                    total_corr = corr_data
                else
                    total_corr .+= corr_data
                end
                count += 1
            end
        end
    end

    return count > 0 ? total_corr / count : total_corr
end

function get_arg(key, default)
    index = findfirst(==(key), ARGS)
    if index !== nothing && index < length(ARGS)
        return ARGS[index + 1]
    end
    return default
end

# Retrieve the main directory from command-line arguments
main_dir = get_arg("--main_dir", "/default/path/to/main_dir")

# Average each type of correlation
avg_cross_corr_E_E = average_correlations(main_dir, "cross_corr_E_E.jld2")
avg_cross_corr_I_I = average_correlations(main_dir, "cross_corr_I_I.jld2")
avg_cross_corr_E_I = average_correlations(main_dir, "cross_corr_E_I.jld2")
avg_cross_corr_I_E = average_correlations(main_dir, "cross_corr_I_E.jld2")
avg_cross_corr_C_C = average_correlations(main_dir, "cross_corr_C_C.jld2")

# Average each type of event-based correlation
avg_cross_corr_E_E_event = average_correlations(main_dir, "cross_corr_E_E_event.jld2")
avg_cross_corr_I_I_event = average_correlations(main_dir, "cross_corr_I_I_event.jld2")
avg_cross_corr_E_I_event = average_correlations(main_dir, "cross_corr_E_I_event.jld2")
avg_cross_corr_I_E_event = average_correlations(main_dir, "cross_corr_I_E_event.jld2")
avg_cross_corr_C_C_event = average_correlations(main_dir, "cross_corr_C_C_event.jld2")

# Average each type of nonevent-based correlation
avg_cross_corr_E_E_nonevent = average_correlations(main_dir, "cross_corr_E_E_nonevent.jld2")
avg_cross_corr_I_I_nonevent = average_correlations(main_dir, "cross_corr_I_I_nonevent.jld2")
avg_cross_corr_E_I_nonevent = average_correlations(main_dir, "cross_corr_E_I_nonevent.jld2")
avg_cross_corr_I_E_nonevent = average_correlations(main_dir, "cross_corr_I_E_nonevent.jld2")
avg_cross_corr_C_C_nonevent = average_correlations(main_dir, "cross_corr_C_C_nonevent.jld2")

# Plot the averaged correlations for overall
plot_correlations(avg_cross_corr_E_E, avg_cross_corr_I_I, avg_cross_corr_E_I, avg_cross_corr_I_E, avg_cross_corr_C_C, joinpath(main_dir, "plot_corr_overall.png"))

# Plot the averaged correlations for event-based
plot_correlations(avg_cross_corr_E_E_event, avg_cross_corr_I_I_event, avg_cross_corr_E_I_event, avg_cross_corr_I_E_event, avg_cross_corr_C_C_event, joinpath(main_dir, "plot_corr_events.png"))

# Plot the averaged correlations for nonevent-based
plot_correlations(avg_cross_corr_E_E_nonevent, avg_cross_corr_I_I_nonevent, avg_cross_corr_E_I_nonevent, avg_cross_corr_I_E_nonevent, avg_cross_corr_C_C_nonevent, joinpath(main_dir, "plot_corr_nonevents.png"))

