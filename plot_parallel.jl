using JLD2
using FilePathsBase
using Plots
include("sim_6.jl") 
using JLD2
using Measures

function average_correlations(main_dir::String, file_name::String)
    total_corr = Dict()
    count = 0

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)

            if isfile(file_path)
                #println("Loading file: $file_path")  # Debug statement
                corr_data = load_object(file_path)
                #println("Loaded data: ", corr_data)  # Debug statement
                # Aggregate data from each file
                for (key, value) in corr_data
                    total_corr[key] = get(total_corr, key, 0) + value
                end
                count += 1
            end
        end
    end

    # Average the aggregated data
    if count > 0
        for key in keys(total_corr)
            total_corr[key] /= count
        end
    end

    return total_corr
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
event_thre = parse(Float64, get_arg("--event_thre", "2.0"))

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
plot_correlations(avg_cross_corr_E_E, avg_cross_corr_I_I, avg_cross_corr_E_I, avg_cross_corr_I_E, avg_cross_corr_C_C, joinpath(main_dir, "plot_corr_overall.png"),event_thre)

# Plot the averaged correlations for event-based
plot_correlations(avg_cross_corr_E_E_event, avg_cross_corr_I_I_event, avg_cross_corr_E_I_event, avg_cross_corr_I_E_event, avg_cross_corr_C_C_event, joinpath(main_dir, "plot_corr_events.png"),event_thre)

# Plot the averaged correlations for nonevent-based
plot_correlations(avg_cross_corr_E_E_nonevent, avg_cross_corr_I_I_nonevent, avg_cross_corr_E_I_nonevent, avg_cross_corr_I_E_nonevent, avg_cross_corr_C_C_nonevent, joinpath(main_dir, "plot_corr_nonevents.png"),event_thre)

value_E_E_event = avg_cross_corr_E_E_event[0]
value_C_C_event = avg_cross_corr_C_C_event[0]
value_I_I_event = avg_cross_corr_I_I_event[0]

value_E_E_nonevent = avg_cross_corr_E_E_nonevent[0]
value_C_C_nonevent = avg_cross_corr_C_C_nonevent[0]
value_I_I_nonevent = avg_cross_corr_I_I_nonevent[0]


# Prepare data for plotting
categories = ["E_E", "C_C", "I_I"]
event_data = [value_E_E_event, value_C_C_event, value_I_I_event]
nonevent_data = [value_E_E_nonevent, value_C_C_nonevent, value_I_I_nonevent]

# Create the plot
p = plot(xticks = (1:3, categories))
scatter!(p, 1:3, event_data, label = "Event", color = :blue)
scatter!(p, 1:3, nonevent_data, label = "Non-event", color = :red)
plot!(p, 1:3, event_data, label = "", color = :blue, line = :line)
plot!(p, 1:3, nonevent_data, label = "", color = :red, line = :line)

# Save the plot to a file
savefig(p, joinpath(main_dir, "event_vs_nonevent.png"))