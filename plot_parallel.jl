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
  
                corr_data = load_object(file_path)

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



if all([haskey(avg_cross_corr_E_E_event, 0), haskey(avg_cross_corr_C_C_event, 0), haskey(avg_cross_corr_I_I_event, 0),
    haskey(avg_cross_corr_E_E_nonevent, 0), haskey(avg_cross_corr_C_C_nonevent, 0), haskey(avg_cross_corr_I_I_nonevent, 0)])
    value_E_E_event = avg_cross_corr_E_E_event[0]
    value_C_C_event = avg_cross_corr_C_C_event[0]
    value_I_I_event = avg_cross_corr_I_I_event[0]

    value_E_E_nonevent = avg_cross_corr_E_E_nonevent[0]
    value_C_C_nonevent = avg_cross_corr_C_C_nonevent[0]
    value_I_I_nonevent = avg_cross_corr_I_I_nonevent[0]
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

end



function load_and_merge_data(main_dir::String)
    e_rates_all = Float64[]
    i_rates_all = Float64[]

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            # Check for the existence of cross_corr_E_E.jld2 in the subdirectory
            cross_corr_file_path = joinpath(sub_dir, "cross_corr_E_E.jld2")
            if isfile(cross_corr_file_path)
                e_file_path = joinpath(sub_dir, "e_rate_after_peak.jld2")
                i_file_path = joinpath(sub_dir, "i_rate_after_peak.jld2")

                if isfile(e_file_path) && isfile(i_file_path)
                    e_data = load(e_file_path, "e_rate_after_peak")
                    i_data = load(i_file_path, "i_rate_after_peak")

                    append!(e_rates_all, e_data)
                    append!(i_rates_all, i_data)
                    println("Processing subdirectory: ", sub_dir)
                    println("E data: ", e_data)
                end
            end
        end
    end

    return e_rates_all, i_rates_all
end



e_rates_merged, i_rates_merged = load_and_merge_data(main_dir)

@save joinpath(main_dir, "merged_e_rates.jld2") e_rates_merged
@save joinpath(main_dir, "merged_i_rates.jld2") i_rates_merged

function plot_rates_with_stats(e_rates, i_rates, output_file)
    e_mean = mean(e_rates)
    e_std = std(e_rates)
    i_mean = mean(i_rates)
    i_std = std(i_rates)

    plot_size = (600, 400)  # Adjust the size as needed
    left_margin = 10mm  # Increase left margin to ensure y-axis labels are visible
    mean_marker_size = 10  # Adjust the size of the marker for the mean

    p = scatter(ones(length(e_rates)), e_rates, label="E data", color=:pink, alpha=0.5, size=plot_size, left_margin=left_margin)
    scatter!(p, 2*ones(length(i_rates)), i_rates, label="I data", color=:lightblue, alpha=0.5)

    # Overlay the mean markers
    scatter!(p, [1], [e_mean], yerr=e_std, label="E mean", color=:red, markershape=:cross, markersize=mean_marker_size)
    scatter!(p, [2], [i_mean], yerr=i_std, label="I mean", color=:blue, markershape=:cross, markersize=mean_marker_size)
    
    xticks!(p, [1, 2], ["E", "I"])
    ylabel!(p, "Rate")
    xlabel!(p, "Neuron Type")
    savefig(p, output_file)
end



plot_rates_with_stats(e_rates_merged, i_rates_merged, joinpath(main_dir, "rates_plot.png"))
