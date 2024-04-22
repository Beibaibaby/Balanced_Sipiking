using FileIO
using JLD2
using Plots
using StatsBase, Statistics

ENV["GKSwstype"] = "100"

# Set the directory name where files are stored
dir_name = "/gpfs/data/doiron-lab/draco/”"

# Loading the files
@load joinpath(dir_name, "E_input.jld2") E_input
@load joinpath(dir_name, "I_input.jld2") I_input





# Step 2A: Define the indices for the time window of interest
time_indices = 5500:10500  # corresponds to 100 ms to 200 ms

using StatsBase, Statistics, Plots

# Indices to skip
skip_indices = Set([1155, 1236, 995, 1417, 140, 208, 512, 775, 883, 969, 1122, 1398, 1950, 2086, 2545, 3552, 3998, 21, 98, 117, 302, 389, 392, 430, 494, 547, 643, 681, 808, 893, 965, 1006, 1184, 1187, 1229, 1401, 1647, 1715, 1735, 1933, 2176, 2183, 2292, 2298, 2379, 2392, 2438, 2465, 2509, 2609, 2672, 2698, 2724, 2891, 2892, 3074, 3092, 3112, 3222, 3310, 3376, 3407, 3429, 3435, 3497, 3910, 12, 69, 100, 114, 136, 144, 162, 223, 238, 303, 349, 358, 369, 385, 411, 426, 445, 451, 504, 563, 580, 595, 646, 711, 776, 786, 800, 849, 876, 879, 897, 940, 1056, 1065, 1081, 1088, 1099, 1135, 1180, 1200, 1201, 1205, 1262, 1300, 1322, 1349, 1443, 1488, 1489, 1491, 1587, 1624, 1691, 1713, 1732, 1749, 1753, 1755, 1769, 1852, 1868, 1935, 1973, 1997, 2008, 2018, 2047, 2077, 2147, 2167, 2193, 2196, 2220, 2222, 2252, 2263, 2301, 2420, 2440, 2488, 2503, 2513, 2588, 2708, 2711, 2740, 2747, 2811, 2843, 2853, 2859, 2886, 3049, 3086, 3167, 3169, 3186, 3208, 3260, 3265, 3271, 3345, 3349, 3350, 3483, 3484, 3502, 3546, 3559, 3572, 3580, 3625, 3631, 3644, 3685, 3689, 3719, 3732, 3754, 3756, 3757, 3764, 3769, 3770, 3841, 3866, 3916, 3922, 3948, 3965, 1, 10, 13, 42])

# Adjusted neuron indices for E and I neurons, excluding skipped ones
E_neurons = setdiff(1:4000, skip_indices)
I_neurons = setdiff(4001:5000, skip_indices)

# Step 2C: Function to calculate average Pearson correlation for random pairs
function average_correlation(input_data, neuron_indices, num_pairs=10000)
    correlations = Float64[]
    for _ = 1:num_pairs
        pair = sample(neuron_indices, 2, replace=false)
        cor = Statistics.cor(input_data[pair[1], time_indices], input_data[pair[2], time_indices])
        push!(correlations, cor)
    end
    mean(correlations)
end

# Step 2D: Compute E-E, I-I, and Total-Total correlations
E_E_correlation = average_correlation(E_input, E_neurons)
I_I_correlation = average_correlation(I_input, E_neurons)
Total_Total_correlation = average_correlation(E_input + I_input,E_neurons)  # assuming broadcasting addition over matrices

# Step 2E: Plot the results
correlation_values = [E_E_correlation, Total_Total_correlation, I_I_correlation]
labels = ["E-E", "Total-Total", "I-I"]
plot(labels, correlation_values, seriestype=:scatter, title="Neuron-Neuron Correlation", xlabel="Pair Type", ylabel="Average Pearson Correlation", legend=false, markersize=8)

# Save the plot
savefig("Neuron_Correlations.png")
