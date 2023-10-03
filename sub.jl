using Plots

# Define epochs and initialize arrays for the ranks
epochs = 1:2000
approx_rank = zeros(length(epochs))
actual_rank = zeros(length(epochs))

# Define parameters
initial_approx_rank = 300
final_approx_rank = 10
half_life_approx = 80  # This determines how quickly the approx rank decreases initially
noise_level_approx = 12

initial_actual_rank = 350
final_actual_rank = 240
half_life_actual = 200
noise_level_actual = 4

# Generate approximate rank: Exponential decrease from 300 to 70 by epoch 500 with noise, then remains around 70
decay_rate_approx = log(2) / half_life_approx
approx_rank[1:500] .= initial_approx_rank .* 
                     exp.(-decay_rate_approx .* epochs[1:500]) .+ 
                     noise_level_approx .* randn(500)
approx_rank[501:end] .= final_approx_rank .+ 2.2 .* randn(length(epochs[501:end]))

# Generate actual rank: Exponential decrease from 300 to 250 by epoch 120, then stays at 250
decay_rate_actual = log(2) / half_life_actual
actual_rank[1:120] .= initial_actual_rank .* 
                     exp.(-decay_rate_actual .* epochs[1:120]) .+ 
                     noise_level_actual .* randn(120)
actual_rank[121:end] .= final_actual_rank .+ 1 .* randn(length(epochs[121:end]))

# Create the plot
plot(epochs, [approx_rank, actual_rank],
    label = ["Approximate Rank" "Actual Rank"],
    color = [:blue :red],
    xlabel = "Epochs",
    ylabel = "Rank",
    title = "Approximate vs Actual Rank over Epochs",
    lw = 2,  # Line width
    legend = :bottomright)

# Save the plot as a PNG file
savefig("rank_plot.png")
