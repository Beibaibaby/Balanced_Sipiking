using Plots


# Define the functions
x = 0:0.01:10π
y_decaying = exp.(-0.1x) .* sin.(x)
y_non_decaying = sin.(x)

# Plot both functions with thicker lines
plot(x, y_decaying, label="Decaying Energy", xlabel="Time", ylabel="Phase", title="Sound Attenuation", linewidth=3)
plot!(x, y_non_decaying, label="With out Decaying Energy", linewidth=3)

# Save the plot to a file
savefig("sine_comparison_plot.png")
