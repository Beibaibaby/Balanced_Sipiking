using Statistics
using Plots  

doplot = false

include("sim_network.jl")
K = 300
Ne = 500 #Ne is the Number of E cells
Ni = 500 #Ni is the Number of I cells
T =  1500 #T is the simulation time (ms)
sqrtK = sqrt(K) 
taue = 10 #membrane time constant for exc. neurons (ms)
taui = 10 #membrane time constant for inh. neurons (ms)

jie = 2.95 / (taui*sqrtK)        #strength from E to I
jei = -3.6 /(taue*sqrtK)      #strength from I to E
jii = -2.9 / (taui*sqrtK)      #strength from I to I 
jee = 3.5 / (taue*sqrtK)        #strength from E to E 

times, ns, Ne, Ncells, T = sim(K,Ne,Ni,T,jie,jei,jii,jee)

println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne+1):Ncells] / T), " Hz")

function compute_sliding_rate(spiketimes, window_size, step_size, T)
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    for i = 1:n_steps
        t_start = (i-1) * step_size + 1  # Ensure the start time is non-zero
        t_end = t_start + window_size - 1  # Adjust end time based on start
        
        # Check for out-of-bounds
        if t_end > T
            break
        end

        n_spikes = sum([sum((t_start .<= spiketimes[j, :]) .& (spiketimes[j, :] .< t_end)) for j=1:n_neurons])
        rates[i] = n_spikes / (n_neurons * window_size) * 1000  # rate in Hz
    end
    #println(rates)
    return rates
end


if doplot 
    println("creating plot")

    # Define plot size: (width, height)
    plot_size = (600, 400) 

    p = scatter(legend=false, xlabel="Time (ms)", ylabel="Neuron", ylims=(0, Ne), xlims=(0, T),
                size=plot_size)
    for ci = 1:Ne
        vals = times[ci, 1:ns[ci]]
        y = ci .* ones(length(vals))
        scatter!(vals, y, markersize=1, markercolor=:black, markerstrokecolor=:black, legend=false, show=false)
    end
    savefig(p, "raster.png")

    # ...

# Parameters for sliding window
# Parameters for sliding window
window_size = 100  # in ms
step_size = 5     # in ms

e_rate = compute_sliding_rate(times[1:Ne, :], window_size, step_size, T)
i_rate = compute_sliding_rate(times[(Ne+1):Ncells, :], window_size, step_size, T)

# Compute the time values based on window_size and step_size
n_steps = length(e_rate)  # or length(i_rate), assuming they have the same length
time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]

p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)",
          label="Excitatory", lw=2, linecolor=:red, size=plot_size)
plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:skyblue)


savefig(p2, "firing_rates_testing.png")

# ...

end

