using Statistics
using Plots  
using Dates  # for generating timestamps
using JSON3  # Use JSON3 package for JSON handling
using Random


struct NetworkParameters
    K::Int
    Ne::Int
    Ni::Int
    T::Int
    sqrtK::Float64
    taue::Int
    taui::Int
    jie::Float64
    jee::Float64
    jei::Float64
    jii::Float64
end

doplot = true

include("network.jl")
include("sim.jl")
K = 800
Ne = 4000 #Ne is the Number of E cells
Ni = 1000 #Ni is the Number of I cells
T =  1000 #T is the simulation time (ms)
sqrtK = sqrt(K) 
taue = 10 #membrane time constant for exc. neurons (ms)
taui = 15 #membrane time constant for inh. neurons (ms)

jie = 4. / (taui*sqrtK)        #strength from E to I
jee = 10. / (taue*sqrtK)        #strength from E to E 

jei = -16. * 1.2 /(taue*sqrtK)      #strength from I to E
jii = -16. / (taui*sqrtK)      #strength from I to I 

#times, ns, Ne, Ncells, T = sim(K,Ne,Ni,T,jie,jei,jii,jee)

#println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
#println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne+1):Ncells] / T), " Hz")

#params = NetworkParameters(800, 4000, 1000, 1500, sqrt(800), 10, 15, 0.2 / (15 * sqrt(800)), 0.4 / (10 * sqrt(800)), -0.8 * 1.2 / (10 * sqrt(800)), -0.8 / (15 * sqrt(800)))
params = NetworkParameters(K, Ne, Ni, T, sqrtK, taue, taui, jie, jee, jei, jii)
#times, ns, Ne, Ncells, T = sim(params.K, params.Ne, params.Ni, params.T, params.jie, params.jei, params.jii, params.jee)
times, ns, Ne, Ncells, T = sim_new()
println("mean excitatory firing rate: ", mean(1000 * ns[1:params.Ne] / params.T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(params.Ne+1):Ncells] / params.T), " Hz")


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

        # Generate a timestamp for unique filenames
        timestamp = Dates.now()
        timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

    println("creating plot")

    # Define plot size: (width, height)
    plot_size = (600, 400) 

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
plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2)

fig_filename = "./figs/$timestamp_str.png"

# Save the figure with the timestamped filename
savefig(p2, fig_filename)

println("Figure saved as $fig_filename")
    


    # Define filenames for the figure and JSON file
    
    json_filename = "./figs/$timestamp_str.json"

    # Save the figure with the timestamped filename
   
    # Store the parameters in a JSON file with the timestamped filename
    param_dict = Dict(
        "K" => params.K,
        "Ne" => params.Ne,
        "Ni" => params.Ni,
        "T" => params.T,
        "sqrtK" => params.sqrtK,
        "taue" => params.taue,
        "taui" => params.taui,
        "jie" => params.jie,
        "jee" => params.jee,
        "jei" => params.jei,
        "jii" => params.jii
    )

    JSON3.open(json_filename, "w") do io
        JSON3.print(io, param_dict)
    end

    println("Parameters saved as $json_filename")
end




function run_simulation_and_tune_params(K, Ne, Ni, T, jie, jee, jei, jii)
    params = NetworkParameters(K, Ne, Ni, T, sqrt(K), 10, 15, jie, jee, jei, jii)
    
    times, ns, Ne, Ncells, T = sim(params.K, params.Ne, params.Ni, params.T, params.jie, params.jei, params.jii, params.jee)
    
    excitatory_rate = mean(1000 * ns[1:params.Ne] / params.T)
    inhibitory_rate = mean(1000 * ns[(params.Ne + 1):Ncells] / params.T)
    
    return times, excitatory_rate, inhibitory_rate
end


function tune_parameters()
    doplot = true
 
    successful_adjustments = 0
    required_successful_adjustments = 50
    # Initial parameter values
    K = 800
    Ne = 4000
    Ni = 1000
    T = 1500
    sqrtK = sqrt(K)
    taue = 10
    taui = 15

    jie = 0.2 / (taui * sqrtK)
    jee = 0.4 / (taue * sqrtK)
    jei = -0.8 * 1.2 / (taue * sqrtK)
    jii = -0.8 / (taui * sqrtK)
    
    adjustment_step = 0.001  # Adjust step size based on sensitivity
    random_factor = 0.001   # Introduce randomness
    while successful_adjustments < required_successful_adjustments
     while true
        println("Round+1")
        times, excitatory_rate, inhibitory_rate = run_simulation_and_tune_params(K, Ne, Ni, T, jie, jee, jei, jii)

        if excitatory_rate <= 5 && inhibitory_rate <= 5 && inhibitory_rate < excitatory_rate && inhibitory_rate >=0.5
            println("Parameters tuned successfully:")
            println("jie = $jie, jee = $jee, jei = $jei, jii = $jii")
            println("Mean excitatory firing rate: $excitatory_rate Hz")
            println("Mean inhibitory firing rate: $inhibitory_rate Hz")
            successful_adjustments=successful_adjustments+1
            if doplot
                timestamp = Dates.now()
                timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

                plot_size = (600, 400)

                window_size = 100
                step_size = 5

                e_rate = compute_sliding_rate(times[1:Ne, :], window_size, step_size, T)
                i_rate = compute_sliding_rate(times[(Ne + 1):Ncells, :], window_size, step_size, T)

                n_steps = length(e_rate)
                time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]

                p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)",
                          label="Excitatory", lw=2, linecolor=:red, size=plot_size)
                plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:skyblue)

                fig_filename = "./figs/ss_$timestamp_str.png"
                savefig(p2, fig_filename)
                println("Figure saved as $fig_filename")

                json_filename = "./figs/ss_$timestamp_str.json"
                param_dict = Dict(
                    "K" => K,
                    "Ne" => Ne,
                    "Ni" => Ni,
                    "T" => T,
                    "sqrtK" => sqrtK,
                    "taue" => taue,
                    "taui" => taui,
                    "jie" => jie,
                    "jee" => jee,
                    "jei" => jei,
                    "jii" => jii
                )
                JSON3.open(json_filename, "w") do io
                    JSON3.print(io, param_dict)
                end
                println("Parameters saved as $json_filename")
            end
        

        else
            # Randomly adjust parameters with noise
            jie += adjustment_step * (-1 + 2*random_factor * randn())
            jee += adjustment_step * (-1 + 2*random_factor * randn())
            jei -= adjustment_step * (-1 + 2*random_factor * randn())
            jii -= adjustment_step * (-1 + 2*random_factor * randn())

            jie = max(0.0001, jie)
            jee = max(0.0001, jee)
    
            # Ensure jei and jii are less than 0
            jei = min(-0.0001, jei)
            jii = min(-0.0001, jii)
        end
    end
 end
end


#tune_parameters()
