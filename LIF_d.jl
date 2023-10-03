using Plots

# LIF neuron parameters
const τ_m = 10.0        # membrane time constant (ms)
const V_thresh = -50.0   # spike threshold (mV)
const V_reset = -65.0    # reset potential (mV)
const V_rest = -65.0     # resting potential (mV)
const R_m = 10.0         # membrane resistance (MΩ)

# Synaptic weight parameters
const A = 10.0          # fixed parameter A
const F = 1.0           # fixed parameter F
const d = 0.5           # depression fraction upon a spike
const tau_d = 5.0      # time constant for D to recover to 1 (ms)

# Simulation parameters
const dt = 1.0          # time step (ms)
const T = 100.0          # total time to simulate (ms)
const time = 0:dt:T      # time vector

# Initialize variables
V = V_rest               # membrane potential
Vs = Float64[]           # membrane potential at each time step
spike_times = []         # vector to record spike times
S = zeros(length(time))  # spike train
D = 1.0                  # dynamic variable D

# Generate an artificial spike train: Spikes at 35 ms, 45 ms, and 55 ms
S[convert(Int, 35/dt)] = 1
S[convert(Int, 45/dt)] = 1
S[convert(Int, 55/dt)] = 1

# Main simulation loop
for (idx, t) in enumerate(time)
    global V, D
    
    # Update synaptic weight W
    W = A * D * F
    
    # LIF equation with synaptic weight
    dV = (-(V - V_rest) + R_m * W * S[idx]) / τ_m
    V += dV * dt
    
    # Check for spike
    if V >= V_thresh
        push!(spike_times, t)
        V = V_reset
    end
    
    # Update D
    if S[idx] == 1
        D *= d  # Depression upon spike
    else
        # Recovery towards 1 when no spike
        dD = (1 - D) / tau_d
        D += dD * dt
    end
    
    push!(Vs, V)
end

# Plotting
p1 = plot(time, Vs, label="Membrane Potential", xlabel="Time (ms)", ylabel="Membrane Potential (mV)", ylims=(-70, -40), legend=:topright, legendfontsize=8)
scatter!(spike_times, V_thresh*ones(length(spike_times)), markershape=:circle, color="red", label="Spikes")
title!("LIF Neuron Spiking Activity")

p2 = plot(time, S, label="Spike Train", xlabel="Time (ms)", ylabel="Amplitude", color="purple", legend=:topright, ylims=(-0.1, 1.1), legendfontsize=8, linewidth=2)
title!("Input Spike Train")

# Combine the plots vertically
plot(p1, p2, layout=(2,1), link=:x)

# Save the plot
savefig("./figs/only_d.png")
