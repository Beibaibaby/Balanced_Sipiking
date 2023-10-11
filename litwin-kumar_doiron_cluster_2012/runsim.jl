#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

using Statistics
using Plots  # use Plots instead of PyPlot

#uncomment the line below and set doplot=true to plot a raster
#using Plots  # this line is not necessary if the line above is active
doplot = true

include("sim.jl")

times, ns, Ne, Ncells, T = sim()

println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne+1):Ncells] / T), " Hz")

if doplot 
    println("creating plot")
    # The plotting might differ slightly due to syntactical differences between Plots and PyPlot.
    p = scatter(legend=false, xlabel="Time", ylabel="Neuron", ylims=(0, Ne), xlims=(0, T))
    for ci = 1:Ne
        vals = times[ci, 1:ns[ci]]
        y = ci .* ones(length(vals))
        scatter!(vals, y, markersize=1, markercolor=:black, markerstrokecolor=:black, legend=false, show=false)
    end
    # Save plot as a PNG file
    savefig(p, "output.png")
end

