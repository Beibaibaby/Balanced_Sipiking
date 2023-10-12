#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information
function sim_old()
    println("Setting up parameters")

    # Network parameters
    Ncells = 5000
    Ne = 4000
    Ni = 1000
    T = 500  # Simulation time (ms)

    taue = 15  # Membrane time constant for excitatory neurons (ms)
    taui = 10

    # Connection probabilities
    pei = 0.5  # I -> E
    pie = 0.5  # E -> I
    pii = 0.5  # I -> I

    K = 800  # Average number of E->E connections per neuron
    sqrtK = sqrt(K)

    # Synaptic weights
    jie = 4.0 / (taui * sqrtK)
    jei = -18.0 * 1.2 / (taue * sqrtK)
    jii = -16 / (taui * sqrtK)

    jeeout = 10.0 / (taue * sqrtK)
    peeout = 0.2  # K / (Nepop * (ratiopee - 1) + Ne)

    # Stimulation
    Nstim = 400
    stimstr = 0 #0.7 / taue
    stimstart = T - 1500
    stimend = T

    # Constants and thresholds
    muemin = 1.1
    muemax = 1.2
    muimin = 1
    muimax = 1.05

    vre = 0.0
    threshe = 1
    threshi = 1

    dt = 0.1
    refrac = 5

    tauerise = 1
    taurise = 1
    tauedecay = 3
    tauidecay = 2

    maxrate = 100  # Maximum average firing rate (Hz)

    # Initialize parameters
    mu = zeros(Ncells)
    thresh = zeros(Ncells)
    tau = zeros(Ncells)

    mu[1:Ne] .= (muemax - muemin) .* rand(Ne) .+ muemin
    mu[(Ne + 1):Ncells] .= (muimax - muimin) .* rand(Ni) .+ muimin

    thresh[1:Ne] .= threshe
    thresh[(Ne + 1):Ncells] .= threshi

    tau[1:Ne] .= taue
    tau[(Ne + 1):Ncells] .= taui

    weights = zeros(Ncells, Ncells)

    # Random connections
    weights[1:Ne, 1:Ne] .= jeeout .* (rand(Ne, Ne) .< peeout)
    weights[1:Ne, (1 + Ne):Ncells] .= jei .* (rand(Ne, Ni) .< pei)
    weights[(1 + Ne):Ncells, 1:Ne] .= jie .* (rand(Ni, Ne) .< pie)
    weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii .* (rand(Ni, Ni) .< pii)

    for ci = 1:Ncells
        weights[ci, ci] = 0
    end

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)

    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    lastSpike = -100 * ones(Ncells)

    Nsteps = round(Int, T / dt)
    println("Starting simulation")

    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0

        for ci = 1:Ncells
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]

            synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end

            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)

                if v[ci] > thresh[ci]
                    v[ci] = vre
                    lastSpike[ci] = t
                    ns[ci] += 1
                    if ns[ci] <= maxTimes
                        times[ci, ns[ci]] = t
                    end

                    for j = 1:Ncells
                        if weights[j, ci] > 0  # E synapse
                            forwardInputsE[j] += weights[j, ci]
                        elseif weights[j, ci] < 0  # I synapse
                            forwardInputsI[j] += weights[j, ci]
                        end
                    end
                end
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
    end

    print("\r")
    times = times[:, 1:maximum(ns)]

    return times, ns, Ne, Ncells, T
end

function sim_new()
    println("Setting up parameters")

    # Network parameters
    Ncells = 5000
    Ne = 4000
    Ni = 1000
    T = 1000  # Simulation time (ms)

    taue = 15  # Membrane time constant for excitatory neurons (ms)
    taui = 10

    # Connection probabilities
    pei = 0.5  # I -> E
    pie = 0.5  # E -> I
    pii = 0.5  # I -> I

    K = 800  # Average number of E->E connections per neuron
    sqrtK = sqrt(K)

    # Synaptic weights
    jie_A = 4.0 / (taui * sqrtK)
    jei_A = -18.0 * 1.2 / (taue * sqrtK)
    jii_A = -16 / (taui * sqrtK)

    jeeout_A = 10.0 / (taue * sqrtK)
    peeout = 0.2  # K / (Nepop * (ratiopee - 1) + Ne)

    # Stimulation
    Nstim = 400
    stimstr = 0 #0.7 / taue
    stimstart = T - 1500
    stimend = T
    d = 0.15       # depression fraction upon a spike
    f = 0.92         # facilitation increment upon a spike
    tau_d = 103.0     # time constant for D to recover to 1 (ms)
    tau_f = 96.0     # time constant for F to recover to 1 (ms)
    # Constants and thresholds
    muemin = 1.1
    muemax = 1.2
    muimin = 1
    muimax = 1.05

    vre = 0.0
    threshe = 1
    threshi = 1

    dt = 0.1
    refrac = 5

    tauerise = 1
    taurise = 1
    tauedecay = 3
    tauidecay = 2

    maxrate = 100  # Maximum average firing rate (Hz)

    # Initialize parameters
    mu = zeros(Ncells)
    thresh = zeros(Ncells)
    tau = zeros(Ncells)

    mu[1:Ne] .= (muemax - muemin) .* rand(Ne) .+ muemin
    mu[(Ne + 1):Ncells] .= (muimax - muimin) .* rand(Ni) .+ muimin

    thresh[1:Ne] .= threshe
    thresh[(Ne + 1):Ncells] .= threshi

    tau[1:Ne] .= taue
    tau[(Ne + 1):Ncells] .= taui

    weights = zeros(Ncells, Ncells)
    weights_D=ones(Ncells , Ncells)
    weights_F=ones(Ncells , Ncells)
    # Random connections
    weights[1:Ne, 1:Ne] .= jeeout_A .* (rand(Ne, Ne) .< peeout)
    weights[1:Ne, (1 + Ne):Ncells] .= jei_A .* (rand(Ne, Ni) .< pei)
    weights[(1 + Ne):Ncells, 1:Ne] .= jie_A .* (rand(Ni, Ne) .< pie)
    weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii_A .* (rand(Ni, Ni) .< pii)

    for ci = 1:Ncells
        weights[ci, ci] = 0
    end

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)

    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)
    
    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    lastSpike = -100 * ones(Ncells)

    Nsteps = round(Int, T / dt)
    println("Starting simulation")
    
    for ti = 1:Nsteps

        if mod(ti,Nsteps/100) == 1  #print percent complete
			print("\r",round(Int,100*ti/Nsteps))
		end

        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0


        #weights_D[(1 + Ne):Ncells, 1:Ne] .+= (1 .- weights_D[(1 + Ne):Ncells, 1:Ne]) ./ tau_d
        #weights_F[(1 + Ne):Ncells, 1:Ne] .+= (1 .- weights_F[(1 + Ne):Ncells, 1:Ne]) ./ tau_d
        #weights[(1 + Ne):Ncells, 1:Ne] .*= weights_D[(1 + Ne):Ncells, 1:Ne] .* weights_F[(1 + Ne):Ncells, 1:Ne] 
        weights_D .+= (1 .- weights_D) ./ tau_d
        weights_F .+= (1 .- weights_F) ./ tau_f
        weights .*= weights_D .* weights_F
        for ci = 1:Ncells
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]

            synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end

            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)

                if v[ci] > thresh[ci]

                    v[ci] = vre
                    lastSpike[ci] = t
                    ns[ci] += 1
                    if ns[ci] <= maxTimes
                        times[ci, ns[ci]] = t
                    end

                    for j = 1:Ncells

                        if weights[j, ci] > 0  # E synapse
                            weights_D[j, ci] = d * weights_D[j, ci]
                            if j > Ne 
                             weights_F[j, ci] = f + weights_F[j, ci]
                            end 
                            forwardInputsE[j] += weights[j, ci]

                        elseif weights[j, ci] < 0  # I synapse
                            forwardInputsI[j] += weights[j, ci]
                        end
                    end
                 
                

                end

            end

        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
    end

    print("\r")
    times = times[:, 1:maximum(ns)]

    return times, ns, Ne, Ncells, T
end

              