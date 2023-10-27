#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information
using Statistics
using ProgressMeter

##oldest version
function sim_old()
    println("Setting up parameters")

    # Network parameters
    Ncells = 5000
    Ne = 4000
    Ni = 1000
    T = 2000  # Simulation time (ms)

    taue = 15  # Membrane time constant for excitatory neurons (ms)
    taui = 10 # Membrane time constant for I neurons (ms)

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
    stimstr = 0 #0.1/taue #0.7 / taue
    stimstart = T - 1100
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
    Nsteps = round(Int, T / dt)
    v_history = zeros(Ncells, Nsteps)  # Nsteps because we're saving at each time step, not just spikes
    E_input=zeros(Ncells, Nsteps) 
    I_input=zeros(Ncells, Nsteps)
    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    v[1:3] .=0
    #v = vre*ones(Ncells)
    lastSpike = -100 * ones(Ncells)

    
    println("Starting simulation")
    

    corr_pairs=100
    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0

        for ci = 1:Ncells
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]

            #synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synInput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            E_input[ci, ti] = synInput_E
            I_input[ci, ti] = synInput_I 
            synInput = synInput_E+synInput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            
            if ci<corr_pairs
                synInput=synInput_E-(mu[ci]+0.2)/tau[ci]
                #synInput = synInput_E
            end

            if ci >= 500 && ci <600
                synInput=synInput+(1-mu[ci])/tau[ci]
            end

            if ci >= 1001 && ci <1100
                synInput=synInput_I+(2-mu[ci])/tau[ci]
            end

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end
           
            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)     
                v_history[ci, ti] = v[ci]
                if v[ci] > thresh[ci]
                    if (ci>=corr_pairs && ci<500 )||(ci>=600 && ci<=1000)  || ci>=1100
                     v[ci] = vre
                    end
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
            else
                v_history[ci, ti] = v[ci]#
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
    end

    print("\r")
    println(maximum(ns))
    println(maximum(ns))
    if maximum(ns)<=maxTimes
         times = times[:, 1:maximum(ns)]
    else
       println("triger")
    end

    return times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights
end

##oldest version + para setting
function sim_working(Ne,Ni,T,taue,taui,pei,pie,pii,pee,K,stimstr_para,Nstim,jie_para,jei_para,jii_para,jee_para)
    println("Setting up parameters")


    # Network parameters

    Ncells = Ne+Ni
    sqrtK = sqrt(K)

    # Synaptic weights
    jie = jie_para / (taui * sqrtK)
    jei = jei_para / (taue * sqrtK)
    jii = jii_para / (taui * sqrtK)
    jee = jee_para / (taue * sqrtK)
   

    # Stimulation
    stimstr = stimstr_para/taue 
    stimstart = T - 1100
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
    weights[1:Ne, 1:Ne] .= jee .* (rand(Ne, Ne) .< pee)
    weights[1:Ne, (1 + Ne):Ncells] .= jei .* (rand(Ne, Ni) .< pei)
    weights[(1 + Ne):Ncells, 1:Ne] .= jie .* (rand(Ni, Ne) .< pie)
    weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii .* (rand(Ni, Ni) .< pii)

    for ci = 1:Ncells
        weights[ci, ci] = 0
    end

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)
    Nsteps = round(Int, T / dt)
    v_history = zeros(Ncells, Nsteps)  # Nsteps because we're saving at each time step, not just spikes
    E_input=zeros(Ncells, Nsteps) 
    I_input=zeros(Ncells, Nsteps)
    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    v[1:3] .=0
    #v = vre*ones(Ncells)
    lastSpike = -100 * ones(Ncells)

    
    println("Starting simulation")
    

    corr_pairs=100
    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0

        for ci = 1:Ncells
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]

            #synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synInput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            E_input[ci, ti] = synInput_E
            I_input[ci, ti] = synInput_I 
            synInput = synInput_E+synInput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            
            if ci<corr_pairs
                synInput=synInput_E-(mu[ci]+0.2)/tau[ci]
                #synInput = synInput_E
            end

            if ci >= 500 && ci <600
                synInput=synInput+(1-mu[ci])/tau[ci]
            end

            if ci >= 1001 && ci <1100
                synInput=synInput_I+(2-mu[ci])/tau[ci]
            end

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end
           
            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)     
                v_history[ci, ti] = v[ci]
                if v[ci] > thresh[ci]
                    if (ci>=corr_pairs && ci<500 )||(ci>=600 && ci<=1000)  || ci>=1100
                     v[ci] = vre
                    end
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
            else
                v_history[ci, ti] = v[ci]#
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
    end

    print("\r")
    println(maximum(ns))
    println(maximum(ns))
    if maximum(ns)<=maxTimes
         times = times[:, 1:maximum(ns)]
    else
       println("triger")
    end

    return times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights
end

##oldest version + para setting+ E->I dyanmic 
function sim_dynamic_EI(Ne,Ni,T,taue,taui,pei,pie,pii,pee,K,stimstr_para,Nstim,jie_para,jei_para,jii_para,jee_para,d,f,stim_duration,stim_start_time)
    println("Setting up parameters")

    # Network parameters

    Ncells = Ne+Ni
    sqrtK = sqrt(K)

    # Synaptic weights
    jie = jie_para / (taui * sqrtK)
    jei = jei_para / (taue * sqrtK)
    jii = jii_para / (taui * sqrtK)
    jee = jee_para / (taue * sqrtK)
   

    # Stimulation
    stimstr = stimstr_para/taue 
    stimstart = stim_start_time
    stimend = stim_start_time+stim_duration

    #d = 0.15       # depression fraction upon a spike
    #f = 0.92  (orignal)      # facilitation increment upon a spike
    #f = 0.5
    #tau_d = 103.0     # time constant for D to recover to 1 (ms)
    #tau_f = 96.0     # time constant for F to recover to 1 (ms)
    
    tau_d = 1030 # time constant for D to recover to 1 (step;not ms)
    tau_f = 960 # time constant for D to recover to 1 (step;not ms)

    #todo: might need adjust the time scale

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
    weights_D = ones(Ne) #the sending D
    weights_F = ones(Ne) #the sending F

    # Here we only need one decay/facilitation factor for one given neuron i, the factors from i to j are all the same

    # Random connections
    weights[1:Ne, 1:Ne] .= jee .* (rand(Ne, Ne) .< pee)
    weights[1:Ne, (1 + Ne):Ncells] .= jei .* (rand(Ne, Ni) .< pei)
    weights[(1 + Ne):Ncells, 1:Ne] .= jie .* (rand(Ni, Ne) .< pie)
    weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii .* (rand(Ni, Ni) .< pii)

    for ci = 1:Ncells
        weights[ci, ci] = 0
    end

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)
    Nsteps = round(Int, T / dt)


    weights_D_mean = zeros(Nsteps)  
    weights_F_mean = zeros(Nsteps)

    v_history = zeros(Ncells, Nsteps)  # Nsteps because we're saving at each time step, not just spikes
    E_input=zeros(Ncells, Nsteps) 
    I_input=zeros(Ncells, Nsteps)
    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)
    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    v[1:3] .=0
    #v = vre*ones(Ncells)
    lastSpike = -100 * ones(Ncells)

    
    println("Starting simulation")
    pl = Progress(Nsteps, 5)
    
    corr_pairs=100

    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0
        
        weights_D .+= (1 .- weights_D) ./ tau_d
        weights_F .+= (1 .- weights_F) ./ tau_f

        weights_D_mean[ti] = mean(weights_D)
        weights_F_mean[ti] = mean(weights_F)

        W_sub_view = @view weights[end-Ne+1:end, 1:Ne]
        W_sub_view .= W_sub_view .* (weights_D .* weights_F)


        for ci = 1:Ncells
            

            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]
   
            #synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synInput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            E_input[ci, ti] = synInput_E
            I_input[ci, ti] = synInput_I 
            synInput = synInput_E+synInput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            
            if ci<corr_pairs
                synInput=synInput_E-(mu[ci]+0.2)/tau[ci]
                #synInput = synInput_E
            end

            if ci >= 500 && ci <600
                synInput=synInput+(1-mu[ci])/tau[ci]
            end

            if ci >= 1001 && ci <1100
                synInput=synInput_I+(2-mu[ci])/tau[ci]
            end

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
            end
           
            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)     
                v_history[ci, ti] = v[ci]
                if v[ci] > thresh[ci]
                    if (ci>=corr_pairs && ci<500 )||(ci>=600 && ci<=1000)  || ci>=1100
                     v[ci] = vre
                    end
                    lastSpike[ci] = t
                    ns[ci] += 1
                    if ns[ci] <= maxTimes
                        times[ci, ns[ci]] = t
                    end
                    
                    #test if remove
                   
                    weights_D[ci] = d * weights_D[ci]
                    weights_F[ci] = f + weights_F[ci]
                                                

                    for j = 1:Ncells
                        if weights[j, ci] > 0  # E synapse
                            forwardInputsE[j] += weights[j, ci]
                        elseif weights[j, ci] < 0  # I synapse
                            forwardInputsI[j] += weights[j, ci]
                        end
                    end
                end
            else
                v_history[ci, ti] = v[ci]#
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
        next!(pl)
    end

    print("\r")
    println(maximum(ns))
    println(maximum(ns))
    if maximum(ns)<=maxTimes
         times = times[:, 1:maximum(ns)]
    else
       println("no over max")
    end
    
    return times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_mean, weights_F_mean
end



function compute_correlation(E_input::Matrix, I_input::Matrix, num_pairs::Int=100)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)

    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells
        error("The number of pairs cannot exceed the number of cells.")
    end

    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    total_correlation = 0.0

    for i = 1:num_pairs
        e_data = E_input[random_e_indices[i], :]
        i_data = I_input[random_i_indices[i], :]

        # Compute Pearson correlation and add to the total
        total_correlation += cor(e_data, i_data)
    end

    # Compute average correlation
    avg_correlation = total_correlation / num_pairs

    return avg_correlation
end

function compute_cross_correlation(E_input::Matrix, I_input::Matrix, tau_range::UnitRange{Int}=-5:5, num_pairs::Int=50)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)

    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells * Ncells
        error("The number of pairs cannot exceed the square of the number of cells.")
    end


    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    avg_correlations = Dict{Int, Float64}()

    for tau in tau_range
        total_correlation = 0.0
        for i = 1:num_pairs
            e_data = E_input[random_e_indices[i], :]
            i_data = I_input[random_i_indices[i], :]
            
            # Shift data according to tau
            if tau > 0
                e_data = e_data[1:end-tau*10]
                i_data = i_data[tau*10+1:end]
            elseif tau < 0
                e_data = e_data[-tau*10+1:end]
                i_data = i_data[1:end+tau*10]
            end
            
            # Compute Pearson correlation for this tau and add to the total
            total_correlation += cor(e_data, i_data)
        end

        # Compute average correlation for this tau
        avg_correlations[tau] = total_correlation / num_pairs
    end

    return avg_correlations
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
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synIput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput = synInput_E+synIput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)

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

function sim_new_new()
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
    
    @Threads.threads for ti = 1:Nsteps
        if mod(ti,Nsteps/100) == 1  #print percent complete
            print("\r", round(Int,100*ti/Nsteps))
        end

        t = dt * ti
        forwardInputsE .= 0.0
        forwardInputsI .= 0.0
        
        weights_D .+= (1 .- weights_D) ./ tau_d
        weights_F .+= (1 .- weights_F) ./ tau_f
        weights .*= weights_D .* weights_F
        
        # Note: This loop could potentially be parallelized further 
        # but needs to ensure no race conditions or false sharing.
        for ci = 1:Ncells
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]
            
            # [Truncated: Remaining ci-loop code...]
            # Ensure in-place operations, pre-compute constants, and
            # be cautious with conditional checks in the loop.
            
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

function compute_sliding_rate_old(spiketimes, window_size, step_size, T)
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

function compute_sliding_rate(spiketimes, window_size, step_size, T)
    println("starting computing firing rate")
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    p = Progress(n_steps, 2) # Initialize the progress bar

    for i = 1:n_steps
        next!(p) # Update the progress bar

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




function plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I,cross_corr_I_E)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    
    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")
    savefig("cross_correlation_plot_PSP_new.png")
 
end

function plot_correlations_mem(cross_corr_E_E, cross_corr_I_I, cross_corr_T_T, cross_corr_E_I, cross_corr_I_E)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_T_T = [cross_corr_T_T[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_T_T, label="T-T", linewidth=2, marker=:circle, color=:green)

    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")
    savefig("cross_correlation_plot_PSP_mem.png")

end

function plot_cells(v_history, cells)
    p = plot(layout=(3,1), size=(600,800))

    for (index, cell) in enumerate(cells)
        plot!(p[index], v_history[cell, :], label="Cell $cell", xlabel="Time", ylabel="Voltage", legend=:topright)
    end
    savefig("v_history.png")
end
