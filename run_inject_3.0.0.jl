using Statistics
using Plots  
using Dates  # for generating timestamps
using JSON3  # Use JSON3 package for JSON handling
using Random
using JLD2
using Distributed
using SharedArrays
using Measures
using Distributions

#print("Number of threads: ")
#println(Sys.CPU_THREADS - 1)
#addprocs(Sys.CPU_THREADS - 1)

#using Profile
#@everywhere using SharedArrays
#@everywhere include("sim_inject_3.0.0.jl") 
using SharedArrays
include("sim_inject_3.0.0.jl") 

struct NetworkParameters
    Ncells::Int
    Ne::Int
    Ni::Int
    T::Int
    taue::Int
    taui::Int
    pei::Float64
    pie::Float64
    pii::Float64
    pee::Float64
    K::Int
    jie_para::Float64
    jei_para::Float64
    jii_para::Float64
    jee_para::Float64
    Nstim::Int
    stimstr_para::Float64
    d::Float64
    f::Float64
    stim_duration::Int
    stim_start_time::Int
    ie_sign::Bool
    ee_sign::Bool
    corr_flag::Bool
    add_noise::Bool
    lambda_noise::Float64
    scale_noise::Float64
end

# Define a function to retrieve a value from ARGS or return a default value if not present.
function get_arg(key, default)
    index = findfirst(==(key), ARGS)
    if index !== nothing && index < length(ARGS)
        return ARGS[index + 1]
    end
    return default
end


function run_experiment(;
    Ncells,
    Ne,
    Ni,
    T,
    taue,
    taui,
    pei,
    pie,
    pii,
    pee,
    K,
    jie,
    jei,
    jii,
    jee,
    Nstim,
    stimstr,
    d,
    f,
    stim_duration,
    stim_start_time,
    ie_sign,
    ee_sign,
    corr_flag,
    low_plot,
    add_noise,
    lambda_noise,
    scale_noise,
    env
)
        
        doplot = true
        do_v_corr = false
        do_repeat_v_corr=false

        #Setting the parameters
        ############################
        # Now, use the provided values to create an instance of the struct:
        params = NetworkParameters(Ncells, Ne, Ni, T, taue, taui, pei, pie, pii, pee, K, jie, jei, jii, jee, Nstim, stimstr,d,f,stim_duration, stim_start_time,ie_sign,ee_sign,corr_flag, add_noise,
        lambda_noise,scale_noise)

        #store it
        #run the stimulus
        #times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights = sim_old()
        global times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_mean, weights_F_mean,weights_IE_mean_history,weights_EE_mean_history
        times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_mean, weights_F_mean,weights_IE_mean_history,weights_EE_mean_history=sim_dynamic(
            params.Ne,params.Ni,params.T,params.taue,params.taui,params.pei,params.pie,params.pii,params.pee,params.K,
            params.stimstr_para,params.Nstim,params.jie_para,params.jei_para,params.jii_para,params.jee_para,params.d,
            params.f,params.stim_duration,params.stim_start_time,params.ie_sign,params.ee_sign,params.corr_flag,
            params.add_noise,
            params.lambda_noise,
            params.scale_noise)
        println("mean excitatory firing rate: ", mean(1000 * ns[1:params.Ne] / params.T), " Hz")
        println("mean inhibitory firing rate: ", mean(1000 * ns[(params.Ne+1):Ncells] / params.T), " Hz")
        product_weights = weights_D_mean .* weights_F_mean




        if doplot 

            # Generate a timestamp for unique filenames
                timestamp = Dates.now()
                timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

                println("creating plot")
           
                # Define plot size: (width, height)
                plot_size = (1000, 600) 
                plot_margin = 50mm

                # Parameters for sliding window
                window_size = 100  # in ms
                step_size = 5     # in ms

                #print(Ne)
                e_rate = compute_sliding_rate(times[1:params.Ne, :], window_size, step_size, params.T)
                i_rate = compute_sliding_rate(times[(params.Ne+1):Ncells, :], window_size, step_size, params.T)


                # Compute the time values based on window_size and step_size
                n_steps = length(e_rate)  # or length(i_rate), assuming they have the same length
                #time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]
                time_values = [i * step_size + window_size  for i in 1:n_steps]

                # Add a code to detect low rate or not 
                
                if low_plot ##this para control whether focus on the zoom in low activity
                    p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="Firing rate (d=$d), (f=$f)", ylim=(0,5))
                    plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,left_margin=plot_margin)
                else
                    p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firinçg rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="Firing rate (d=$d), (f=$f)")
                    plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,left_margin=plot_margin)
                end 
                
                if env == 1
                    dir_name ="/root/autodl-tmp/$timestamp_str"
                elseif env == 2
                   dir_name = "../figs_paras/$timestamp_str"
                else
                    dir_name = "../figs_paras/$timestamp_str"
                end

                
                #Create the directory for results
                if !isdir(dir_name)
                    mkdir(dir_name)
                end
                
                @save "$dir_name/times.jld2" times
                @save "$dir_name/ns.jld2" ns
                @save "$dir_name/Ne.jld2" Ne
                @save "$dir_name/Ncells.jld2" Ncells
                @save "$dir_name/T.jld2" T
                @save "$dir_name/v_history.jld2" v_history
                @save "$dir_name/E_input.jld2" E_input
                @save "$dir_name/I_input.jld2" I_input
                @save "$dir_name/weights.jld2" weights
                @save "$dir_name/weights_D_mean.jld2" weights_D_mean
                @save "$dir_name/weights_F_mean.jld2" weights_F_mean
                @save "$dir_name/weights_IE_mean_history.jld2" weights_IE_mean_history
                @save "$dir_name/weights_EE_mean_history.jld2" weights_EE_mean_history
                


                fig_filename = "$dir_name/plot_FR_$timestamp_str.png"
                savefig(p2, fig_filename)

                # Generate scaled x-values
                x_values = 0:0.1:(length(product_weights) - 1) * 0.1

                # Create individual subplots
                # Adjust title font size and margins for the plots
                title_fontsize = 12  # You can adjust this value as needed
                left_margin = 10   # Adjust this value as needed

                # Create individual subplots with the modifications
                fig_width = 1800  # Width in pixels
                fig_height = 600  # Height in pixels
                title_fontsize = 14  # Title font size
                

                p3 = plot(
                    x_values, weights_EE_mean_history, 
                    label="EE Mean History", 
                    xlabel="Time (ms)", 
                    ylabel="Value", 
                    title="E->E strength (d=$d), (f=$f)", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # And also here in pixels
                )
                
                p4 = plot(
                    x_values, weights_IE_mean_history, 
                    label="IE Mean History", 
                    xlabel="Time (ms)", 
                    ylabel="Value", 
                    title="E->I strength (d=$d), (f=$f)", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # And also here in pixels
                )
                
                # Now combine the subplots
                p5 = plot(p3, p4, layout=(2,1), size=(fig_width, fig_height*2))

                # Save the plot as an image
                savefig(p5,"$dir_name/plot_syn_$timestamp_str.png") 

                println("Combined figure saved as $dir_name/plot_SYN_$timestamp_str.png") 

                p6 = plot(
                    x_values, product_weights, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="Product Value", 
                    title="Product of D and F (d=$d), (f=$f)", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p7 = plot(
                    x_values, weights_D_mean, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="D Value", 
                    title="D factor (d=$d), (f=$f)", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p8 = plot(
                    x_values, weights_F_mean, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="F Value", 
                    title="F factor (d=$d), (f=$f)", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                combined_plot = plot(p6, p7, p8, layout = (3, 1), size = (fig_width*1.2, fig_height*3.2))

                savefig(combined_plot,"$dir_name/plot_DF_$timestamp_str.png") 


                println("Figure saved as $fig_filename")  
                json_filename = "$dir_name/$timestamp_str.json"




                param_dict = Dict(
                    "Ncells" => params.Ncells,
                    "Ne" => params.Ne,
                    "Ni" => params.Ni,
                    "T" => params.T,
                    "taue" => params.taue,
                    "taui" => params.taui,
                    "pei" => params.pei,
                    "pie" => params.pie,
                    "pii" => params.pii,
                    "pee" => params.pee,
                    "K" => params.K,
                    "jie_para" => params.jie_para,
                    "jee_para" => params.jee_para,
                    "jei_para" => params.jei_para,
                    "jii_para" => params.jii_para,
                    "Nstim" => params.Nstim,
                    "stimstr_para" => params.stimstr_para,
                    "d" => params.d,
                    "f" => params.f,
                    "stim_duration" => params.stim_duration,
                    "stim_start_time" => params.stim_start_time,
                    "ie_sign"=> params.ie_sign,
                    "ee_sign"=> params.ee_sign,
                    "corr_flag" =>  params.corr_flag,
                    "low_plot" => params.low_plot,
                    "lambda_noise" => params.lambda_noise,
                    "add_noise" => params.add_noise,
                    "scale_noise" => params.scale_noise
                )

                # Now, you can access any of these values using the dictionary's keys, e.g., param_dict["Ne"] or param_dict["jie"].


                JSON3.open(json_filename, "w") do io
                JSON3.print(io, param_dict)
                end

                println("Parameters saved as $json_filename")
                println("done")
                #println(cross_EPSP_EPSP)
                #println("cross E-E: ", cross_corr_E_E)

        end



        if do_v_corr #if compute the correlation

            EPSP_EPSP_pool=v_history[1:100,end-9999:end]
            TT_pool=v_history[501:600,end-9999:end]
            IPSP_IPSP_pool=v_history[1001:1100,end-9999:end]

            IPSP_IPSP_pool= IPSP_IPSP_pool .- 2
            TT_pool=TT_pool .- 1
            EPSP_EPSP_pool=EPSP_EPSP_pool .+ 0.2

            # Initialize accumulators
            EPSP_EPSP_accumulator = zeros(100, 10000)
            TT_accumulator = zeros(100, 10000)
            IPSP_IPSP_accumulator = zeros(100, 10000)

            if do_repeat_v_corr

            n_run=20
            for iter in 1:n_run
                println("Current iteration: $iter")  # This line prints the current iteration number
                local times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, EPSP_EPSP_pool, TT_pool, IPSP_IPSP_pool
                times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights=sim_working(params.Ne,params.Ni,params.T,params.taue,params.taui,params.pei,params.pie,params.pii,params.pee,params.K,params.stimstr_para,params.Nstim,params.jie_para,params.jei_para,params.jii_para,params.jee_para)

                EPSP_EPSP_pool = v_history[1:100, end-9999:end]
                TT_pool = v_history[501:600, end-9999:end]
                IPSP_IPSP_pool = v_history[1001:1100, end-9999:end]
                
                IPSP_IPSP_pool = IPSP_IPSP_pool .- 1 #hard?
                TT_pool = TT_pool .- 1
                EPSP_EPSP_pool = EPSP_EPSP_pool .+ 0.2
                
                # Add to accumulator
                EPSP_EPSP_accumulator .+= EPSP_EPSP_pool
                TT_accumulator .+= TT_pool
                IPSP_IPSP_accumulator .+= IPSP_IPSP_pool
            end

            # Compute the average
            EPSP_EPSP_avg = EPSP_EPSP_accumulator ./ n_run
            TT_avg = TT_accumulator ./ n_run
            IPSP_IPSP_avg = IPSP_IPSP_accumulator ./ n_run

            IPSP_IPSP_pool= EPSP_EPSP_avg
            TT_pool=TT_avg
            EPSP_EPSP_pool=IPSP_IPSP_avg

            end


            avg_correlation_E_I=compute_correlation(E_input, I_input)
            avg_correlation_E_E=compute_correlation(E_input, E_input)
            avg_correlation_I_I=compute_correlation(I_input, I_input)
            cross_corr_E_E=compute_cross_correlation(E_input[:, end-999:end], E_input[:, end-999:end])
            cross_corr_I_I=compute_cross_correlation(I_input[:, end-999:end], I_input[:, end-999:end])
            cross_corr_E_I=compute_cross_correlation(E_input[:, end-999:end], I_input[:, end-999:end])
            cross_corr_I_E=compute_cross_correlation(I_input[:, end-999:end], E_input[:, end-999:end])


            println("avg correlation(E-I): ", avg_correlation_E_I)
            println("avg correlation(E-E): ", avg_correlation_E_E)
            println("avg correlation(I-I): ", avg_correlation_I_I)

            cross_EPSP_EPSP=compute_cross_correlation(EPSP_EPSP_pool,EPSP_EPSP_pool)
            cross_IPSP_IPSP=compute_cross_correlation(IPSP_IPSP_pool,IPSP_IPSP_pool)
            cross_EPSP_IPSP=compute_cross_correlation(EPSP_EPSP_pool,IPSP_IPSP_pool)
            cross_IPSP_EPSP=compute_cross_correlation(IPSP_IPSP_pool,EPSP_EPSP_pool)
            cross_T_T=compute_cross_correlation(TT_pool,TT_pool)

            #plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I,cross_corr_I_E)
            plot_correlations_mem(cross_EPSP_EPSP, cross_IPSP_IPSP, cross_T_T, cross_EPSP_IPSP, cross_IPSP_EPSP)

            plot_cells(v_history, [1, 505, 1005])

            println("finish")

        end


end



# Extract parameters from ARGS or use default values
Ncells = parse(Int, get_arg("--Ncells", "5000"))
Ne = parse(Int, get_arg("--Ne", "4000"))
Ni = parse(Int, get_arg("--Ni", "1000"))
T = parse(Int, get_arg("--T", "3000"))
taue = parse(Int, get_arg("--taue", "15"))
taui = parse(Int, get_arg("--taui", "10"))
pei = parse(Float64, get_arg("--pei", "0.5"))
pie = parse(Float64, get_arg("--pie", "0.5"))
pii = parse(Float64, get_arg("--pii", "0.5"))
pee = parse(Float64, get_arg("--pee", "0.2"))
K = parse(Int, get_arg("--K", "800"))
jie = parse(Float64, get_arg("--jie", "4.0"))
jei = parse(Float64, get_arg("--jei", string(-18.0 * 1.2)))
jii = parse(Float64, get_arg("--jii", "-16.0"))
jee = parse(Float64, get_arg("--jee", "10.0"))
Nstim = parse(Int, get_arg("--Nstim", "4000"))
stimstr = parse(Float64, get_arg("--stimstr", "0.0"))
d = parse(Float64, get_arg("--d", "0.15"))
f = parse(Float64, get_arg("--f", "0.92"))
stim_duration= parse(Int, get_arg("--stim_duration", "2"))
stim_start_time= parse(Int, get_arg("--stim_start_time", "1000"))

ie_sign = parse(Bool, get_arg("--ie_sign", "true")) #controal E->I is dynamic or not 
ee_sign = parse(Bool, get_arg("--ee_sign", "true")) #controal E->E is dynamic or not 
corr_flag = parse(Bool, get_arg("--corr_flag", "false")) ##wether compute and plot EPSP and IPSP
low_plot = parse(Bool, get_arg("--low_plot", "false")) #contronl whether manully plot a low ativity regime
lambda_noise = parse(Float64, get_arg("--lambda_noise", "2.0"))
add_noise = parse(Bool, get_arg("--add_noise", "true"))
scale_noise = parse(Float64, get_arg("--scale_noise", "0.7"))
env = scale_noise = parse(Int, get_arg("--env", "1"))

run_experiment(;Ncells,
    Ne,
    Ni,
    T,
    taue,
    taui,
    pei,
    pie,
    pii,
    pee,
    K,
    jie,
    jei,
    jii,
    jee,
    Nstim,
    stimstr,
    d,
    f,
    stim_duration,
    stim_start_time,
    ie_sign,
    ee_sign,
    corr_flag,
    low_plot,
    add_noise,
    lambda_noise,
    scale_noise,
    env
)