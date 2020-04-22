# Introduce packages
using SpecialFunctions # For generating gamma distribution
using LinearAlgebra
using Distributions
using StatsBase # To use sample
using Statistics
using Optim # To use optimize
using Roots
using CSV
using DataFrames
using DelimitedFiles
using Plots

include("sbm_fn.jl")
include("trial_fn.jl")

function fn_update_R0(Gc,N,T_arr)
    k = sum(sum(Gc))/N # mean degree of the network
    T = mean(T_arr[T_arr .> 0]) # Average probability that an infectious individual will transmit the disease to a susceptible individual with whom they haVE_true contact
    R0 = T * (k^2/k - 1)
    return k, T, R0
end

function fn_update_vars(samplesize, samplesize_mat, n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_infectious_people_mat, TTE_mat, VE_true_mat, R0_mat, treatment_gp, timestep_fn, alpha, power, nstatus, tstatus, trial_begintime, trial_endtime, gamma_infectperiod_maxduration, Gc, N, T_arr)
    #println("within function (samplesize): ",samplesize)
    #println("within function (samplesize_mat): ",samplesize_mat)
    #println("n_infectious_people_mat 1: ", n_infectious_people_mat)
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    #println("n_infectious_control, n_infectious_treatment: ", n_infectious_control, n_infectious_treatment)
    n_infectious_people = n_infectious_control + n_infectious_treatment
    #println("n_infectious_people: ", n_infectious_people)
    TTE = fn_TTE(nstatus, tstatus, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)
    (k, T, R0) = fn_update_R0(Gc, N, T_arr)
    #println("within function (samplesize) after 2nd update: ",samplesize)
    samplesize_mat = vcat(samplesize_mat, samplesize)
    println("samplesize_mat: ", samplesize_mat)
    n_infectious_people_mat = vcat(n_infectious_people_mat, n_infectious_people)
    println("n_infectious_people_mat: ", n_infectious_people_mat)
    TTE_mat = vcat(TTE_mat, TTE)
    println("TTE_mat: ", TTE_mat)
    VE_true_mat = vcat(VE_true_mat, VE_true)
    println("VE_true_mat: ", VE_true_mat)
    R0_mat = vcat(R0_mat, R0)
    println("R0_mat: ", R0_mat)
    #println("within function (samplesize_mat): ",samplesize_mat)
    #samplesize_mat = samplesize_mat[2:end] #Remove the first element from array
    #println("within function (samplesize) after 3rd update: ",samplesize)
    #println("within function (samplesize_mat): ",samplesize_mat)

    return samplesize_mat, n_infectious_people_mat, TTE_mat, VE_true_mat, R0_mat
end

# Set parameter values
N = 2000 # Total population size
begintime = 1.0
S0 = N # Number of suspectible people in the population at day 0
V0 = 0 # Number of vaccinated individuals in the population at day 0
casenum0 = 5 # Number of infectious individuals introduced into the community to begin the outbreak
immunenum0 = 0 # Number of people who are immune to the disease at the beginning of the outbreak
import_lambda = 1 # Number of occurrences variable for imported cases timeseries, assumed to follow Poisson Distribution
lambda0 = 0.000052 # Per-time-step probability of infection for a susceptible nodes from an infectious neighbour
prop_in_trial = 0.6 # Proportion of the population/ cluster will be enrolled in the trial
prop_in_highrisk = 0.2 # Proportion of the population/ in the cluster are at high risk
prop_in_hcw = (2.8/1000*N + N/8)/N # Based on data in England, 2.8 doctor per patient and 8 nurses per patient
prop_in_keyvec = 0.3 # Proportion of children in the population/ in the cluster who are key transmisison vector
## The households
# hhnum: Number of households in the network
# hhsize_avg: AVE_truerage size of one household
# hhsize_range: Range of household sizes (sizes are uniformly distributed)
par_hh = DataFrame(hhnum=650, hhsize_avg=3, hhsize_range=2)
#par_hh = DataFrame(hhnum=1, hhsize_avg=500, hhsize_range=0)

## The clustered network
# communitynum: Number of communities in the clustered network
# communitysize_avg: AVE_truerage size of one community
# communitysize_range: Range of community sizes (sizes are uniformly distributed)
par_community = DataFrame(communitynum=3, communitysize_avg=650, communitysize_range=10)
#par_community = DataFrame(communitynum=1, communitysize_avg=500, communitysize_range=0)

prop_in_eachcompartment = prop_in_highrisk
par_compartment = DataFrame(compartmentnum=2, compartmentsize_avg=(N*prop_in_eachcompartment)/2, compartmentsize_range=25)

## Contact and transmission probabilities between nodes
# cprob_hhwithin_cwithin: Probability of contacts of an edge between two nodes in the same household and the same community
# cprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# cprob_hhbetween_cbetween: Probability of transmission of an edge between two nodes in different communities
# tprob_hhwithin_cwithin: Probability of transmission of an edge between two nodes in the same household and the same community
# tprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# tprob_hhbetween_cbetween: Probability of contacts of an edge between two nodes in different communities
par_prob = DataFrame(cprob_hhwithin_cwithin=1, cprob_hhbetween_cwithin=1, cprob_hhbetween_cbetween=1, tprob_hhwithin_cwithin=lambda0, tprob_hhbetween_cwithin=lambda0, tprob_hhbetween_cbetween=lambda0)

## Disease properties
# Use Ebola-like parameters (from Hitchings (2018)) - Gamma-distributed
#par_disease = DataFrame(incubperiod_shape=3.11, incubperiod_rate=0.32, infectperiod_shape=1.13, infectperiod_rate=0.226)
# Use COVID-19-like parameters (from Li et al (2020)) - Gamma-distributed
#par_disease = DataFrame(incubperiod_shape=2.11, incubperiod_rate=0.4, infectperiod_shape=3.0, infectperiod_rate=0.35)
par_disease = DataFrame(incubperiod_shape=3.45, incubperiod_rate=0.66, infectperiod_shape=5.0, infectperiod_rate=0.8)

nsim = 25 # Number of simulations
trial_begintime = 10.0
endtime = trial_begintime + 190.0 # Duration of simulation (timestep is in a unit of days)
trial_endtime = endtime
trial_posthighrisk = 7.0 # The number of days after trial_begintime we begin to recruit the rest of the population
trial_communitynum = [1]
treatment_gp = 1 # it needs to be an integer, not an array
vac_efficacy = [0.6]
protection_threshold = 0.5
stage = 1 # Number of interim analyses will be done in the trial
allocation_ratio = [0.5 0.5]
#allocation_ratio = [0.9/(0.9+1) 1/(1+0.9)]
#allocation_ratio = [0.78/(0.0018+0.78) 0.0018/(0.0018+0.78)]
alpha = 0.9 # Desired type I error
power = 0.8
#delta = 0.5 # The clinical relevant different to power for
#sigma = 1 # The assumed value for the SD of patient responses
#two_sided = false # A logical variable indicating whether to use a two-sided null hypothesis
p_es = 0.65 # Effect size
p0_te = 0.55 # Treatment Effect
e = 1.64 # Efficacy stopping boundary
f = 1.80 # Futility stopping boundary

#Initialization
samplesize_mat = zeros(1)
n_infectious_people_mat = zeros(1)
TTE_mat = zeros(1)
VE_true_mat = zeros(1)
R0_mat = zeros(1)

gamma_infectperiod_duration = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),1000)
gamma_infectperiod_maxduration = maximum(gamma_infectperiod_duration)
gamma_incubperiod_duration = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),1000)
gamma_incubperiod_maxduration = maximum(gamma_incubperiod_duration)

sbm_sol_mat = zeros(Int, round(Int,endtime), 5) # Same as sbm_sol, except it is a matrix, not a DataFrame

for isim in 1:nsim
    # Initializations
    nstatus = fill('S', N, round(Int,endtime)) # health statuses
    nstatus = hcat([1:1:N;], nstatus) # Put indexes on the 1st column
    tstatus = fill(-1, N) # trial statuses: (-1=not in trial, 0=in control, 1=vaccine candidate number)
    tstatus = hcat([1:1:N;], tstatus) # Put indexes on the 1st column
    sbm_sol = DataFrame(S=fill(0,round(Int,endtime)), E=fill(0,round(Int,endtime)), I=fill(0,round(Int,endtime)), R=fill(0,round(Int,endtime)), V=fill(0,round(Int,endtime))) #Initialize the matrix which holds SEIR incidence of all timestep

    # Divide the population into compartments
    #(hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, compartmentnum_highrisksize_arr, compartmentnum_hcw_arr, compartmentnum_keyVE_truec_arr) = fn_cluster_partition(N, prop_in_highrisk, prop_in_hcw, prop_in_keyvec)

    # For high-risk group
    prop_in_eachcompartment = prop_in_highrisk
    par_compartment_highrisk = DataFrame(compartmentnum=2, compartmentsize_avg=(N*prop_in_eachcompartment)/2, compartmentsize_range=300)
    par_compartment = par_compartment_highrisk
    compartment_highrisksize_arr = fn_par_cluster(N, par_hh, par_community, par_compartment, "compartment")
    compartmentnum_highrisksize_arr = sample(1:N, round(Int,compartment_highrisksize_arr[1]), replace=false, ordered=true) # Assign compartment number to each individual in the population

    # For high-exposure group
    prop_in_eachcompartment = prop_in_hcw
    par_compartment_hcw = DataFrame(compartmentnum=2, compartmentsize_avg=(N*prop_in_eachcompartment)/2, compartmentsize_range=50)
    par_compartment = par_compartment_hcw
    compartmentsize_hcw_arr = fn_par_cluster(N, par_hh, par_community, par_compartment, "compartment")
    compartmentnum_hcw_arr = sample(1:N, round(Int,compartmentsize_hcw_arr[1]), replace=false, ordered=true) # Assign compartment number to each individual in the population

    # For key transmission VE_true
    prop_in_eachcompartment = prop_in_keyvec
    par_compartment_keyvec = DataFrame(compartmentnum=2, compartmentsize_avg=(N*prop_in_eachcompartment)/2, compartmentsize_range=200)
    par_compartment = par_compartment_keyvec
    compartmentsize_keyvec_arr = fn_par_cluster(N, par_hh, par_community, par_compartment, "compartment")
    compartmentnum_keyvec_arr = sample(1:N, round(Int,compartmentsize_keyvec_arr[1]), replace=false, ordered=true)  # Assign compartment number to each individual in the population

    # Compute the parameters of the clusters
    hhsize_arr = fn_par_cluster(N, par_hh, par_community, par_compartment, "household") # Define the sizes of each household
    hhsize_arr = round.(Int, hhsize_arr)
    hhnum_arr = fn_partition(hhsize_arr) # Assign household number to each individual in the population
    communitysize_arr = fn_par_cluster(N, par_hh, par_community, par_compartment, "community") # Define the sizes of each community
    communitysize_arr = round.(Int, communitysize_arr)
    communitynum_arr = fn_partition(communitysize_arr) # Assign community number to each individual in the population

    # Generate the number of imported cases for the duration of the outbreak
    importcasenum_timeseries = rand(Poisson(import_lambda), round(Int,endtime))
    importcasenum_timeseries[1] = casenum0

    nstatus_fn202 = nstatus
    tstatus_fn202 = tstatus
    sbm_sol_fn202 = sbm_sol

    # From importcasesnum_timeseries, obtain number of import cases at t=timestep
    for timestep202 in 1:round(Int,trial_endtime)
        casenum = importcasenum_timeseries[timestep202]

        # Generate node_names of imported cases at t=timestep
        s_elements = findall(x->x=='S', nstatus_fn202[:,timestep202+1])
        u_elements = findall(x->x=='U', nstatus_fn202[:,timestep202+1])
        elements = union(s_elements, u_elements)
        casenum_draw = minimum([length(elements), casenum])
        importcases = sample(elements, floor(Int,casenum_draw), replace=false) # Sampling without replacement

        # Compute the parameters for disease properties
        #incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
        #infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
        incubperiod = zeros(length(importcases)) # Incubation period of the disease
        infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),length(importcases)) # Infectious period of the disease

        for index1 in 1:(length(importcases))

            # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
            tbound1 = ceil(Int, min(timestep202 + ceil(incubperiod[index1,1]), round(Int,endtime)))
            tbound2 = ceil(Int, min(timestep202 + ceil(incubperiod[index1,1]) + ceil(infectperiod[index1,1]), round(Int,endtime)))

            # column_index start at 2 because nstatus_fn[:,1] is nodes_name
            # for index2 in timestep:(round(Int,tbound1))
            #    nstatus_fn[importcases[index1],index2+1] = "E"
            # end

            for index3 in tbound1:tbound2
                nstatus_fn202[importcases[index1],index3+1] = 'I'
            end

            for index4 in (tbound2+1):(round(Int,endtime))
                nstatus_fn202[importcases[index1],index4+1] = 'R'
            end
        end
    end

    # Compute who-contact-whom network graphs
    Gc = fn_contact_network(par_prob, hhsize_arr, communitysize_arr, hhnum_arr, communitynum_arr) # Construct a who-contact-whom stochastic block network

    timestep3 = 1
    nstatus_fn = fn_partialimmune(immunenum0, nstatus_fn202) # Generate immune people
    nstatus_fn = fn_importcases(par_disease, importcasenum_timeseries, nstatus_fn, timestep3) # Import infectious cases at t-timestep3

    D = fn_countelements(nstatus_fn[:,timestep3+1]) # Count number of occurrences of SEIRV at a particular t=timestep3

    # Initilizations
    T_arr = zeros(round(Int,endtime))
    V = zeros(Int8,V0) # V contains nodes_name of the vaccinated individuals

    Gc_fn202 = Gc
    sbm_sol = sbm_sol_fn202
    nstatus202 = nstatus_fn
    intermediatetime1 = begintime+1
    intermediatetime2 = endtime

    for timestep1 in (round(Int,intermediatetime1)):(round(Int,intermediatetime2))

        # Set local variables
        nstatus_fn = nstatus_fn202
        tstatus_fn = tstatus_fn202
        sbm_sol_fn = sbm_sol_fn202

        nstatus_fn = fn_importcases(par_disease, importcasenum_timeseries, nstatus_fn, timestep1) # Import cases
        Gt = fn_transmit_network(Gc_fn202, par_prob, hhnum_arr, communitynum_arr, nstatus_fn, timestep1) # Construct a who-infect-whom stochastic block network based on the contact network Gc
        potential_transmit_indexes = fn_findnonzeros(Gt) # The index numbers that will haVE_true disease transmission according to the stochastic block network model
        transmit_indexes = fn_uniqueS(potential_transmit_indexes, nstatus_fn, timestep1) # Check if potential_transmit_indexes are susceptibles
        T_arr[timestep1] = fn_computeT(N, par_prob, Gc_fn202, nstatus_fn, hhnum_arr, communitynum_arr, timestep1)

        if size(transmit_indexes,1)>0 # Check if there are infectees

            nstatus_fn = fn_spread(par_disease, nstatus_fn, transmit_indexes, tstatus_fn, timestep1) # Spread the diseae within the network and update nstatus

            # Count number of occurrences of SEIRV at a particular timestep
            for timestep2 in timestep1:(round(Int,intermediatetime2))

                D = fn_countelements(nstatus_fn[:,timestep2+1]) # Count number of occurrences of SEIRV at a particular t=timestep2

                # Put these SEIRV incidence values into a DataFrame sbm_sol_fn
                sbm_sol_fn[timestep2,:S] = D[1,:S] # Put S value into the appropriate cell
                sbm_sol_fn[timestep2,:E] = D[1,:E] # Put E value into the appropriate cell
                sbm_sol_fn[timestep2,:I] = D[1,:I] # Put I value into the appropriate cell
                sbm_sol_fn[timestep2,:R] = D[1,:R] # Put R value into the appropriate cell
                sbm_sol_fn[timestep2,:V] = D[1,:V] # Put V value into the appropriate cell
            end

            #sbm_sol = sbm_sol_fn
            nstatus = nstatus_fn
            tstatus = tstatus_fn
            Gc = Gc_fn202

            # Compute R0 in a network based on connectivity, without considering intervention
            if timestep1 == round(Int,intermediatetime2)
                (k,T,R0) = fn_update_R0(Gc, N, T_arr)
                println("k = ", k, ", T = ",T, ", and R0 = ",R0)
                global R0
            end
        end
    end

    # Same as sbm_sol, except it is a matrix, not a DataFrame
    sbm_sol_mat[:,1] = sbm_sol[:,1]
    sbm_sol_mat[:,2] = sbm_sol[:,2]
    sbm_sol_mat[:,3] = sbm_sol[:,3]
    sbm_sol_mat[:,4] = sbm_sol[:,4]
    sbm_sol_mat[:,5] = sbm_sol[:,5]

    # Set local variable (for this function)
    # tstatus_fn = tstatus

    # Find those who are not in trial
    potential_nodes_notrial_index = findall(x->x==-1, tstatus[:,2]) # Nodes name of those not in trial
    keyvec_arr = setdiff(potential_nodes_notrial_index, compartmentnum_keyvec_arr) # Nodes name of those key transmission vector - children - not in trial

    # Enrollment size
    enrolsize = ceil(Int,length(keyvec_arr)*prop_in_trial)
    #println("Enrollment size = ", enrolsize)

    if length(keyvec_arr)>0
        # Select those who will be enrolled into the trial
        nodes_in_trial = sample(keyvec_arr, enrolsize, replace=false, ordered=true)
        while length(unique(nodes_in_trial)) != enrolsize # Found a bug in sample. See duplicate elements despite replace=false.
            nodes_in_trial = sample(keyvec_arr, enrolsize, replace=false, ordered=true)
        end
        n_treatment_enroll = floor(Int, enrolsize*allocation_ratio[2])
        elements_draw = minimum([length(nodes_in_trial), n_treatment_enroll])
        nodes_in_treatment = sample(nodes_in_trial, elements_draw, replace=false, ordered=true)
        nodes_in_control = setdiff(nodes_in_trial, nodes_in_treatment)
        #println("Individuals in trial = ", length(nodes_in_trial))
        println("Number of individuals in treatment = ", length(nodes_in_treatment))
        #println("Individuals in control = ", length(nodes_in_control))

        # Assign trial statuses and update tstatus_fn
        for index3 in 1:length(nodes_in_treatment)
            tstatus[nodes_in_treatment[index3],2] = 1 # In treatment group
        end

        for index4 in 1:length(nodes_in_control)
            tstatus[nodes_in_control[index4],2] = 0 # In control group
        end
    end

    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE = fn_TTE(nstatus, tstatus, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)


    println("Results:")
    println("Number of iteration: ", isim)
    if isim==1
        samplesize_mat = deepcopy(samplesize)
        n_infectious_people = n_infectious_control + n_infectious_treatment
        n_infectious_people_mat = deepcopy(n_infectious_people)
        TTE_mat = deepcopy(TTE)
        VE_true_mat = deepcopy(VE_true)
        R0_mat = deepcopy(R0)
        global samplesize_mat, n_infectious_people_mat, TTE_mat, VE_true_mat, R0_mat

    else
        (samplesize_mat, n_infectious_people_mat, TTE_mat, VE_true_mat, R0_mat) = fn_update_vars(samplesize, samplesize_mat, n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_infectious_people_mat, TTE_mat, VE_true_mat, R0_mat, treatment_gp, timestep_fn, alpha, power, nstatus, tstatus, trial_begintime, trial_endtime, gamma_infectperiod_maxduration, Gc, N, T_arr)

        if isim==nsim
            samplesize_CI = quantile(samplesize_mat, [0.5, 0.05, 0.95])
            println("Sample size (from true cases): mean: ", mean(filter(isfinite, samplesize_mat)), ", CI: ", samplesize_CI)

            n_infectious_people_CI = quantile(n_infectious_people_mat, [0.5, 0.05, 0.95])
            println("Average number of infectious people in the trial: mean: ", mean(filter(isfinite, n_infectious_people_mat)), ", CI: ", n_infectious_people_CI)

            TTE_CI = quantile(TTE_mat[:,2], [0.5, 0.05, 0.95])
            println("Time-to-Event: mean: ", mean(filter(isfinite, TTE_mat)), ", CI: ", TTE_CI)

            VE_CI = quantile(VE_true_mat, [0.5, 0.05, 0.95])
            println("Vaccine efficacy: mean: ", mean(filter(isfinite, VE_true_mat)), ", CI: ", VE_CI)

            R0_CI = quantile(R0_mat, [0.5, 0.05, 0.95])
            println("Reproductive number without intervention: mean: ", mean(filter(isfinite, R0_mat)), ", CI: ", R0_CI)

            #Y = fn_divide(soln_mat3, endtime, nsim, 1)
            #lowerCI = zeros(round(Int,endtime))
            #upperCI = zeros(round(Int,endtime))
            #medianCI = zeros(round(Int,endtime))
            #for index = 1:round(Int,endtime)
            #    lowerCI[index] = quantile!(Y[index,:], 0.05)
            #    upperCI[index] = quantile!(Y[index,:], 0.95)
            #    medianCI[index] = median(Y[index,:])
            #end
            #X = [lowerCI, medianCI, upperCI]
            #plot1 = plot(1:Int(endtime), Y, legend=false)
            #plot2 = plot(1:Int(endtime), X, legend=[lowerCI, median, upperCI])
            #plot(plot1, plot2, layout = (2, 1))
        end
    end
end
