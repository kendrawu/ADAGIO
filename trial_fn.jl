# Function to generate imported cases and randomly distribute them into different clusters
function fn_importcases(par_disease, importcasenum_timeseries, nstatus, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # importcasenum_timeseries: Array that holds number of infectious individuals introduced into the community throughout the duration of the outbreak
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # nstatus_fn: Health statuses of all individuals in the population at each time step

    # Create a local variable for nstatus
    nstatus_fn = nstatus

    # From importcasesnum_timeseries, obtain number of import cases at t=timestep
    casenum = importcasenum_timeseries[timestep]

    # Generate node_names of imported cases at t=timestep
    s_elements = findall(x->x=='S', nstatus_fn[:,timestep+1])
    u_elements = findall(x->x=='U', nstatus_fn[:,timestep+1])
    elements = union(s_elements, u_elements)
    importcases = sample(elements, casenum, replace=false) # Sampling without replacement

    # Compute the parameters for disease properties
    #incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
    #infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
    incubperiod = zeros(length(importcases)) # Incubation period of the disease
    infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),length(importcases)) # Infectious period of the disease

    for index1 in 1:(length(importcases))

        # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
        tbound1 = ceil(Int, min(timestep + ceil(incubperiod[index1,1]), round(Int,endtime)))
        tbound2 = ceil(Int, min(timestep + ceil(incubperiod[index1,1]) + ceil(infectperiod[index1,1]), round(Int,endtime)))

        # column_index start at 2 because nstatus_fn[:,1] is nodes_name
        # for index2 in timestep:(round(Int,tbound1))
        #    nstatus_fn[importcases[index1],index2+1] = "E"
        # end

        for index3 in tbound1:tbound2
            nstatus_fn[importcases[index1],index3+1] = 'I'
        end

        for index4 in (tbound2+1):(round(Int,endtime))
            nstatus_fn[importcases[index1],index4+1] = 'R'
        end

    end

    return nstatus_fn
end

# Function to simulate the spread of the disease and return the statuses of each nodes_name at all timesteps
function fn_spread(par_disease, nstatus, infectees, tstatus, timestep)

   # Inputs:
   # par_disease: User-defined parameters of the disease
   # nstate: Health statuses of all the individuals in the population at each time step
   # infectees: A list of nodes_names that are to be infected at t=timestep according to SBM
   # tstatus: The trial status of nodes_names at timestep
   # timestep: The time t of the outbreak

   # Output:
   # nstatus_fn: Health statuses of all the individuals in the population at each time step

   # Create a local variables
   nstatus_fn = nstatus

   # Compute the parameters for disease properties
   incubperiod = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),length(infectees)) # Incubation period of the disease
   infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),length(infectees)) # Infectious period of the disease

   # From the infected cases in infectees, allow health statuses change as time progresses
   for index1 in 1:(length(infectees))

       # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
       tbound1 = ceil(Int, min(timestep + incubperiod[index1,1], round(Int,endtime)))
       tbound2 = ceil(Int, min(timestep + incubperiod[index1,1] + infectperiod[index1,1], round(Int,endtime)))

       # column_index start at 2 because nstatus[:,1] is nodes_name
       for index2 in timestep:tbound1
           nstatus_fn[infectees[index1],index2+1] = 'E'
       end

       for index3 in (tbound1+1):tbound2
           nstatus_fn[infectees[index1],index3+1] = 'I'
       end

       for index4 in (tbound2+1):(round(Int,endtime))
           nstatus_fn[infectees[index1],index4+1] = 'R'
       end
   end


   # Obtain node names that are in treatment groups
   V = findall(x->x>0,tstatus[:,2])
   if length(V)>0
       for index5 in 1:length(V)
           v_element = V[index5]
           if nstatus_fn[v_element,timestep+1] == 'S' # Check if the individual has been infected
               for t in timestep:(round(Int,endtime))
                   nstatus_fn[v_element,t+1] = 'V'
               end

               # Determine which treatment group has vaccine efficacy below threshold
               treatmentgp_unprotect = findall(x->x<protection_threshold, vac_efficacy[1,:])

               if length(treatmentgp_unprotect)>0
                   for index6 in 1:length(treatmentgp_unprotect)
                       vnodes_unprotected = findall(x->x==index6, tstatus_fn[:,2]) # Check who is in this unprotected treatment group
                       for index7 in 1:length(vnodes_unprotected)
                           v_unprotected_element = vnodes_unprotected[index7]
                           if nstatus_fn[v_unprotected_element,timestep+1] == 'V'
                               for t in timestep:(round(Int,endtime))
                                   nstatus_fn[v_unprotected_element,t+1] = 'U'
                               end
                           end
                       end
                   end
               end
           end
       end
   end

   return nstatus_fn
end

# From nonzeros_indexes, this function finds and return the node_names of those who are susceptible at t=(timestep+1)
function fn_uniqueS(nonzeros_indexes, nstatus, timestep)

    # Inputs:
    # nonzeros_indexes: A list of indexes from the SBM
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # unique_indexes: A list of nodes_names that can potentially be infected at t=(timestep+1) according to SBM

    # Initialization
    s_elements = Int[]
    sizehint!(s_elements, size(nstatus,1))
    u_elements = Int[]
    sizehint!(u_elements, size(nstatus,1))

    if (length(nonzeros_indexes))>0
        for i in 1:(length(nonzeros_indexes))
            # Obtain nodes_name if nonzeros_indexes[i]'s status is susceptible at time=(timestep+1)
            if nstatus[nonzeros_indexes[i],timestep+1] == 'S'
                push!(s_elements::Array{Int,1},nstatus[nonzeros_indexes[i],1])
            elseif nstatus[nonzeros_indexes[i],timestep+1] == 'U'
                push!(u_elements::Array{Int,1},nstatus[nonzeros_indexes[i],1])
            end
        end

        update_s_elements = union(s_elements, u_elements)

        # Take into account susceptible deplection to find the minimum between the potential infectees and the available susceptibles
        bound = min(length(nonzeros_indexes), length(update_s_elements)) # Find the minimum among the variables

        if bound> 0
            unique_indexes = sample(update_s_elements, bound, replace=false)
            return unique_indexes
        else
            unique_indexes = Int8[]
            return unique_indexes
        end
    else
        unique_indexes = Int8[]
        return unique_indexes
    end
end

function fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)

    # Inputs:
    # N: Total population size
    # par_hh: User-defined parameters of households
    # par_community: User-defined parameters of communities
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # par_disease: User-defined parameters of the disease
    # import_lambda: Number of occurrences variable for imported cases timeseries, assumed to follow Poisson Distribution
    # casenum0: Number of infectious individuals introduced into the community to begin the outbreak
    # immunenum0: Number of people who are immune to the disease since day 0 of the outbreak
    # endtime: The time to end the simulation of the entire outbreak

    # Outputs:
    # nstatus: Health statuses of all the individuals in the population at each time step
    # tstatus: The trial statuses of each individual in the population
    # sbm_sol: The incidences of SEIRV of all timesteps
    # hhsize_arr: Sizes of each household
    # communitysize_arr: Sizes of each community
    # hhnum_arr: Household number of each individuals in the population
    # communitynum_arr: Community number of each individuals in the population
    # importcasenum_timeseries: node_names of imported cases on a time series for the duration of the outbreak
    # Gc: The who-contact-whom network map

    # Initializations
    nstatus = fill('S', N, round(Int,endtime)) # health statuses
    nstatus = hcat([1:1:N;], nstatus) # Put indexes on the 1st column
    tstatus = fill(-1, N) # trial statuses: (-1=not in trial, 0=in control, 1=vaccine candidate number)
    tstatus = hcat([1:1:N;], tstatus) # Put indexes on the 1st column
    sbm_sol = DataFrame(S=fill(0,round(Int,endtime)), E=fill(0,round(Int,endtime)), I=fill(0,round(Int,endtime)), R=fill(0,round(Int,endtime)), V=fill(0,round(Int,endtime))) #Initialize the matrix which holds SEIR incidence of all timestep

    # Compute the parameters of the clusters
    hhsize_arr = fn_par_cluster(N, par_hh, par_community, "household") # Define the sizes of each household
    hhnum_arr = fn_partition(hhsize_arr) # Assign household number to each individual in the population
    communitysize_arr = fn_par_cluster(N, par_hh, par_community, "community") # Define the sizes of each community
    communitynum_arr = fn_partition(communitysize_arr) # Assign community number to each individual in the population

    # Generate the number of imported cases for the duration of the outbreak
    importcasenum_timeseries = fn_importcases_timeseries(import_lambda, casenum0, endtime)

    # Compute who-contact-whom network graphs
    Gc = fn_contact_network(par_prob, hhsize_arr, communitysize_arr, hhnum_arr, communitynum_arr) # Construct a who-contact-whom stochastic block network

    timestep3 = 1
    nstatus = fn_partialimmune(immunenum0, nstatus) # Generate immune people
    nstatus = fn_importcases(par_disease, importcasenum_timeseries, nstatus, timestep3) # Import infectious cases at t-timestep3

    D = fn_countelements(nstatus[:,timestep3+1]) # Count number of occurrences of SEIRV at a particular t=timestep3

    # Put these SEIRV incidence values into a DataFrame sbm_sol
    sbm_sol[timestep3,:S] = D[1,:S] # Put S value into the appropriate cell
    sbm_sol[timestep3,:E] = D[1,:E] # Put E value into the appropriate cell
    sbm_sol[timestep3,:I] = D[1,:I] # Put I value into the appropriate cell
    sbm_sol[timestep3,:R] = D[1,:R] # Put R value into the appropriate cell
    sbm_sol[timestep3,:V] = D[1,:V] # Put V value into the appropriate cell

    return nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc
end

function fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, intermediatetime1, intermediatetime2, endtime)

    # Inputs:
    # nstatus: Health statuses of all the individuals in the population at each time step
    # tstatus: The trial statuses of each individual in the population
    # sbm_sol: The incidences of SEIRV of all timesteps
    # par_hh: User-defined parameters of the households
    # par_community: User-defined parameters of the network clusters
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # par_disease: User-defined parameters of the disease
    # hhsize_arr: Sizes of each household
    # communitysize_arr: Sizes of each community
    # hhnum_arr: Household number of each individuals in the population
    # communitynum_arr: Community number of each individuals in the population
    # importcasenum_timeseries: node_names of imported cases on a time series for the duration of the outbreak
    # Gc: The who-contact-whom network map
    # intermediatetime1: The immediate start time of a segment of the outbreak
    # intermediatetime2: The immediate end time of a segment of the outbreak
    # endtime: The end simulation of the entire outbreak

    # Outputs:
    # nstatus: Health statuses of all the individuals in the population at each time step
    # tstatus: The trial statuses of each individual in the population
    # sbm_sol: The incidences of SEIRV of all timesteps
    # T: Transmissibility
    # R0: Basic reproductive number

    # Initializations
    T_arr = zeros(round(Int,endtime))
    V = zeros(Int8,V0) # V contains nodes_name of the vaccinated individuals

    # Begin transmission
    N = size(nstatus)[1]

    for timestep1 in (round(Int,intermediatetime1)):(round(Int,intermediatetime2))

        # Set local variables
        nstatus_fn = nstatus
        tstatus_fn = tstatus
        sbm_sol_fn = sbm_sol
        Gc_fn = Gc

        nstatus_fn = fn_importcases(par_disease, importcasenum_timeseries, nstatus_fn, timestep1) # Import cases
        Gt = fn_transmit_network(Gc_fn, par_prob, hhnum_arr, communitynum_arr, nstatus_fn, timestep1) # Construct a who-infect-whom stochastic block network based on the contact network Gc
        potential_transmit_indexes = fn_findnonzeros(Gt) # The index numbers that will have disease transmission according to the stochastic block network model
        transmit_indexes = fn_uniqueS(potential_transmit_indexes, nstatus_fn, timestep1) # Check if potential_transmit_indexes are susceptibles
        T_arr[timestep1] = fn_computeT(N, par_prob, Gc_fn, nstatus_fn, hhnum_arr, communitynum_arr, timestep1)

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
        end
        sbm_sol = sbm_sol_fn
        nstatus = nstatus_fn
        tstatus = tstatus_fn
        Gc = Gc_fn

        # Compute R0 in a network
        if timestep1 == round(Int,intermediatetime2)
            k = sum(sum(Gc))/N # mean degree of the network
            T = mean(T_arr[T_arr .> 0]) # Average probability that an infectious individual will transmit the disease to a susceptible individual with whom they have contact
            R0 = T * (k^2/k - 1)
            println("k = ", k, ", T = ",T, ", and R0 = ",R0)

            global T
            global R0
        end
    end
    return nstatus, tstatus, sbm_sol, T, R0
end

# Function to enroll the population into trial by individuals
function fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)

    # Inputs:
    # N: Total size of population
    # tstatus: The trial statuses of each individual in the population
    # prop_in_trial: Proportion of a cluster will be enrolled into the trial
    # allocation_ratio: Allocation ratio of the randomization trial
    # vac_efficacy: Vaccine efficacy of vaccine candidate(s)
    # protection_threshold: A value to determine if a vaccinated person will be infected if exposed

    # Output:
    # tstatus_fn: The trial statuses of each individual in the population

    # Set local variable for this function
    tstatus_fn = tstatus

    # Find those who are not in trial
    potential_nodes_notrial_index = findall(x->x==-1, tstatus_fn[:,2]) # Nodes name of those not in trial

    # Enrollment size
    enrolsize = ceil(Int,length(potential_nodes_notrial_index)*prop_in_trial)
    #println("Enrollment size = ", enrolsize)

    if length(potential_nodes_notrial_index)>0
        # Select those who will be enrolled into the trial
        nodes_in_trial = sample(potential_nodes_notrial_index, enrolsize, replace=false, ordered=true)
        while length(unique(nodes_in_trial)) != enrolsize # Found a bug in sample. See duplicate elements despite replace=false.
            nodes_in_trial = sample(potential_nodes_notrial_index, enrolsize, replace=false, ordered=true)
        end
        nodes_in_treatment = sample(nodes_in_trial, ceil(Int,length(potential_nodes_notrial_index)*prop_in_trial*allocation_ratio[2]), replace=false, ordered=true)
        nodes_in_control = setdiff(nodes_in_trial, nodes_in_treatment)
        #println("Individuals in trial = ", length(nodes_in_trial))
        #println("Individuals in treatment = ", length(nodes_in_treatment))
        #println("Individuals in control = ", length(nodes_in_control))

        # Assign trial statuses and update tstatus_fn
        for index3 in 1:length(nodes_in_treatment)
            tstatus_fn[nodes_in_treatment[index3],2] = 1 # In treatment group
        end
        for index4 in 1:length(nodes_in_control)
            tstatus_fn[nodes_in_control[index4],2] = 0 # In control group
        end
    end

    return tstatus_fn, nodes_in_control, nodes_in_treatment
end

# Function to generate imported cases and distribute them into a specific clusters
function fn_importcases_cRCT(par_disease, importcasenum_timeseries, nstatus, trial_communitynum, communitynum_arr, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # importcasenum_timeseries: Array that holds number of infectious individuals introduced into the community throughout the duration of the outbreak
    # nstatus: Health statuses of all individuals in the population at each time step
    # trial_communitynum: The communitynum that will be recruited
    # communitynum_arr: Holds the community number of each individuals in the population
    # timestep: The time t of the outbreak

    # Output:
    # nstatus_fn: Health statuses of all individuals in the population at each time step

    # Create a local variable for nstatus
    nstatus_fn = nstatus

    # From importcasesnum_timeseries, obtain number of import cases at t=timestep
    casenum = importcasenum_timeseries[round(Int,timestep)]

    # Generate node_names of imported cases at t=timestep
    s_elements = findall(x->x=='S', nstatus_fn[:,round(Int,timestep)+1])
    communitynum_elements = findall(x->x==trial_communitynum[1], communitynum_arr) # for now
    s_elements_communitynum = intersect(s_elements, communitynum_elements)
    importcases = sample(s_elements_communitynum, casenum, replace=false) # Sampling without replacement

    # Compute the parameters for disease properties
    #incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
    #infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
    incubperiod = zeros(length(importcases)) # Incubation period of the disease
    infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),length(importcases)) # Infectious period of the disease

    for index1 in 1:(length(importcases))

        # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
        tbound1 = ceil(Int, min(round(Int,timestep) + ceil(incubperiod[index1,1]), round(Int,endtime)))
        tbound2 = ceil(Int, min(round(Int,timestep) + ceil(incubperiod[index1,1]) + ceil(infectperiod[index1,1]), round(Int,endtime)))

        # column_index start at 2 because nstatus_fn[:,1] is nodes_name
        # for index2 in timestep:(round(Int,tbound1))
        #    nstatus_fn[importcases[index1],index2+1] = "E"
        # end

        for index3 in tbound1:tbound2
            nstatus_fn[importcases[index1],index3+1] = 'I'
        end

        for index4 in (tbound2+1):(round(Int,endtime))
            nstatus_fn[importcases[index1],index4+1] = 'R'
        end

    end

    return nstatus_fn
end

function fn_transmodel_cRCT(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, trial_communitynum, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, intermediatetime1, intermediatetime2, endtime)

    # Inputs:
    # nstatus: Health statuses of all the individuals in the population at each time step
    # tstatus: The trial statuses of each individual in the population
    # sbm_sol: The incidences of SEIRV of all timesteps
    # par_hh: User-defined parameters of the households
    # par_community: User-defined parameters of the network clusters
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # par_disease: User-defined parameters of the disease
    # hhsize_arr: Sizes of each household
    # communitysize_arr: Sizes of each community
    # hhnum_arr: Household number of each individuals in the population
    # trial_communitynum: Communitynum that will be recruited
    # communitynum_arr: Community number of each individuals in the population
    # importcasenum_timeseries: node_names of imported cases on a time series for the duration of the outbreak
    # Gc: The who-contact-whom network map
    # intermediatetime1: The immediate start time of a segment of the outbreak
    # intermediatetime2: The immediate end time of a segment of the outbreak
    # endtime: The end simulation of the entire outbreak

    # Outputs:
    # nstatus: Health statuses of all the individuals in the population at each time step
    # tstatus: The trial statuses of each individual in the population
    # sbm_sol: The incidences of SEIRV of all timesteps
    # T: Transmissibility
    # R0: Basic reproductive number

    # Initializations
    T_arr = zeros(round(Int,endtime))
    V = zeros(Int8,V0) # V contains nodes_name of the vaccinated individuals

    # Begin transmission
    N = size(nstatus)[1]

    for timestep1 in (round(Int,intermediatetime1)):(round(Int,intermediatetime2))

        # Set local variables
        nstatus_fn = nstatus
        tstatus_fn = tstatus
        sbm_sol_fn = sbm_sol
        Gc_fn = Gc

        nstatus_fn = fn_importcases_cRCT(par_disease, importcasenum_timeseries, nstatus_fn, trial_communitynum, communitynum_arr, timestep1) # Import cases
        Gt = fn_transmit_network(Gc_fn, par_prob, hhnum_arr, communitynum_arr, nstatus_fn, timestep1) # Construct a who-infect-whom stochastic block network based on the contact network Gc
        potential_transmit_indexes = fn_findnonzeros(Gt) # The index numbers that will have disease transmission according to the stochastic block network model
        transmit_indexes = fn_uniqueS(potential_transmit_indexes, nstatus_fn, timestep1) # Check if potential_transmit_indexes are susceptibles
        T_arr[timestep1] = fn_computeT(N, par_prob, Gc_fn, nstatus_fn, hhnum_arr, communitynum_arr, timestep1)

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
        end
        sbm_sol = sbm_sol_fn
        nstatus = nstatus_fn
        tstatus = tstatus_fn
        Gc = Gc_fn

        # Compute R0 in a network
        if timestep1 == round(Int,intermediatetime2)
            k = sum(sum(Gc))/N # mean degree of the network
            T = mean(T_arr[T_arr .> 0]) # Average probability that an infectious individual will transmit the disease to a susceptible individual with whom they have contact
            R0 = T * (k^2/k - 1)
            println("k = ", k, ", T = ",T, ", and R0 = ",R0)

            global T
            global R0
        end
    end
    return nstatus, tstatus, sbm_sol, T, R0
end

# Function to enroll the population into trial by community cluster(s)
function fn_trialsetup_cRCT(N, par_disease, tstatus, communitysize_arr, communitynum_arr, trial_communitynum, nstatus, allocation_ratio, vac_efficacy, protection_threshold)

    # Inputs:
    # N: Total size of population
    # par_disease: User-defined parameters of the disease
    # tstatus: The trial statuses of each individual in the population
    # communitysize_arr: Sizes of each community
    # communitynum_arr: Holds the community number that each individual in the population belongs
    # trial_communitynum: The cluster numbers that are in the trial
    # nstatus: The health statuses of each individual in the population
    # prop_in_trial: Proportion of a cluster will be enrolled into the trial
    # allocation_ratio: Allocation ratio of the randomization trial
    # vac_efficacy: Vaccine efficacy of vaccine candidate(s)
    # protection_threshold: A value to determine if a vaccinated person will be infected if exposed

    # Output:
    # tstatus_fn: The trial statuses of each individual in the population

    # Initialization
    potential_nodes_in_trial_index = Int8[]
    sizehint!(potential_nodes_in_trial_index, N)
    potential_nodes_in_trial = Int16[]
    sizehint!(potential_nodes_in_trial, N)

    # Set local variable for this function
    tstatus_fn = tstatus

    # Combine tstatus and communitynum_arr for info searching
    tstatus_fn = hcat(tstatus, communitynum_arr)

    for index1 in 1:length(trial_communitynum)
        communitysize = communitysize_arr[trial_communitynum[index1]] # Obtain size of community selected for trial
        enrolsize = ceil(Int,communitysize*prop_in_trial) # Size of enrollment in the community
        #println("Enrollment size = ", enrolsize)
        #println("index1 = ", index1)
        potential_nodes_in_trial_index = findall(x->x==index1, tstatus_fn[:,3]) # Find everyone who are in this community

        if (size(potential_nodes_in_trial_index)[1])>0
            potential_nodes_in_trial = Int16[]
            for index2 in 1:length(potential_nodes_in_trial_index)
                potential_nodes_in_trial_index_local = potential_nodes_in_trial_index
                push!(potential_nodes_in_trial::Array{Int16,1},potential_nodes_in_trial_index_local[CartesianIndex(index2,1)][1]) # Extract node names from potential_nodes_in_trial_index
            end
        end

        # Select those who will be enrolled into the trial (assume two arms, for now)
        nodes_in_trial = sample(potential_nodes_in_trial, enrolsize, replace=false, ordered=true)
        while length(unique(nodes_in_trial)) != enrolsize # Found a bug in sample. See duplicate elements despite replace=false.
            nodes_in_trial = sample(potential_nodes_in_trial, enrolsize, replace=false, ordered=true)
        end
        nodes_in_treatment = sample(nodes_in_trial, ceil(Int,communitysize*prop_in_trial*allocation_ratio[2]), replace=false, ordered=true)
        nodes_in_control = setdiff(nodes_in_trial, nodes_in_treatment)
        #println("Individuals in trial = ", length(nodes_in_trial))
        #println("Individuals in treatment = ", length(nodes_in_treatment))
        #println("Individuals in control = ", length(nodes_in_control))

        # Assign trial statuses and update tstatus_fn
        for index3 in 1:length(nodes_in_treatment)
            tstatus_fn[nodes_in_treatment[index3],2] = 1 # In treatment group
        end
        if length(nodes_in_control)>0
            for index4 in 1:length(nodes_in_control)
                tstatus_fn[nodes_in_control[index4],2] = 0 # In control group
            end
        end
    end

    tstatus_fn = tstatus_fn[:,1:2] # Only keep the first 2 columns

    return tstatus_fn
end

function fn_iternation_cRCT_non_adpative(nsim, soln1, tstatus1, VE_true1, samplesize1, N, par_hh, par_community, par_prob, par_disease, prop_in_trial, import_lambda, casenum0, immunenum0, trial_communitynum, allocation_ratio, vac_efficacy, protection_threshold, treatment_gp, gamma_infectperiod_maxduration, trial_begintime, trial_endtime, endtime)

    soln1_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln1_mat[:,1] = soln1[:,1]
    soln1_mat[:,2] = soln1[:,2]
    soln1_mat[:,3] = soln1[:,3]
    soln1_mat[:,4] = soln1[:,4]
    soln1_mat[:,5] = soln1[:,5]

    (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel_cRCT(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, trial_communitynum, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
    tstatus = fn_trialsetup_cRCT(N, par_disease, tstatus, communitysize_arr, communitynum_arr, trial_communitynum, nstatus, allocation_ratio, vac_efficacy, protection_threshold)
    (nstatus2, tstatus2, soln2, T2, R02) = fn_transmodel_cRCT(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, trial_communitynum, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true2) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize2 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE2 = fn_TTE(nstatus2, tstatus2, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    soln2_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln2_mat[:,1] = soln2[:,1]
    soln2_mat[:,2] = soln2[:,2]
    soln2_mat[:,3] = soln2[:,3]
    soln2_mat[:,4] = soln2[:,4]
    soln2_mat[:,5] = soln2[:,5]
    soln_mat = vcat(soln1_mat, soln2_mat)
    nstatus_mat = vcat(nstatus1, nstatus2)
    tstatus_mat = vcat(tstatus1, tstatus2)
    VE_true_mat = vcat(VE_true1, VE_true2)
    samplesize_mat = vcat(samplesize1, samplesize2)
    TTE_mat = vcat(TTE1, TTE2)
    T = vcat(T1, T2)
    R0 = vcat(R01, R02)

    for itern = 3:nsim

        (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
        (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel_cRCT(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, trial_communitynum, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
        tstatus = fn_trialsetup_cRCT(N, par_disease, tstatus, communitysize_arr, communitynum_arr, trial_communitynum, nstatus, allocation_ratio, vac_efficacy, protection_threshold)
        (nstatus3, tstatus3, soln3, T3, R03) = fn_transmodel_cRCT(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, trial_communitynum, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

        timestep_fn = endtime
        (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true3) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
        samplesize3 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
        TTE3 = fn_TTE(nstatus3, tstatus3, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

        soln3_mat = zeros(Int, round(Int,endtime), 5)
        soln3_mat[:,1] = soln3[:,1]
        soln3_mat[:,2] = soln3[:,2]
        soln3_mat[:,3] = soln3[:,3]
        soln3_mat[:,4] = soln3[:,4]
        soln3_mat[:,5] = soln3[:,5]
        soln_mat = vcat(soln_mat, soln3_mat)
        nstatus_mat = vcat(nstatus_mat, nstatus3)
        tstatus_mat = vcat(tstatus_mat, tstatus3)
        VE_true_mat = vcat(VE_true_mat, VE_true3)
        samplesize_mat = vcat(samplesize_mat, samplesize3)
        TTE_mat = vcat(TTE_mat, TTE3)
        T = vcat(T, T3)
        R0 = vcat(R0, R03)
    end

    return soln_mat, nstatus_mat, tstatus_mat, VE_true_mat, samplesize_mat, TTE_mat, communitysize_arr, communitynum_arr, T, R0
end

# Function to enroll the contact list of a infector into the trial
function fn_trialsetup_contacts(Gc, nstatus, infector_node_name, communitynum_arr, timestep)

    # Inputs:
    # Gc: The who-contact-whom stochastic block matrix graph
    # nstatus: Health statuses of all individuals in the population at each time step
    # infector_node_name: Holds the node names of an infector
    # communitynum_arr: Holds the community number that each individual in the population belongs
    # timestep: The time t of outbreak

    # Output:
    # tstatus_fn: The trial statuses of each individual in the population

    # Set local variable for this function
    tstatus_fn = tstatus

    # Combine tstatus and communitynum_arr for info searching
    tstatus_fn = hcat(tstatus, communitynum_arr)

    # Initialization
    nodes_in_trial = zeros(Int8, N)

    # Find the contact list of infector_node_name
    if nstatus[infector_node_name, timestep+1] != 'I'
        println("wrong infector_node_name")
    else
        nodes_in_trial = findall(x->x!=0, Gc[infector_node_name,:])
    end

    # Select those who will be enrolled into the trial (assume two arms, for now)
    nodes_in_treatment = sample(nodes_in_trial, ceil(Int,length(nodes_in_trial)*prop_in_trial*allocation_ratio[2]), replace=false, ordered=true)
    nodes_in_control = setdiff(nodes_in_trial, nodes_in_treatment)
    println("nodes in trial = ", length(nodes_in_trial) , nodes_in_trial)
    println("nodes in treatment = ", length(nodes_in_treatment) , nodes_in_treatment)
    println("nodes in control = ", length(nodes_in_control) , nodes_in_control)

    # Assign trial statuses and update tstatus_fn
    for index3 in 1:length(nodes_in_treatment)
        tstatus_fn[nodes_in_treatment[index3],2] = 1 # In treatment group
    end
    for index4 in 1:length(nodes_in_control)
        tstatus_fn[nodes_in_control[index4],2] = 0 # In control group
    end

    return tstatus_fn
end

function fn_nsim1_iternation(nsim, soln1, N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)

    soln = zeros(Int, nsim*round(Int,endtime), 5)
    (soln2, nstatus2, communitysize_arr2, communitynum_arr2, T2, R02) = fn_transmodel_nsim1(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    soln = vcat(soln1, soln2)
    T = vcat(T1, T2)
    R0 = vcat(R01, R02)

    for itern = 3:nsim
        (soln2, nstatus2, communitysize_arr2, communitynum_arr2, T2, R02) = fn_transmodel_nsim1(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
        soln = vcat(soln, soln2)
        T = vcat(T, T2)
        R0 = vcat(R0, R02)
    end

    return soln, nstatus2, communitysize_arr2, communitynum_arr2, T, R0
end

function fn_iternation_iRCT_non_adpative(nsim, soln1, tstatus1, VE_true1, samplesize1, N, par_hh, par_community, par_prob, par_disease, prop_in_trial, import_lambda, casenum0, immunenum0, allocation_ratio, vac_efficacy, protection_threshold, treatment_gp, gamma_infectperiod_maxduration, trial_begintime, trial_endtime, endtime)

    soln1_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln1_mat[:,1] = soln1[:,1]
    soln1_mat[:,2] = soln1[:,2]
    soln1_mat[:,3] = soln1[:,3]
    soln1_mat[:,4] = soln1[:,4]
    soln1_mat[:,5] = soln1[:,5]

    (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
    (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)
    (nstatus2, tstatus2, soln2, T2, R02) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true2) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize2 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE2 = fn_TTE(nstatus2, tstatus2, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    soln2_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln2_mat[:,1] = soln2[:,1]
    soln2_mat[:,2] = soln2[:,2]
    soln2_mat[:,3] = soln2[:,3]
    soln2_mat[:,4] = soln2[:,4]
    soln2_mat[:,5] = soln2[:,5]
    soln_mat = vcat(soln1_mat, soln2_mat)
    nstatus_mat = vcat(nstatus1, nstatus2)
    tstatus_mat = vcat(tstatus1, tstatus2)
    VE_true_mat = vcat(VE_true1, VE_true2)
    samplesize_mat = vcat(samplesize1, samplesize2)
    TTE_mat = vcat(TTE1, TTE2)
    T = vcat(T1, T2)
    R0 = vcat(R01, R02)

    for itern = 3:nsim
        (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
        (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
        (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)
        (nstatus3, tstatus3, soln3, T3, R03) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

        timestep_fn = endtime
        (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true3) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
        samplesize3 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
        TTE3 = fn_TTE(nstatus3, tstatus3, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

        soln3_mat = zeros(Int, round(Int,endtime), 5)
        soln3_mat[:,1] = soln3[:,1]
        soln3_mat[:,2] = soln3[:,2]
        soln3_mat[:,3] = soln3[:,3]
        soln3_mat[:,4] = soln3[:,4]
        soln3_mat[:,5] = soln3[:,5]
        soln_mat = vcat(soln_mat, soln3_mat)
        nstatus_mat = vcat(nstatus_mat, nstatus3)
        tstatus_mat = vcat(tstatus_mat, tstatus3)
        VE_true_mat = vcat(VE_true_mat, VE_true3)
        samplesize_mat = vcat(samplesize_mat, samplesize3)
        TTE_mat = vcat(TTE_mat, TTE3)
        T = vcat(T, T3)
        R0 = vcat(R0, R03)
    end

    return soln_mat, nstatus_mat, tstatus_mat, VE_true_mat, samplesize_mat, TTE_mat, communitysize_arr, communitynum_arr, T, R0
end

function fn_iternation_iRCT_MLE(nsim, soln1, tstatus1, VE_true1, samplesize1, N, par_hh, par_community, par_prob, par_disease, prop_in_trial, import_lambda, casenum0, immunenum0, allocation_ratio, vac_efficacy, protection_threshold, treatment_gp, gamma_infectperiod_maxduration, trial_begintime, trial_endtime, endtime)

    soln1_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln1_mat[:,1] = soln1[:,1]
    soln1_mat[:,2] = soln1[:,2]
    soln1_mat[:,3] = soln1[:,3]
    soln1_mat[:,4] = soln1[:,4]
    soln1_mat[:,5] = soln1[:,5]

    (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
    (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)
    (nstatus2, tstatus2, soln2, T2, R02) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true2) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize2 = fn_samplesize_truecases(n_control, n_treatment, n_infected_control, n_infected_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE2 = fn_TTE(nstatus2, tstatus2, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    soln2_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
    soln2_mat[:,1] = soln2[:,1]
    soln2_mat[:,2] = soln2[:,2]
    soln2_mat[:,3] = soln2[:,3]
    soln2_mat[:,4] = soln2[:,4]
    soln2_mat[:,5] = soln2[:,5]
    soln_mat = vcat(soln1_mat, soln2_mat)
    nstatus_mat = vcat(nstatus1, nstatus2)
    tstatus_mat = vcat(tstatus1, tstatus2)
    VE_true_mat = vcat(VE_true1, VE_true2)
    samplesize_mat = vcat(samplesize1, samplesize2)
    TTE_mat = vcat(TTE1, TTE2)
    T = vcat(T1, T2)
    R0 = vcat(R01, R02)

    for itern = 3:nsim
        (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
        (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
        (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)
        (nstatus3, tstatus3, soln3, T3, R03) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

        timestep_fn = endtime
        (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true3) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
        samplesize3 = fn_samplesize_truecases(n_control, n_treatment, n_infected_control, n_infected_treatment, treatment_gp, timestep_fn, alpha, power)
        TTE3 = fn_TTE(nstatus3, tstatus3, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

        soln3_mat = zeros(Int, round(Int,endtime), 5)
        soln3_mat[:,1] = soln3[:,1]
        soln3_mat[:,2] = soln3[:,2]
        soln3_mat[:,3] = soln3[:,3]
        soln3_mat[:,4] = soln3[:,4]
        soln3_mat[:,5] = soln3[:,5]
        soln_mat = vcat(soln_mat, soln3_mat)
        nstatus_mat = vcat(nstatus_mat, nstatus3)
        tstatus_mat = vcat(tstatus_mat, tstatus3)
        VE_true_mat = vcat(VE_true_mat, VE_true3)
        samplesize_mat = vcat(samplesize_mat, samplesize3)
        TTE_mat = vcat(TTE_mat, TTE3)
        T = vcat(T, T3)
        R0 = vcat(R0, R03)
    end

    return soln_mat, nstatus_mat, tstatus_mat, VE_true_mat, samplesize_mat, TTE_mat, communitysize_arr, communitynum_arr, T, R0
end

function fn_divide(soln_mat, endtime, nsim, colnum)
    Y = zeros(Int, round(Int,endtime), nsim)
    for index in 1:nsim
        Y[:,index] = soln_mat[round(Int,endtime*(index-1)+1):round(Int,endtime*index),colnum]
    end
    return Y
end

# Function to count the occurrence of unique elements in v,
# in this case, an array with health and trial statuses
function fn_countelements2(v)

    # Input:
    # v: A vector that holds values which one wants to count the number of occurrences from

    # Output:
    # D: Number of occurrences of each unique value in v

    v_convert = fill("T", length(v))

    for i in 1:length(v)
        if v[i] == 1
            v_convert[i] = "S"
        elseif v[i] == 2
            v_convert[i] = "E"
        elseif v[i] == 3
            v_convert[i] = "I"
        elseif v[i] == 4
            v_convert[i] = "R"
        elseif v[i] == 5
            v_convert[i] = "ES"
        elseif v[i] == 6
            v_convert[i] = "EE"
        elseif v[i] == 7
            v_convert[i] = "EI"
        elseif v[i] == 8
            v_convert[i] = "ER"
        elseif v[i] == 9
            v_convert[i] = "V"
        elseif v[i] == 10
            v_convert[i] = "U"
        elseif v[i] == 11
            v_convert[i] = "VS"
        elseif v[i] == 12
            v_convert[i] = "VE"
        elseif v[i] == 13
            v_convert[i] = "VI"
        elseif v[i] == 14
            v_convert[i] = "VR"
        else
            throw(ArgumentError("Unrecognized status."))
        end
    end

    u = unique(v_convert)
    d = Dict([(i,count(x->x==i,v_convert)) for i in u])
    D = convert(DataFrame, d)

    D = insertcols!(D, 1, S=0, makeunique=true) # Set S as zero if column does not exist
    D = insertcols!(D, 1, E=0, makeunique=true) # Set E as zero if column does not exist
    D = insertcols!(D, 1, I=0, makeunique=true) # Set I as zero if column does not exist
    D = insertcols!(D, 1, R=0, makeunique=true) # Set R as zero if column does not exist
    D = insertcols!(D, 1, ES=0, makeunique=true) # ES = enrolled but not vaccinated, susceptible
    D = insertcols!(D, 1, EE=0, makeunique=true) # EE = enrolled but not vaccinated, exposed
    D = insertcols!(D, 1, EI=0, makeunique=true) # EI = enrolled but not vaccinated, infectious
    D = insertcols!(D, 1, ER=0, makeunique=true) # ER = enrolled but not vaccinated, removed
    D = insertcols!(D, 1, V=0, makeunique=true) # Set V as zero if column does not exist
    D = insertcols!(D, 1, VS=0, makeunique=true) # VS = enrolled and vaccinated, susceptible
    D = insertcols!(D, 1, VE=0, makeunique=true) # VE = enrolled and vaccinated, exposed
    D = insertcols!(D, 1, VI=0, makeunique=true) # VI = enrolled and vaccinated, infectious
    D = insertcols!(D, 1, VR=0, makeunique=true) # VR = enrolled and vaccinated, removed

    return D
end

function fn_translate_nstatus(N, nstatus_mat, tstatus_mat, endtime, nsim, colnum)

    # Initializations
    nstatus_mat_num = zeros(Int,length(nstatus_mat[:,colnum]))
    X = zeros(Int, N, nsim)

    for index1 in 1:length(nstatus_mat[:,colnum])
        if nstatus_mat[index1,colnum] == 'S' && tstatus_mat[index1,2] < 0
            nstatus_mat_num[index1] = 1
        elseif nstatus_mat[index1,colnum] == 'E' && tstatus_mat[index1,2] < 0
            nstatus_mat_num[index1] = 2
        elseif nstatus_mat[index1,colnum] == 'I' && tstatus_mat[index1,2] < 0
            nstatus_mat_num[index1] = 3
        elseif nstatus_mat[index1,colnum] == 'R' && tstatus_mat[index1,2] < 0
            nstatus_mat_num[index1] = 4
        elseif nstatus_mat[index1,colnum] == 'S' && tstatus_mat[index1,2] == 0
            nstatus_mat_num[index1] = 5
        elseif nstatus_mat[index1,colnum] == 'E' && tstatus_mat[index1,2] == 0
            nstatus_mat_num[index1] = 6
        elseif nstatus_mat[index1,colnum] == 'I' && tstatus_mat[index1,2] == 0
            nstatus_mat_num[index1] = 7
        elseif nstatus_mat[index1,colnum] == 'R' && tstatus_mat[index1,2] == 0
            nstatus_mat_num[index1] = 8
        elseif nstatus_mat[index1,colnum] == 'V'
            nstatus_mat_num[index1] = 9
        elseif nstatus_mat[index1,colnum] == 'U'
            nstatus_mat_num[index1] = 10
        elseif nstatus_mat[index1,colnum] == 'S' && tstatus_mat[index1,2] > 0
            nstatus_mat_num[index1] = 11
        elseif nstatus_mat[index1,colnum] == 'E' && tstatus_mat[index1,2] > 0
            nstatus_mat_num[index1] = 12
        elseif nstatus_mat[index1,colnum] == 'I' && tstatus_mat[index1,2] > 0
            nstatus_mat_num[index1] = 13
        elseif nstatus_mat[index1,colnum] == 'R' && tstatus_mat[index1,2] > 0
            nstatus_mat_num[index1] = 14
        else
            println(index1)
        end
    end
    for index2 in 1:nsim
        X[:,index2] = nstatus_mat_num[round(Int,N*(index2-1)+1):round(Int,N*index2)]
    end
    return X
end

# Function to return number of people in control and treatment_gp and those who are infectious in these groups
# Assume true observed cases are those infectious for now
function fn_vaccine_efficacy(tstatus, nstatus, timestep, treatment_gp)

    # Inputs:
    # tstatus: The trial statuses of everyone in the population
    # nstatus: The health statuses of everyone in the population at timestep
    # timestep: Time that end observation (from day 1 of outbreak)
    # treatment_gp: The treatment group number that we are comparing with control group

    # Outputs:
    # n_control: Number of people in the control group
    # n_treatment: Number of people in the treatment group
    # n_infectious_control: Number of infectious people in the control group
    # n_infectious_treatment: Number of infectious people in the treatment group
    # n_expoed_control: Number of exposed people in the control group
    # n_expoed_treatment: Number of exposed people in the treatment group
    # VE_true: Vaccine efficacy between control and treatment_gp

    # Determine the number of people who are in the control group
    nodes_control = findall(x->x==0, tstatus[:,2])
    n_control = length(nodes_control)

    # Determine the ones who have been infectious (I) from day 1 to t=timestep
    nodes_infectious1 = findall(x->x=='I', nstatus[:,1+1])
    nodes_infectious2 = findall(x->x=='I', nstatus[:,2+1])
    nodes_infectious = vcat(nodes_infectious1, nodes_infectious2)
    for t in 3:(round(Int,timestep))
        nodes_infectious3 = findall(x->x=='I', nstatus[:,t+1])
        nodes_infectious = vcat(nodes_infectious, nodes_infectious3)
    end
    nodes_infectious = unique(nodes_infectious) # Determine the cumulative incidence

    # Determine the cumulative incidence from day 1 to t=timestep of the ones in control group
    nodes_infectious_control = intersect(nodes_control, nodes_infectious)
    n_infectious_control = length(nodes_infectious_control)

    # Determine the ones who have been exposed (E) from day 1 to t=timestep
    nodes_exposed1 = findall(x->x=='E', nstatus[:,1+1])
    nodes_exposed2 = findall(x->x=='E', nstatus[:,2+1])
    nodes_exposed = vcat(nodes_exposed1, nodes_exposed2)
    for t in 3:(round(Int,timestep))
        nodes_exposed3 = findall(x->x=='E', nstatus[:,t+1])
        nodes_exposed = vcat(nodes_exposed, nodes_exposed3)
    end
    nodes_exposed = unique(nodes_exposed) # Determine the cumulative incidence

    # Determine the cumulative incidence from day 1 to t=timestep of the ones in control group
    nodes_exposed_control = intersect(nodes_control, nodes_exposed)
    n_exposed_control = length(nodes_exposed_control)

    # Determine the cumulative incidence of the ones in treatment_gp
    nodes_treatment = findall(x->x==treatment_gp, tstatus[:,2])
    n_treatment = length(nodes_treatment)
    nodes_infectious_treatment = intersect(nodes_treatment, nodes_infectious)
    n_infectious_treatment = length(nodes_infectious_treatment)
    nodes_exposed_treatment = intersect(nodes_treatment, nodes_exposed)
    n_exposed_treatment = length(nodes_exposed_treatment)

    # Compute vaccine efficacy based on true cases
    VE_true = 1-((n_infectious_treatment/n_treatment)/ (n_infectious_control/n_control))

    return n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true
end

# Function to return sample size of one treatment group based on number of true observed cumulative cases
function fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep, alpha, power)

    # Inputs:
    # n_control: Number of people in the control group
    # n_treatment: Number of people in the treatment group
    # n_infectious_control: Number of infectious people in the control group
    # n_infectious_treatment: Number of infectious people in the treatment group
    # treatment_gp: The treatment group number that we are comparing with control group
    # timestep: Time that end observation (from day 1 of outbreak)
    # alpha: Desired type I error
    # power: Power

    # Output:
    # samplesize_truecases: Sample size of one treatment group based on true observed cases, according to Sakpal et al (2010) Perspect Clin Res.

    # Compute critical value for alpha
    d = Normal()
    z_alpha = quantile(d, alpha)
    z_beta = quantile(d, power)

    p1 = (n_treatment - n_infectious_treatment)/ n_treatment # Proportion of people not infected in treatment group
    p2 = (n_control - n_infectious_control)/ n_control # Proportion of people not infected in control group

    # Compute sample size based on true cases
    samplesize_truecases = ((z_alpha+z_beta)^2 * ((p1*(1-p1) + p2*(1-p2))))/ (p1-p2)^2

    return samplesize_truecases
end

function fn_TTE(nstatus, tstatus, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    nodes_control = findall(x->x==0, tstatus[:,2])
    elements_control = findall(x->x=='I', nstatus[:,round(Int,trial_begintime+1)])

    nodes_treatment = findall(x->x==treatment_gp, tstatus[:,2])
    elements_treatment = findall(x->x=='I', nstatus[:,round(Int,trial_begintime+1)])

    elements = union(elements_control, elements_treatment)
    TTE_mat = zeros(Int, length(elements), 2) # Initialization

    bound_begin = trial_begintime + 1
    infectperiod_duration = trial_begintime + gamma_infectperiod_maxduration
    bound_end = minimum([infectperiod_duration, trial_endtime, endtime])+1

    for i in 1: length(elements)
        nstatus_timesegment = nstatus[elements[i], round(Int,bound_begin): round(Int,bound_end)]
        recovered_tmp = findall(x->x=='R', nstatus_timesegment)
        recovered_time = minimum(recovered_tmp)
        TTE_mat[i,1] = elements[i]
        TTE_mat[i,2] = recovered_time - 1
    end

    return TTE_mat
end

function fn_MLE_Rosenberger(p0,p1)
    allocation_prob_Rosenberger = sqrt(p1/p0)/((sqrt(p1/p0))+1)
end
function fn_MLE_Rosenberger_for_adapt(p)
    allocation_prob_Rosenberger = sqrt(p1/p0)/((sqrt(p1/p0))+1)
end

function fn_MLE_Neyman(p)
    allocation_prob_Neyman = (sqrt(p1[j]*(1-p1[j]))/(sqrt(p0[i]*(1-p0[i])) + sqrt(p1[j]*(1-p1[j]))))
    if p0[i]*(1-p0[i]) == 0 || p1[j]*(1-p1[j]) == 0
        allocation_prob = 0.5
        p0[i] = 1
        p1[j] = allocation_prob
    end
    return allocation_prob_Neyman
end

# Function to return estimated sample size required by the trial
function fn_testStat(tstatus, sigma, allocation_ratio, J, alpha, power, delta, e, f)

    # Inputs:
    # tstatus: The trial statuses of everyone in the population
    # sigma: The assumed value for the SD of the patient
    # allocation_ratio: The allocation ratio between the arms
    # J: Stage number of interim analysis of the trial
    # p: Effect size
    # p0: Treatment effect
    # alpha: Desired type I error
    # power: Power
    # delta: The clinical relevant different to power for
    # e: Efficacy stopping boundary
    # f: Futility stopping boundary

    # Outputs:
    # samplesize_est: Estimated sample size
    # Z: Wald test statistics
    # decision: decision of an arm to whether stop or continue

    # Initializations
    p = zeros(length(allocation_ratio)) # Proportion of population in different arms of the trial
    SD = zeros(length(allocation_ratio)) # Standard deviation
    SD_sq = zeros(length(allocation_ratio)) # Square of standard deviation
    n_arms = zeros(length(allocation_ratio)) # Number of population in each arm
    vac_prop = zeros(Int, length(allocation_ratio)) # Number of people who are in active arms
    r = zeros(length(allocation_ratio)) # Allocation ratio relative to the control group

    for i in 2: length(allocation_ratio)
        r[i] = allocation_ratio[i]/ allocation_ratio[1]
        if r[i] < 0
            r[i] = 1/r[i]
        end
    end

    p = allocation_ratio
    for i in 1:length(allocation_ratio)
        SD[i] = sqrt(p[i]*(1-p[i]))
        SD_sq[i] = (SD[i])^2
        #SE_sq[i] = (SD[i])^2/(N_est*p[i])
    end
    #SE = sqrt(sum(SE_sq))

    if alpha == 0.025
        z_alpha = 1.96
    else
        throw(ArgumentError("Need a z_alpha value."))
    end

    if power == 0.9
        z_beta = 1.28
    else
        throw(ArgumentError("Need a z_beta value."))
    end

    # Determine sample sizes
    n = ((z_alpha + z_beta)^2 * sum(SD_sq))/ delta^2
    n_arms[1] = length(allocation_ratio) * n * (1+1)
    for i in 2: length(allocation_ratio)
        n_arms[i] = length(allocation_ratio) * n * (1+1/r[i])
    end

    samplesize_est = sum(n_arms)
    samplesize_est_mean = mean(n_arms)
    samplesize_est_var = sqrt(var(n_arms))

    # Determine the number of people who are in active arms
    for i in 1:length(allocation_ratio)
        vac_prop = findall(x->x==(i-1), tstatus[:,2])
    end

    # Compute Wald test statistics by comparing treatment group(s) and control group
    for i in 1:length(allocation_ratio)
        Z[J,i] = (vac_prop[1] - vac_prop[i])/ sigma * sqrt((allocation_ratio[i] + allocation_ratio[1])/ (allocation_ratio[i]*allocation_ratio[1]*vac_prop[1]))
        if Z[J,i] <= e
            decision = "Stop for efficacy"
            println("Arm ", i-1," :", decision, " at stage ", J)
        elseif Z[J,i] > f
            decision = "Stop for futility"
            println("Arm ", i-1," :", decision, " at stage ", J)
        elseif e < Z[J,i] <= f
            decision = "Continue with adaption"
            println("Arm ", i-1," :", decision, " at stage ", J)
        end
    end

    return samplesize_est, Z, decision
end

function fn_adapt_freq_MLE(par0, method, n_control, n_treatment, n_infectious_control, n_infectious_treatment)

    # The allocation_prob is of treatment arm relative to control group

    # Inputs:
    # method: Method to use for frequentist adaptive
    # par0: initial start point for p=[p0,p1]
    # n_control: Number of people in the control group
    # n_treatment: Number of people in the treatment group
    # n_infected_control: Number of infectious people in the control group
    # n_infected_treatment: Number of infectious people in the treatment group

    # Output:
    # MLE: Maximum Likelihood Estimation of function fn_MLE(p)

    success_control = n_control - n_infectious_control
    success_treatment = n_treatment - n_infectious_treatment

    p0 = success_control/ n_control
    p1 = success_treatment/ n_treatment

    #if method == "Rosenberger"
        #P = sqrt(p1/p0)
        #allocation_prob = P/(P+1)
        optimum = optimize(fn_MLE_Rosenberger_for_adapt, par0)
        MLE = optimum.minimum
    #elseif method == "Neyman"
        #optimum = optimize(fn_MLE_Neyman, par0)
        #MLE = optimum.minium
    #else
        #throw(ArgumentError("Incorrect method defined."))
    #end

    return MLE
end

function fn_adapt_freq(method, timestep)

    # The allocation_prob is of treatment arm relative to control group

    # Inputs:
    # method: Method to use for frequentist adaptive
    # timestep: The time t of the outbreak

    # Output:
    # allocation_prob: The best allocation probability of treatment arm relative to control arm

    # Parameter values
    p0 = [0.1:0.1:1;]
    p1 = [0.1:0.1:1;]
    par0 = [p1[1],p0[1]]
    MLE = zeros(length(p0)*length(p1))
    element_selected_p0 = Float32[]
    element_selected_p1 = Float32[]
    sizehint!(element_selected_p0, length(p0)*length(p1))
    sizehint!(element_selected_p1, length(p0)*length(p1))

    if method == "Rosenberger"
        for i in 1:length(p0), j in 1:length(p1)
            P = sqrt(p1[j]/p0[i])
            allocation_prob = P/(P+1)
            fn_MLE_Rosenberger_for_adapt(p) = sqrt(p1/p0)/((sqrt(p1/p0))+1)
            optimum = optimize(fn_MLE_Rosenberger_for_adapt, [0.1, 0.1])
            MLE[counter1] = optimum.minimum
            if MLE[counter1] == minimum(MLE)
                push!(element_selected_p0::Array{Float64,1},p0[i])
                push!(element_selected_p1::Array{Float64,1},p1[j])
            end
            #global counter1 = counter1 + 1
        end
    elseif method == "Neyman"
        for i in 1:length(p0), j in 1:length(p1)
            optimum = optimize(fn_MLE_Neyman, par0)
            MLE[counter2] = optimum.minimum
            if MLE[counter2] == minimum(MLE)
                push!(element_selected_p0::Array{Float64,1},p0[i])
                push!(element_selected_p1::Array{Float64,1},p1[j])
            end
            #global counter2 = counter2 + 1
        end
    else
        throw(ArgumentError("Incorrect method defined."))
    end

    p0_selected = element_selected_p0[end]
    p1_selected = element_selected_p1[end]
    return p0_selected, p1_selected
end

function fn_adapt_Bayes(itern)

    fn(p0,p1) = -sqrt(p1/p0)/((sqrt(p1/p0))+1)

    # Parameters
    rewardnum0 = 0
    rewardnum1 = 0

    # Initializations
    soln = zeros(itern)
    posterior = zeros(2)
    element_selected3_p0 = Float64[]
    element_selected3_p1 = Float64[]
    sizehint!(element_selected3_p0, itern)
    sizehint!(element_selected3_p1, itern)

    for i in 1:itern
        p0 = rand(Uniform(), 1)[1]
        p1 = rand(Uniform(), 1)[1]
        while p0==0.0
            p0 = rand(Uniform(), 1)[1]
        end
        while p1==0.0
            p1 = rand(Uniform(), 1)[1]
        end
        soln[i] = fn(p0, p1)
        if soln[i] == minimum(soln)
            push!(element_selected3_p0::Array{Float64,1},p0)
            push!(element_selected3_p1::Array{Float64,1},p1)
            posterior = [p0, p1]
            rewardnum1 = rewardnum1 + 1
        else
            rewardnum0 = rewardnum0 + 1
        end
    end

    return posterior
end

function fn_casecounting(X, N, prop_in_trial, allocation_ratio, sim_num)
    XD = fn_countelements2(X[:,sim_num])
    unvac_uninfected = XD[:,:ES][1]
    unvac_infected = round(Int, N*prop_in_trial*allocation_ratio[1] - unvac_uninfected)
    vac_uninfected = XD[:,:VS][1] + XD[:,:V][1]
    vac_infected = round(Int, N*(1-prop_in_trial*allocation_ratio[1]) - vac_uninfected)
    infected = unvac_infected + vac_infected
    #p1 = vac_uninfected/ (N*prop_in_trial*allocation_ratio[treatment_gp])
    #p2 = unvac_uninfected/ (N*prop_in_trial*allocation_ratio[1])
    return infected
end

function fn_plot(X,Y)
    #plot0=([soln1[:,:S], soln1[:,:E], soln1[:,:I], soln1[:,:R], soln1[:,:V]], label=['S', 'E', 'I', 'R', 'V'])
    plot1=plot(Y, color=["grey"], ylabel="Number of infected individuals", legend=false)
    plot2=plot([lowerCI, medianCI, upperCI], label=["5% quantile", "Median", "95% quantile"], xlabel="Time (days)")
    return (plot1, plot2)
end
