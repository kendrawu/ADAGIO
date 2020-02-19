# Function to return the sizes of households or communities
function fn_par_cluster(N, par_hh, par_community, clustertype)

    # Inputs:
    # N: Total population size
    # par_hh: User-defined parameters of the households
    # par_community: User-defined parameters of the network clusters
    # clustertype: to decide whether we compute the sizes of households or communities in the network

    # Output:
    # clustersize_arr: Sizes of each cluster, whether they are households (hh) or communities

    # Assign parameter values
    if clustertype == "household"
        clusternum = par_hh[1,:hhnum]
        clustersize_arr = zeros(Int, par_hh[1,:hhnum]) # Holds the values of the sizes of the households
        clustersize_avg_mat = [(par_hh[1,:hhsize_avg]) for i=1:par_hh[1,:hhnum]]
        if par_hh[1,:hhsize_range] != 0
            clusterrandom_mat = rand(Uniform(-par_hh[1,:hhsize_range]/2,par_hh[1,:hhsize_range]/2), (par_hh[1,:hhnum], par_hh[1,:hhnum])) # Obtain a random household size value based on uniform distribution
        else
            clusterrandom_mat = zeros(Int, par_hh[1,:hhnum])
        end
        clusterrandom_mat = round.(Int8, clusterrandom_mat)
    elseif clustertype == "community"
        clusternum = par_community[1,:communitynum]
        clustersize_arr = zeros(Int, par_community[1,:communitynum]) # Holds the values of the sizes of the clusters
        clustersize_avg_mat = [(par_community[1,:communitysize_avg]) for i=1:par_community[1,:communitynum]]
        if par_community[1,:communitysize_range] != 0
            clusterrandom_mat = rand(Uniform(-par_community[1,:communitysize_range]/2,par_community[1,:communitysize_range]/2), (par_community[1,:communitynum], par_community[1,:communitynum])) # Obtain a random community size value based on uniform distribution
        else
            clusterrandom_mat = zeros(Int, par_community[1,:communitynum])
        end
        clusterrandom_mat = round.(Int8, clusterrandom_mat)
    else
        throw(ArgumentError("The decision variable needs to be either household or community."))
    end

    # Check if inputs are valid
    if clusternum<1
        throw(ArgumentError("The number of cluster needs to be at least 1."))
    end

    # Determine the size of each cluster
    for i in 1:clusternum
        if i != clusternum
            clustersize_arr[i] = clustersize_avg_mat[i] + clusterrandom_mat[i]
        else
            if (sum(clustersize_arr))>N
                throw(ArgumentError("Not enough people to partition them into all the clusters. Please rerun the simulation or consider lowering the number of clusters/ increasing N."))
            else
                clustersize_arr[i] = N - sum(clustersize_arr[1:(i-1)])
            end
        end
    end

    return clustersize_arr
end

# Function to define which node goes to which cluster and return a list of the cluster number
function fn_partition(clustersize_arr)

    # Input:
    # clustersize_arr: Sizes of each cluster

    # Output:
    # clusternum_arr: Cluster number of each individual in the population

    # Initialization
    cluster_partition_arr = zeros(Int, length(clustersize_arr))

    for i = 1: (length(clustersize_arr))
        cluster_partition_arr[i] = clustersize_arr[i]
    end
    clusternum_arr = ones(Int, sum(clustersize_arr))

    # Define the node numbers where partitions occur
    if size(clustersize_arr,1) >= 2
        for q1 in 2:(length(clustersize_arr))
            cluster_partition_arr[q1] = cluster_partition_arr[q1-1] + clustersize_arr[q1]
        end
    end

    # Distribute the individuals into clusters according to their node number
    if length(cluster_partition_arr) >= 2
        for nodenum in 1:(sum(clustersize_arr))
            for q2 in 2:(length(clustersize_arr))
                if (nodenum > cluster_partition_arr[q2-1]) && (nodenum <= cluster_partition_arr[q2])
                    clusternum_arr[nodenum] = q2
                end
            end
        end
    end

    return clusternum_arr
end

# Function to generate node_names of imported cases on a time series, assumed to be poisson distributed
function fn_importcases_timeseries(import_lambda, casenum0, endtime)

    # Inputs:
    # import_lambda: Number of occurrences variable for imported cases timeseries, assumed to follow Poisson Distribution
    # casenum0: Number of infectious individuals introduced into the community to begin the outbreak
    # endtime: Duration of simulation (timestep is in a unit of days)

    # Output:
    # importcasenum_timeseries: node_names of imported cases on a time series for the duration of the outbreak

    importcasenum_timeseries = rand(Poisson(import_lambda), round(Int,endtime))
    importcasenum_timeseries[1] = casenum0

    return importcasenum_timeseries
end

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
    importcases = sample(s_elements, casenum, replace=false) # Sampling without replacement

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

# Function to generate immune individuals and randomly distribute them into different clusters
# at the beginning of the outbreak
function fn_partialimmune(immunenum0, nstatus)

    # Inputs:
    # immunenum0: Number of people who are immune to the disease since day 0 of the outbreak
    # nstatus: Health statuses of all individuals in the population at each time step

    # Output:
    # nstatus_fn: Health statuses of all individuals in the population at each time step

    # Create a local variable for nstatus
    nstatus_fn = nstatus

    # Generate node_names of immuned people at t=timestep
    immuneppl = sample(1:N, immunenum0, replace=false) # Sampling without replacement

    for index1 in 1:(length(immuneppl))
        for index4 in 1:(round(Int,endtime))
            nstatus_fn[immuneppl[index1],index4+1] = 'R'
        end
    end

    return nstatus_fn
end

# Function to construct and return the who-contact-whom using stochastic block model
function fn_contact_network(par_prob, hhsize_arr, communitysize_arr, hhnum_arr, communitynum_arr)

    # Inputs:
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # hhsize_arr: Sizes of each household
    # communitysize_arr: Sizes of each community
    # hhnum_arr: Household number of each individuals in the population
    # communitynum_arr: Community number of each individuals in the population

    # Output:
    # Gc: The who-contact-whom stochastic block matrix graph

    # Initialization
    Gc = zeros(Int8, sum(communitysize_arr), sum(communitysize_arr))

    # Construct a who-contact-whom stochastic block matrix graph
    for i = 1:sum(communitysize_arr), j = 1:sum(communitysize_arr)
        if i != j # Check if infector and infectee are the same individual
            if communitynum_arr[i] == communitynum_arr[j] # Check if infector and infectee are from the same community
                if hhnum_arr[i] == hhnum_arr[j] # Check if infector and infectee are from the same houseshold
                    Gc[i,j] = rand(Bernoulli(par_prob[1,:cprob_hhwithin_cwithin]),1)[1]
                else
                    Gc[i,j] = rand(Bernoulli(par_prob[1,:cprob_hhbetween_cwithin]),1)[1]
                end
            else
                # Infector and infectee are from different communities, so they cannot be from the same household
                Gc[i,j] = rand(Bernoulli(par_prob[1,:cprob_hhbetween_cbetween]),1)[1]
            end
        end
    end

    return Gc
end

# Function to return the contact list of a infector
function fn_contactlist(Gc, nstatus, infector_node_name, timestep)

    # Inputs:
    # Gc: The who-contact-whom stochastic block matrix graph
    # nstatus: Health statuses of all individuals in the population at each time step
    # infector_node_name: Holds the node names of an infector

    # Output:
    # contact_arr: Contact list of an infector

    # Initialization
    contact_arr = zeros(Int8, N)

    # Find the contact list of infector_node_name
    if nstatus[infector_node_name, timestep+1] != 'I'
        println("wrong infector_node_name")
    else
        contact_arr = findall(x->x!=0, Gc[infector_node_name,:])
    end

    return contact_arr
end

# Function to return the contact and contac-of-contacts of an infectious person
function fn_computeT(N, par_prob, Gc, nstatus, hhnum_arr, communitynum_arr, timestep)

    # Inputs:
    # N: Total population size
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # Gc: The who-contact-whom stochastic block matrix graph
    # nstatus: Health statuses of all individuals in the population at each time step
    # hhnum_arr: Holds the household number of each individuals in the population
    # communitynum_arr: Holds the community number of each individuals in the population
    # timestep: The time t of the outbreak

    # Output:
    # T: Average probability that an infectious individual will transmit the disease to a susceptible individual with whom they have contact at t=timestep

    prob = zeros(N,N)
    prob_per_infector = zeros(N)

    for i = 1:size(Gc,1), j =1:size(Gc,2)
        if Gc[i,j] != 0 && nstatus[i,timestep+1] == 'I' # Check if there is a contact edge between the pair and if the infector is infectious
            if communitynum_arr[i] == communitynum_arr[j] # Check if infector and infectee are from the same community
                if hhnum_arr[i] == hhnum_arr[j] # Check if infector and infectee are from the same household
                    prob[i,j] = par_prob[1,:tprob_hhwithin_cwithin]
                else
                    prob[i,j] = par_prob[1,:tprob_hhbetween_cwithin]
                end
            else
                # Infector and infectee are from different communities, so they cannot be from the same household
                prob[i,j] = par_prob[1,:tprob_hhbetween_cbetween]
            end
        end
    end

    for q = 1:N
        prob_per_infector[q] = sum(prob[q,:])
    end

    filter!(x->~isnan(x),prob_per_infector)
    T = mean(prob_per_infector)
    return T
end

# Function to construct and return the who-infect-whom stochastic block matrix
function fn_transmit_network(Gc, par_prob, hhnum_arr, communitynum_arr, nstatus, timestep)

    # Inputs:
    # Gc: The who-contact-whom stochastic block matrix graph
    # par_prob: User-defined contact and transmission probabilities between two nodes
    # hhnum_arr: Holds the household number of each individuals in the population
    # communitynum_arr: Holds the community number of each individuals in the population
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # Gt: The who-infect-whom stochastic block matrix graph

    # Initialization
    Gt = zeros(Int8, size(Gc,1), size(Gc,2))

    # Construct a who-infect-whom stochastic block matrix graph
    for i = 1:size(Gc,1)
        if Gc[i,:] != 0 && nstatus[i,timestep+1] == 'I' # Check if the infector is infectious
            for j = 1:size(Gc,2) # Check if there is a contact edge between the infector-infectee pair
                if communitynum_arr[i] == communitynum_arr[j] # Check if infector and infectee are from the same community
                    if hhnum_arr[i] == hhnum_arr[j] # Check if infector and infectee are from the same household
                        Gt[i,j] = rand(Bernoulli(par_prob[1,:tprob_hhwithin_cwithin]),1)[1]
                    else
                        Gt[i,j] = rand(Bernoulli(par_prob[1,:tprob_hhbetween_cwithin]),1)[1]
                    end
                else
                    # Infector and infectee are from different communities, so they cannot be from the same household
                    Gt[i,j] = rand(Bernoulli(par_prob[1,:tprob_hhbetween_cbetween]),1)[1]
                end
            end
        end
    end

    return Gt
end

# Function to find the indexes that are non-zeros
function fn_findnonzeros(M)

    # Input:
    # M: A matrix

    # Output:
    # indexes: A list of indexes that are non-zeros in the matrix

    # Define the type of arrays
    network_index_arr1 = Int[]
    sizehint!(network_index_arr1, size(M,1))
    network_index_arr2 = Int[]
    sizehint!(network_index_arr2, size(M,2))

    # Find index numbers in M that are non-zeros, which indicates the probability of transmission between indexes i and j are non-zeros
    for i = 1:size(M,1), j = 1:size(M,2)
        if M[i,j] != 0
            push!(network_index_arr1::Array{Int,1},i)
            push!(network_index_arr2::Array{Int,1},j)
        end
    end

    # These lists return the unique indices, so that an infector or an infectee will not be counted twice
    x = unique(sort(network_index_arr1))
    nonzeros_indexes = unique(sort(network_index_arr2))

    return nonzeros_indexes
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

    if (length(nonzeros_indexes))>0
        for i in 1:(length(nonzeros_indexes))
            # Obtain nodes_name if nonzeros_indexes[i]'s status is susceptible at time=(timestep+1)
            if nstatus[nonzeros_indexes[i],timestep+1] == 'S'
                push!(s_elements::Array{Int,1},nstatus[nonzeros_indexes[i],1])
            end
        end

        # Take into account susceptible deplection to find the minimum between the potential infectees and the available susceptibles
        bound = min(length(nonzeros_indexes), length(s_elements)) # Find the minimum among the variables

        if bound> 0
            unique_indexes = sample(s_elements, bound, replace=false)
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

# Function to simulate the spread of the disease and return the statuses of each nodes_name at all timesteps
function fn_spread(par_disease, nstatus, infectees, tstatus, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # nstate: Health statuses of all the individuals in the population at each time step
    # infectees: A list of nodes_names that are to be infected at t=timestep according to SBM
    # tstatus: Trial statuses of everyone in the population at timestep
    # timestep: The time t of the outbreak

    # Output:
    # nstatus_fn: Health statuses of all the individuals in the population at each time step

    # Create a local variable for nstatus
    nstatus_fn = nstatus

    # Compute the parameters for disease properties
    #incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
    #infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
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

    V = findall(x->x>0, tstatus[:,2])
    if length(V)>0
        println("V = ", V, length(V))
        for index5 in 1:(length(V))
            for t in timestep:(round(Int,endtime))
                #println(index5, t)
                nstatus_fn[V[index5],t+1] = 'V'
            end
        end
    end

    return nstatus_fn
end

# Function to count the occurrence of unique elements in nstatus[:,timestep],
# in this case, the function counts the incidence of SEIRV nodes
function fn_countelements(v)

    # Input:
    # v: A vector that holds values which one wants to count the number of occurrences from

    # Output:
    # D: Number of occurrences of each unique value in v

    u = unique(v)
    d = Dict([(i,count(x->x==i,v)) for i in u])
    D = convert(DataFrame, d)

    D = insertcols!(D, 1, S=0, makeunique=true) # Set S as zero if column does not exist
    D = insertcols!(D, 1, E=0, makeunique=true) # Set E as zero if column does not exist
    D = insertcols!(D, 1, I=0, makeunique=true) # Set I as zero if column does not exist
    D = insertcols!(D, 1, R=0, makeunique=true) # Set R as zero if column does not exist
    D = insertcols!(D, 1, V=0, makeunique=true) # Set V as zero if column does not exist

    return D
end

function fn_transmodel_nsim1(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)

    # Initializations
    nstatus = fill('S', N, round(Int,endtime)) # health statuses
    nstatus = hcat([1:1:N;], nstatus) # Put indexes on the 1st column
    tstatus = fill(-1, N) # trial statuses: (-1=not in trial, 0=in control, 1=vaccine candidate number)
    tstatus = hcat([1:1:N;], tstatus) # Put indexes on the 1st column
    sbm_sol = DataFrame(S=fill(0,round(Int,endtime)), E=fill(0,round(Int,endtime)), I=fill(0,round(Int,endtime)), R=fill(0,round(Int,endtime)), V=fill(0,round(Int,endtime))) #Initialize the matrix which holds SEIR incidence of all timestep
    T_arr = zeros(round(Int,endtime))
    V = zeros(Int8, V0)

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

    # Begin transmission
    for timestep1 in 2:(round(Int,endtime))

        # Set local variables
        nstatus_fn = nstatus
        tstatus_fn = tstatus
        sbm_sol_fn = sbm_sol

        nstatus_fn = fn_importcases(par_disease, importcasenum_timeseries, nstatus_fn, timestep1) # Import cases
        Gt = fn_transmit_network(Gc, par_prob, hhnum_arr, communitynum_arr, nstatus_fn, timestep1) # Construct a who-infect-whom stochastic block network based on the contact network Gc
        potential_transmit_indexes = fn_findnonzeros(Gt) # The index numbers that will have disease transmission according to the stochastic block network model
        transmit_indexes = fn_uniqueS(potential_transmit_indexes, nstatus_fn, timestep1) # Check if potential_transmit_indexes are susceptibles
        T_arr[timestep1] = fn_computeT(N, par_prob, Gc, nstatus_fn, hhnum_arr, communitynum_arr, timestep1)

        if size(transmit_indexes,1)>0 # Check if there are infectees

            nstatus_fn = fn_spread(par_disease, nstatus_fn, transmit_indexes, tstatus_fn, timestep1) # Spread the diseae within the network and update nstatus

            # Count number of occurrences of SEIRV at a particular timestep
            for timestep2 in timestep1:(round(Int,endtime))

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

        # Compute R0 in a network
        if timestep1 == round(Int,endtime)
            k = sum(sum(Gc))/N # mean degree of the network
            T = mean(T_arr[T_arr .> 0]) # Average probability that an infectious individual will transmit the disease to a susceptible individual with whom they have contact
            R0 = T * (k^2/k - 1)
            println("k = ", k, ", T = ",T, ", and R0 = ",R0)

            global T
            global R0
        end
    end

    sbm_sol_mat = zeros(Int, round(Int,endtime), 5) # Same as sbm_sol, except it is a matrix, not a DataFrame
    sbm_sol_mat[:,1] = sbm_sol[:,1]
    sbm_sol_mat[:,2] = sbm_sol[:,2]
    sbm_sol_mat[:,3] = sbm_sol[:,3]
    sbm_sol_mat[:,4] = sbm_sol[:,4]
    sbm_sol_mat[:,5] = sbm_sol[:,5]
    return sbm_sol_mat, nstatus, tstatus, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, T, R0
end
