# Stochastic block network model
# Author: Kendra Wu
# Date: 20 January 2020

# Simulates disease transmission of a 3-level clustered network with contact structure and imported cases at random time and clusters based on stochastic block model (SBM) using Ebola-like parameters from Hitchings et al (2018).

# Begin timing the processing time
@time begin #begin timing the processing time

# Introduce packages
using SpecialFunctions # For generating gamma distribution
using LinearAlgebra
using Distributions
using StatsBase # To use sample
#using Graphs
using CSV, DataFrames
using DelimitedFiles
using Plots

# Set parameter values
N = 500 # Total population size
endtime = 200.0 # Duration of simulation (timestep is in a unit of days)
S0 = N # Number of suspectible people in the population at day 0
V0 = 0 # Number of vaccinated individuals in the population at day 0
casenum0 = 10 # Number of infectious individuals introduced into the community to begin the outbreak
import_lambda = .05 # Number of occurrences variable for imported cases timeseries, assumed to follow Poisson Distribution

## The households
# hhnum: Number of households in the network
# hhsize_avg: Average size of one household
# hhsize_range: Range of household sizes (sizes are uniformly distributed)
par_hh = DataFrame(hhnum=150, hhsize_avg=3, hhsize_range=2)

## The clustered network
# communitynum: Number of communities in the clustered network
# communitysize_avg: Average size of one clustered community
# communitysize_range: Range of community sizes (sizes are uniformly distributed)
par_community = DataFrame(communitynum=3, communitysize_avg=100, communitysize_range=20)

## Contact and transmission probabilities between nodes
# cprob_hhwithin_cwithin: Probability of contacts of an edge between two nodes in the same household and the same community
# cprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# cprob_hhbetween_cbetween: Probability of transmission of an edge between two nodes in different communities
# tprob_hhwithin_cwithin: Probability of transmission of an edge between two nodes in the same household and the same community
# tprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# tprob_hhbetween_cbetween: Probability of contacts of an edge between two nodes in different communities
par_prob = DataFrame(cprob_hhwithin_cwithin=1, cprob_hhbetween_cwithin=0.6, cprob_hhbetween_cbetween=0.8, tprob_hhwithin_cwithin=0.6, tprob_hhbetween_cwithin=0.003, tprob_hhbetween_cbetween=0.0001)

## Disease properties
# Use Ebola-like parameters (from Hitchings (2018)) - Gamma-distributed
par_disease = DataFrame(incubperiod_shape=3.11, incubperiod_rate=0.32, infectperiod_shape=1.13, infectperiod_rate=0.226)

# Initializations
nstatus = fill("S", N, round(Int,endtime))
nstatus = hcat([1:1:N;], nstatus) # Put indexes on the 1st column
sbm_sol = DataFrame(S=fill(0,round(Int,endtime)), E=fill(0,round(Int,endtime)), I=fill(0,round(Int,endtime)), R=fill(0,round(Int,endtime)), V=fill(0,round(Int,endtime))) #Initialize the matrix which holds SEIR incidence of all timestep
V = zeros(Int8,V0) # V contains nodes_name of the vaccinated individuals, to be obtained from Cambridge

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
        clusterrandom_mat = rand(Uniform(-par_hh[1,:hhsize_range]/2,par_hh[1,:hhsize_range]/2), (par_hh[1,:hhnum], par_hh[1,:hhnum])) # Obtain a random household size value based on uniform distribution
        clusterrandom_mat = round.(Int8, clusterrandom_mat)
    elseif clustertype == "community"
        clusternum = par_community[1,:communitynum]
        clustersize_arr = zeros(Int, par_community[1,:communitynum]) # Holds the values of the sizes of the clusters
        clustersize_avg_mat = [(par_community[1,:communitysize_avg]) for i=1:par_community[1,:communitynum]]
        clusterrandom_mat = rand(Uniform(-par_community[1,:communitysize_range]/2,par_community[1,:communitysize_range]/2), (par_community[1,:communitynum], par_community[1,:communitynum])) # Obtain a random community size value based on uniform distribution
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
            clustersize_arr[i] = clustersize_avg_mat[i] + clusterrandom_mat[i,i]
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
    cluster_partition_arr = zeros(Int, size(clustersize_arr,1))

    for i = 1: (size(clustersize_arr,1))
        cluster_partition_arr[i] = clustersize_arr[i]
    end
    clusternum_arr = ones(Int, sum(clustersize_arr))

    # Define the node numbers where partitions occur
    if size(clustersize_arr,1) >= 2
        for q1 in 2:(size(clustersize_arr,1))
            cluster_partition_arr[q1] = cluster_partition_arr[q1-1] + clustersize_arr[q1]
        end
    end

    # Distribute the individuals into clusters according to their node number
    if size(cluster_partition_arr,1) >= 2
        for nodenum in 1:(sum(clustersize_arr))
            for q2 in 2:(size(clustersize_arr,1))
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

    # Initialization
    importcasenum_timeseries = zeros(Int8, (round(Int,endtime)))

    for t = 1:(round(Int,endtime))
        if t == 1
            importcasenum_timeseries[t] = casenum0
        else
            importcasenum_timeseries[t] = rand(Poisson(import_lambda), 1)[1]
        end
    end
    return importcasenum_timeseries
end

# Function to generate imported cases and randomly distribute them into different clusters (at timestep=1, for now)
function fn_importcases(par_disease, importcasenum_timeseries, nstatus, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # importcasenum_timeseries: Array that holds number of infectious individuals introduced into the community throughout the duration of the outbreak
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # nstatus: Health statuses of all individuals in the population at each time step

    # From importcasesnum_timeseries, obtain number of import cases at t=timestep
    casenum = importcasenum_timeseries[timestep]

    # Generate node_names of imported cases at t=timestep
    importcases = sample(1:N, casenum, replace=false) # Sampling without replacement

    for index1 in 1:(size(importcases,1))

        # Compute the parameters for disease properties
        incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
        infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
        incubperiod = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),1) #Incubation period of the disease
        infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),1) #Infectious period of the disease

        # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
        incubperiod[1] = 0
        tbound1 = min(timestep + ceil(incubperiod[1]), endtime)
        tbound2 = min(timestep + ceil(incubperiod[1]) + ceil(infectperiod[1]), endtime)

        # column_index start at 2 because nstatus[:,1] is nodes_name
        #for index2 in timestep:(round(Int,tbound1))
        #    nstatus[importcases[index1],index2+1] = "E"
        #end

        for index3 in (round(Int,tbound1)):(round(Int,tbound2))
            nstatus[importcases[index1],index3+1] = "I"
        end

        for index4 in (round(Int,tbound2)+1):(round(Int,endtime))
            nstatus[importcases[index1],index4+1] = "R"
        end
    end

    return nstatus
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

    return Gc
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
    for i = 1:size(Gc,1), j =1:size(Gc,2)
        if Gc[i,j] != 0 && nstatus[i,timestep+1] == "I" # Check if there is a contact edge between the pair and if the infector is infectious
            if communitynum_arr[i] == communitynum_arr[j] # Check if infector and infectee are from the same cluster
                if hhnum_arr[i] == hhnum_arr[j]
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
    network_index_arr2 = Int[]

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
    indexes_tmp = fill(0,size(nonzeros_indexes,1)) # To holds nodes_names of the individuals who are susceptible at (timestep+1) and they will be infected according to SBM model
    r = 1 # A counter

    for i in 1:(size(nonzeros_indexes,1))
        # Obtain nodes_name if nonzeros_indexes[i]'s status is susceptible at time=(timestep+1)
        if nstatus[nonzeros_indexes[i],timestep+1] == "S"
            indexes_tmp[r] = nstatus[nonzeros_indexes[i],:1]
            r += 1
        end
    end

    filter!(x->x≠0,indexes_tmp) # To remove the zero elements in indexes_tmp
    # Take into account susceptible deplection to find the minimum between the potential infectees and the available susceptibles
    bound = min(r-1, size(nonzeros_indexes,1), size(indexes_tmp,1)) # Find the minimum among the variables

    if bound> 0
        unique_indexes = sample(indexes_tmp, bound, replace=false)
        return unique_indexes
    else
        unique_indexes = Int8[]
        return unique_indexes
    end
end

# Function to simulate the spread of the disease and return the statuses of each nodes_name at all timesteps
function fn_spread(par_disease, nstatus, infectees, V, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # infectees: A list of nodes_names that are to be infected at t=timestep according to SBM
    # V: A list of nodes_names that are vaccinated at t=timestep
    # timestep: The time t of the outbreak

    # Output:
    # nstate: Health statuses of all the individuals in the population at each time step

    # From the infected cases in infectees, allow health statuses change as time progresses
    for index1 in 1:(size(infectees,1))

        # Compute the parameters for disease properties
        incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
        infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
        incubperiod = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),1) # Incubation period of the disease
        infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),1) # Infectious period of the disease

        # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
        tbound1 = min(timestep + ceil(incubperiod[1]), endtime)
        tbound2 = min(timestep + ceil(incubperiod[1]) + ceil(infectperiod[1]), endtime)

        # column_index start at 2 because nstatus[:,1] is nodes_name
        for index2 in timestep:(round(Int,tbound1))
            nstatus[infectees[index1],index2+1] = "E"
        end

        for index3 in (round(Int,tbound1)+1):(round(Int,tbound2))
            nstatus[infectees[index1],index3+1] = "I"
        end

        for index4 in (round(Int,tbound2)+1):(round(Int,endtime))
            nstatus[infectees[index1],index4+1] = "R"
        end
    end

    if size(V,1)>0
        for index5 in 1:(size(V,1))
            nstatus[V[index5],timestep:(round(Int,endtime)+1)] = "V"
        end
    end

    return nstatus
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
    return D
end

# Main algorithm
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
    nstatus = fn_importcases(par_disease, importcasenum_timeseries, nstatus, timestep3) # Import infectious cases at t-timestep3

    D = fn_countelements(nstatus[:,timestep3+1]) # Count number of occurrences of SEIRV at a particular t=timestep3
    D = insertcols!(D, 1, S=0, makeunique=true) # Set S as zero if column does not exist
    D = insertcols!(D, 1, E=0, makeunique=true) # Set E as zero if column does not exist
    D = insertcols!(D, 1, I=0, makeunique=true) # Set I as zero if column does not exist
    D = insertcols!(D, 1, R=0, makeunique=true) # Set R as zero if column does not exist
    D = insertcols!(D, 1, V=0, makeunique=true) # Set V as zero if column does not exist

    # Put these SEIRV incidence values into a DataFrame sbm_sol
    sbm_sol[timestep3,:S] = D[1,:S] # Put S value into the appropriate cell
    sbm_sol[timestep3,:E] = D[1,:E] # Put E value into the appropriate cell
    sbm_sol[timestep3,:I] = D[1,:I] # Put I value into the appropriate cell
    sbm_sol[timestep3,:R] = D[1,:R] # Put R value into the appropriate cell
    sbm_sol[timestep3,:V] = D[1,:V] # Put V value into the appropriate cell

# Begin transmission
for timestep1 in 2:(round(Int,endtime))

    if timestep1 != 2
        nstatus = fn_importcases(par_disease, importcasenum_timeseries, nstatus, timestep1)
    end

    Gt = fn_transmit_network(Gc, par_prob, hhnum_arr, communitynum_arr, nstatus, timestep1) # Construct a who-infect-whom stochastic block network based on the contact network Gc
    potential_transmit_indexes = fn_findnonzeros(Gt) # The index numbers that will have disease transmission according to the stochastic block network model
    transmit_indexes = fn_uniqueS(potential_transmit_indexes, nstatus, timestep1) # Check if potential_transmit_indexes are susceptibles

    if size(transmit_indexes,1)>0 # Check if there are infectees

        global nstatus
        global sbm_sol

        nstatus = fn_spread(par_disease, nstatus, transmit_indexes, V, timestep1) # Spread the diseae within the network and update nstatus

        # Count number of occurrences of SEIRV at a particular timestep
        for timestep2 in 2:(round(Int,endtime))

            D = fn_countelements(nstatus[:,timestep2]) # Count number of occurrences of SEIRV at a particular t=timestep2
            D = insertcols!(D, 1, S=0, makeunique=true) # Set S as zero if column does not exist
            D = insertcols!(D, 1, E=0, makeunique=true) # Set E as zero if column does not exist
            D = insertcols!(D, 1, I=0, makeunique=true) # Set I as zero if column does not exist
            D = insertcols!(D, 1, R=0, makeunique=true) # Set R as zero if column does not exist
            D = insertcols!(D, 1, V=0, makeunique=true) # Set V as zero if column does not exist

            # Put these SEIRV incidence values into a DataFrame sbm_sol
            sbm_sol[timestep2,:S] = D[1,:S] # Put S value into the appropriate cell
            sbm_sol[timestep2,:E] = D[1,:E] # Put E value into the appropriate cell
            sbm_sol[timestep2,:I] = D[1,:I] # Put I value into the appropriate cell
            sbm_sol[timestep2,:R] = D[1,:R] # Put R value into the appropriate cell
            sbm_sol[timestep2,:V] = D[1,:V] # Put V value into the appropriate cell
        end
    end
end

sbm_sol = sbm_sol[1:size(sbm_sol,1) .!= 1,: ] # Delete the dummy row 1 from sbm_sol
CSV.write("./data/sbm_sol.csv", sbm_sol, writeheader=true) # Write key results into files
#println(first(sbm_sol,10))

print("Processing time:")
end #stop timeing the processing time

# Plot graphs
plot([sbm_sol.S, sbm_sol.E, sbm_sol.I, sbm_sol.R, sbm_sol.V], label=["S","E","I","R","V"], title="Stochastic Block Model", xlabel="Time (days)", ylabel="Number of Population")
