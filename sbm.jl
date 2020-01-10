# Stochastic block network model
# Authors: Kendra Wu
# Date: 27 November 2019

# Simulates disease transmission of a 2-level clustered network based on stochastic block model (SBM) using Ebola-like parameters from Hitchings et al (2018).

# Begin timing the processing time
@time begin #begin timing the processing time

# Introduce packages
using SpecialFunctions # For generating gamma distribution
using LinearAlgebra
using Distributions
using CSV, DataFrames
using DelimitedFiles
using Plots

# Set parameter values
N = 500 # Total population size
endtime = 100.0 # Duration of simulation (timestep is in a unit of days)
S0 = N # Number of suspectible people in the population at day 0
V0 = 0 # Number of vaccinated individuals in the population at day 0
casenum0 = 10 # Number of infectious individuals introduced into the community to begin the outbreak
init_seir = DataFrame(s=S0-casenum0, e=0, i=casenum0, r=0, v=V0) # Initial conditions of SEIRV

## The clustered network
#communitynum: Number of communities in the clustered network
#communitysize_avg: Average size of one clustered community
#communitysize_range: Range of community sizes (sizes are uniformly distributed)
#prob_within: Probability of an edge between two nodes in the same community
#prob_between: Probability of an edge between two nodes in different communities
par_cluster = DataFrame(communitynum=2, communitysize_avg=100, communitysize_range=40, prob_within=0.003, prob_between=0.0001)

## Disease properties
lambda = 0.01 # Per-time-step harzard of infection for a sucsptible nodes from an infectious neighbour
# Use Ebola-like parameters (from Hitchings (2018)) - Gamma-distributed
par_disease = DataFrame(incubperiod_shape=3.11, incubperiod_rate=0.32, infectperiod_shape=1.13, infectperiod_rate=0.226)

# Initializations
beta_t = ones(round(Int,endtime+1)) #Initialize beta_t
nstatus = convert(DataFrame,Array{Union{Missing, String}}(missing, N, round(Int,endtime))) # Holds the health statuses of all the individuals in the population at each time step
for r1=1:N, r2=1:round(Int,endtime)
    nstatus[r1,r2] = "S" # Set the initial health status as "S" for everyone in the population at all time steps
end
nstatus = insertcols!(nstatus, 1, nodes_name=1:N) # Insert a column with column name nodes_name and values as 1:N
sbm_sol = DataFrame(S=fill(0,round(Int,endtime)), E=fill(0,round(Int,endtime)), I=fill(0,round(Int,endtime)), R=fill(0,round(Int,endtime)), V=fill(0,round(Int,endtime)), N=fill(0,round(Int,endtime))) #Initialize the matrix which holds SEIR incidence of all timestep
V = zeros(Int,V0) # V contains nodes_name of the vaccinated individuals, to be obtained from Cambridge

# Function to return rate of transmission, beta
function func_beta(beta,p,t)
    for t = 1: (trunc(Int,tspan[2])+1)
        # p[1]=betahat; p[2]=a1; p[3]=a2; p[4]=atau; p[5]=sigma; p[6]=gamma
        beta[t] = p[1] * (1 - p[3]/(1 + exp(-p[2] * (t - p[4]))))
    end
end

# Function to return the sizes of clusters
function func_par_cluster(par_cluster)

    # Input:
    # par_cluster: User-defined parameters of the network clusters

    # Output:
    # communitysize: Sizes of each cluster

    communitysize = zeros(Int, par_cluster[1,:communitynum]) # Holds the values of the sizes of the clusters

    communitysize_avg_mat = [(par_cluster[1,:communitysize_avg]) for i=1:par_cluster[1,:communitynum]]
    random_mat = rand(Uniform(-par_cluster[1,:communitysize_range]/2,par_cluster[1,:communitysize_range]/2), (par_cluster[1,:communitynum], par_cluster[1,:communitynum])) #Obtain a random community size value based on uniform distribution
    random_mat = round.(Int, random_mat)

    # Determine the sizes of each cluster
    for i in 1:(par_cluster[1,:communitynum])
        if i != par_cluster[1,:communitynum]
            communitysize[i] = communitysize_avg_mat[i] + random_mat[i,i]
        else
            if (sum(communitysize))>N
                throw(ArgumentError("Not enough people to partition them into all the clusters. Please rerun the simulation or consider lowering the number of clusters/ increasing N."))
            else
                communitysize[i] = N - sum(communitysize[1:(i-1)])
            end
        end
    end

    return communitysize
end

# Function to define which node goes to which cluster and return a list of the cluster number
function func_partition(communitysize)

    # Input:
    # communitysize: Sizes of each cluster

    # Output:
    # clusternum: Holds the cluster number of each individuals in the population

    # Initialization
    community_partition = zeros(Int, size(communitysize,1), size(communitysize,1))

    for i = 1: (size(communitysize,1))
        community_partition[i] = communitysize[i]
    end
    clusternum = ones(Int, sum(communitysize))

    # Define the node numbers where partitions occur
    if size(communitysize,1) >= 2
        for q1 in 2:(size(communitysize,1))
            community_partition[q1] = community_partition[q1-1] + communitysize[q1]
        end
    end

    # Distribute the individuals into clusters according to their node number
    if size(community_partition,1) >= 2
        for nodenum in 1:(sum(communitysize))
            for q2 in 2:(size(communitysize,1))
                if (nodenum >= community_partition[q2-1]) && (nodenum <= community_partition[q2])
                    clusternum[nodenum] = q2
                end
            end
        end
    end

    return clusternum
end

# Function to generate imported cases and randomly distribute them into different clusters (at timestep=1, for now)
function func_importcases(casenum0, par_disease, nstatus, timestep)

    # Inputs:
    # casenum0: Number of infectious individuals introduced into the community to begin the outbreak
    # par_disease: User-defined parameters of the disease
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # nstatus: Health statuses of all individuals in the population at each time step

    # Generate node_names of imported cases
    importcases = rand(1:N, casenum0)

    # From the infected cases in G, allow health statuses change as time progresses
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
        #    nstatus[G[index1],index2+1] = "E"
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

# Function to construct and return the who-infect-whom stochastic block matrix
function func_network(par_cluster, communitysize, clusternum, nstatus, timestep)

    # Inputs:
    # par_cluster: User-defined parameters of the network clusters
    # communitysize: Sizes of each cluster
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # P: The who-infect-whom stochastic block matrix graph

    # Initializations
    P = zeros(Int64, sum(communitysize), sum(communitysize))

    # Construct a who-infect-whom stochastic block matrix graph
    for i = 1:sum(communitysize)
        if nstatus[i,timestep+1] != "I" # (timestep+1) because nstatus[,1] holds the nodes_names
            for j = 1:sum(communitysize)
                P[i,j] = 0
            end
        else
            for j = 1:sum(communitysize)
                if clusternum[i] == clusternum[j] # Check if infector and infectee are from the same cluster
                    P[i,j] = rand(Bernoulli(par_cluster[1,:prob_within]),1)[1] # To obtain the value from the array
                else
                    P[i,j] = rand(Bernoulli(par_cluster[1,:prob_between]),1)[1] # To obtain the value from the array
                end
            end
        end
    end

    return P
end

# Function to find the indexes in P that are non-zeros which indicates a transmission and return a list of indexes of potential infectees
function func_find(P)

    # Input:
    # P: The who-infect-whom stochastic block matrix

    # Output:
    # y: A list of potential infectees from the who-infect-whom stochastic block model

    # Define the type of arrays
    network_index_arr1 = Int[]
    network_index_arr2 = Int[]

    # Find index numbers in P that are non-zeros, which indicates the probability of transmission between indexes i and j are non-zeros
    for i = 1:size(P,1), j = 1:size(P,2)
        if P[i,j] != 0
            push!(network_index_arr1::Array{Int64,1},i)
            push!(network_index_arr2::Array{Int64,1},j)
        end
    end

    # These lists return the unique indices, so that an infector or an infectee will not be counted twice
    x = unique(sort(network_index_arr1))
    y = unique(sort(network_index_arr2))

    return y
end

# From y, this function finds and return the node_names of those who are susceptible at t=(timestep+1)
function func_uniqueS(y, nstatus, timestep)

    # Inputs:
    # y: A list of potential infectees from the who-infect-whom SBM
    # nstatus: Health statuses of all individuals in the population at each time step
    # timestep: The time t of the outbreak

    # Output:
    # G: A list of nodes_names that are to be infected at t=(timestep+1) according to SBM

    # Initialization
    G1 = fill(0,size(y,1)) # To holds nodes_names of the individuals who are susceptible at (timestep+1) and they will be infected according to SBM model
    r = 1 # A counter

    for i in 1:(size(y,1))
        # Obtain nodes_name if y[i]'s status is susceptible at time=(timestep+1)
        if nstatus[y[i],timestep+1] == "S"
            G1[r] = nstatus[y[i],:nodes_name]
            r += 1
        end
    end

    filter!(x->x≠0,G1) # To remove the zero elements in G1
    # Take into account susceptible deplection
    bound = min(r-1, size(y,1), size(G1,1)) # Find the minimum among the variables

    if bound> 0
        G = rand(G1, bound)
        return G
    else
        G = Int[]
        return G
    end
end

# Function to simulate the spread of the disease and return the statuses of each nodes_name at all timesteps
function func_spread(par_disease, nstatus, G, V, timestep)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # G: A list of nodes_names that are to be infected at t=timestep according to SBM
    # V: A list of nodes_names that are vaccinated at t=timestep
    # timestep: The time t of the outbreak

    # Output:
    # nstate: Health statuses of all the individuals in the population at each time step

    # From the infected cases in G, allow health statuses change as time progresses
    for index1 in 1:(size(G,1))

        # Compute the parameters for disease properties
        incubperiod_avg = ceil(par_disease[1,:incubperiod_shape]/par_disease[1,:incubperiod_rate])
        infectperiod_avg = ceil(par_disease[1,:infectperiod_shape]/par_disease[1,:infectperiod_rate])
        incubperiod = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),1) #Incubation period of the disease
        infectperiod = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),1) #Infectious period of the disease

        # Set time boundaries according to incubation and infectious periods, and time should not exceed endtime.
        tbound1 = min(timestep + ceil(incubperiod[1]), endtime)
        tbound2 = min(timestep + ceil(incubperiod[1]) + ceil(infectperiod[1]), endtime)

        # column_index start at 2 because nstatus[:,1] is nodes_name
        for index2 in timestep:(round(Int,tbound1))
            nstatus[G[index1],index2+1] = "E"
        end

        for index3 in (round(Int,tbound1)+1):(round(Int,tbound2))
            nstatus[G[index1],index3+1] = "I"
        end

        for index4 in (round(Int,tbound2)+1):(round(Int,endtime))
            nstatus[G[index1],index4+1] = "R"
        end
    end

    if size(V,1)>0
        for index5 in 1:(size(V,1))
            nstatus[V[index5],timestep:(round(Int,endtime)+1)] = "V"
        end
    end

    return nstatus
end

# Function to count the occurrence of unique elements in vector sbm.nodes_status,
# in this case, the function counts the incidence of SEIRV nodes
function func_countelements(v)

    # Input:
    # v: A vector that holds values which one wants to count the number of occurrences from

    # Output:
    # D: Number of occurrences of each unique value in v

    u = unique(v)
    d = Dict([(i,count(x->x==i,v)) for i in u])
    D = convert(DataFrame, d)
    return D
end

# Compute the parameters of the clusters
communitysize = func_par_cluster(par_cluster) # Define the sizes of each cluster
clusternum = func_partition(communitysize) # Assign cluster number to each individual in the population

# Import infectious cases (at timestep=1, for now)

    nstatus = func_importcases(casenum0, par_disease, nstatus, 1)

    timestep3=1
    D = func_countelements(nstatus[:,timestep3+1]) # Count number of occurrences of SEIRV at a particular t=timestep2
    D = insertcols!(D, 1, S=0, makeunique=true) #Set S as zero if column does not exist
    D = insertcols!(D, 1, E=0, makeunique=true) #Set E as zero if column does not exist
    D = insertcols!(D, 1, I=0, makeunique=true) #Set I as zero if column does not exist
    D = insertcols!(D, 1, R=0, makeunique=true) #Set R as zero if column does not exist
    D = insertcols!(D, 1, V=0, makeunique=true) #Set V as zero if column does not exist

    # Put these SEIRV incidence values into a DataFrame sbm_sol
    sbm_sol[timestep3,:S] = D[1,:S] #Put S value into the appropriate cell
    sbm_sol[timestep3,:E] = D[1,:E] #Put E value into the appropriate cell
    sbm_sol[timestep3,:I] = D[1,:I] #Put I value into the appropriate cell
    sbm_sol[timestep3,:R] = D[1,:R] #Put R value into the appropriate cell
    sbm_sol[timestep3,:V] = D[1,:V] #Put V value into the appropriate cell
    sbm_sol[timestep3,:N] = sbm_sol[timestep3,:S] + sbm_sol[timestep3,:E] + sbm_sol[timestep3,:I] + sbm_sol[timestep3,:R] + sbm_sol[timestep3,:V] # The total number of the population, for accurancy check


# Main algorithm
for timestep1 in 2:(round(Int,endtime))

    P = func_network(par_cluster, communitysize, clusternum, nstatus, timestep1) # Construct a who-infect-whom stochastic block network
    y = func_find(P) # The index numbers that will have disease transmission according to the stochastic block network model
    G = func_uniqueS(y, nstatus, timestep1) # Make sure the infectees by nodes_names from y are susceptibles

    if size(y,1)>0 && size(G,1)>0

        global nstatus
        global sbm_sol

        nstatus = func_spread(par_disease, nstatus, G, V, timestep1) # Spread the diseae within the network and update nstatus

        # Count number of occurrences of SEIRV at a particular timestep
        for timestep2 in 2:(round(Int,endtime))

            #global sbm_sol

            D = func_countelements(nstatus[:,timestep2]) # Count number of occurrences of SEIRV at a particular t=timestep2
            D = insertcols!(D, 1, S=0, makeunique=true) #Set S as zero if column does not exist
            D = insertcols!(D, 1, E=0, makeunique=true) #Set E as zero if column does not exist
            D = insertcols!(D, 1, I=0, makeunique=true) #Set I as zero if column does not exist
            D = insertcols!(D, 1, R=0, makeunique=true) #Set R as zero if column does not exist
            D = insertcols!(D, 1, V=0, makeunique=true) #Set V as zero if column does not exist

            # Put these SEIRV incidence values into a DataFrame sbm_sol
            sbm_sol[timestep2,:S] = D[1,:S] #Put S value into the appropriate cell
            sbm_sol[timestep2,:E] = D[1,:E] #Put E value into the appropriate cell
            sbm_sol[timestep2,:I] = D[1,:I] #Put I value into the appropriate cell
            sbm_sol[timestep2,:R] = D[1,:R] #Put R value into the appropriate cell
            sbm_sol[timestep2,:V] = D[1,:V] #Put V value into the appropriate cell
            sbm_sol[timestep2,:N] = sbm_sol[timestep2,:S] + sbm_sol[timestep2,:E] + sbm_sol[timestep2,:I] + sbm_sol[timestep2,:R] + sbm_sol[timestep2,:V]#The total number of the population, for accurancy check
        end
    end
end

sbm_sol = sbm_sol[1:size(sbm_sol,1) .!= 1,: ] # Delete the dummy row 1 from sbm_sol
CSV.write("./data/sbm_sol.csv", sbm_sol, writeheader=true) # Write key results into files
#println(last(sbm_sol,10))

# Plot graphs
plot([sbm_sol.S, sbm_sol.E, sbm_sol.I, sbm_sol.R, sbm_sol.V], label=["S","E","I","R","V"], title="Stochastic Block Model", xlabel="Time (days)", ylabel="Number of Population")

print("Processing time:")
end #stop timeing the processing time
