# Block network model like Hitchings using DifferentialEquations
# Authors: Kendra Wu
# Date: 12 November 2019

# This is a stochastic block model (SBM) of a clustered network using Ebola-like
# parameters from Hitchings et al (2018). Vaccinated individuals are not included.

#Begin timing the processing time
@time begin #begin timing the processing time

#Introduce packages
using SpecialFunctions #For generating gamma distribution
using LinearAlgebra
using Distributions
using LightGraphs #For generating Stochastic Block Model (SBM)
using CSV, DataFrames
using DelimitedFiles
using Plots
#using IJulia # For deprecation warnings suppression

# To suppression deprecation warnings
# IJulia.installkernel("Julia nodeps", "--depwarn=no")

# Set parameter values
N = 500 #Total population size
endtime = 100.0 #Duration of simulation (in days)
par_ode = [0.94,0.19,0.6,27.79,0.14,0.33] #Parameter values of ODE system for beta
tspan = (0.0,endtime)
case0num = 10 #Number of infectious individuals introduced into the community to begin the outbreak
## The clustered network
communitynum = 1 #Number of communities in the clustered network
communitysize_avg = 100 #Average size of one clustered community
communitysize_range = 40 #Range of community sizes (sizes are uniformly distributed)
rates_within = 0.15 #Probability of an edge between two nodes in the same community
rates_between = 0.001 #Probability of an edge between two nodes in different communities
## Disease properties - Gamma-distributed
incubperiod_shape = 3.11
incubperiod_rate = 0.32
infectperiod_shape = 1.13
infectperiod_rate = 0.226
incubperiod_avg = ceil(incubperiod_shape/incubperiod_rate)
infectperiod_avg = ceil(infectperiod_shape/infectperiod_rate)

# Initializations
beta_t = ones(round(Int,endtime+1)) #Initialize beta_t
InfectSize = case0num #Initialize the maximum number of infected individuals at each timestep
sbm = hcat(DataFrame(nodes_name=1:N, nodes_status=Vector{Union{Missing, String}}(missing, N), e_nodes_infectday=zeros(N), e_nodes_incubperiod=zeros(N), i_nodes_infectperiod=zeros(N)), makeunique=true)
sbm.nodes_status = "S" #Set initial status of each node as susceptible
nstatus = convert(DataFrame,Array{Union{Missing, String}}(missing, N, round(Int,endtime+1)))
for r1=1:N, r2=1:round(Int,endtime)
    nstatus[r1,r2] = "S"
end
nstatus = insertcols!(nstatus, 1, nodes_name=1:N)
sbm_sol = DataFrame(S=fill(0,round(Int,endtime+1)), E=fill(0,round(Int,endtime+1)), I=fill(0,round(Int,endtime+1)), R=fill(0,round(Int,endtime+1)), N=fill(0,round(Int,endtime+1))) #Initialize the matrix which holds SEIR incidence of all timestep

# Function to return rate of transmission, beta
function func_beta(beta,p,t)
    for t = 1: (trunc(Int,tspan[2])+1)
        # p[1]=betahat; p[2]=a1; p[3]=a2; p[4]=atau; p[5]=sigma; p[6]=gamma
        beta[t] = p[1] * (1 - p[3]/(1 + exp(-p[2] * (t - p[4]))))
    end
end

# Function to construct the who-infect-whom network and return it as a DataFrame
function func_network(N)
    # g = stochastic_block_model(studypopsize, rates_mat)
    # Construct a who-infect-whom network
    g = SimpleGraph{Int16}(N, N)

    # Format the who-infect-whom network into a DataFrame
    savegraph("./data/sbm.csv", g) #Write results into a file
    #For PC
    x = CSV.read(joinpath(dirname(pathof(DataFrames)), "C:/Users/kendraw/OneDrive - Nexus365/ADAGIO/Scripts/data/sbm.csv")) #Convert g into DataFrame and name it as x
    #For Mac                                           
    #x = CSV.read(joinpath(dirname(pathof(DataFrames)), "/Users/kendramwu/OneDrive - Nexus365/ADAGIO/Scripts/data/sbm.csv")) #Convert g into DataFrame and name it as x
    x = select!(x, Not(3:7)) #Delete the junk columns
    x = DataFrame(Infector=x[:,1], Infectee=x[:,2]) #Add column names to x

    #Find the unique values in x, so that no node is infected twice
    y = unique(x.Infectee)

    #G = dummy[1:InfectSize] #Take the first InfectSize number from dummy
    return y
end

#From the list y, this function finds the node_names who are susceptible at
#time = timestep
function func_uniqueS(y, nstatus, InfectSize, timestep)

    # Initialization
    G = fill(0,size(y,1))
    r = 1

    for i in 1:(size(y,1))
        # Obtain nodes_name if y[i]'s status is susceptible at time=timestep
        if nstatus[y[i],timestep] == "S"
            G[r] = nstatus[y[i],:nodes_name]
            r = r + 1
        end
    end

    #Set an upper bound of how many susceptible individuals are infected at each
    #timestep. Take into account susceptible deplection.
    bound = min(r-1, size(y,1), InfectSize)
    G = G[1:bound]
    return G
end

# Function to simulate the spread of disease within the network and return the nodes' names
# Note the listed status of a particular node is at time = timestep
function func_S2E(G, timestep)
    for index1 in 1:(size(G,1))
        sbm.nodes_status[G[index1]] = "E"
        sbm.e_nodes_infectday[G[index1]] = sbm.e_nodes_infectday[G[index1]] + 1
        sbm.e_nodes_incubperiod[G[index1]] = timestep + incubperiod_avg
        sbm.i_nodes_infectperiod[G[index1]] = timestep + incubperiod_avg + infectperiod_avg
    end
    return sbm
end

function func_E2I(sbm, G, timestep)
    # From the exposed cases in sbm, allow health statuses change as time progresses
    for index1 in 1:(size(G,1))

        #Set time boundaries. Time should not exceed endtime.
        tbound1 = min(timestep + incubperiod_avg, endtime)
        tbound2 = min(timestep + incubperiod_avg + infectperiod_avg, endtime)

        # colnum doesn't start at 2 because i) nstatus[:,1] is nodes_name
        # ii) timestep=1 is when infected cases were introduced.
        # Thus, the transmission begins at colnum=3
        for index2 in (timestep+2):(round(Int,tbound1)+2)
            nstatus[G[index1],index2] = "E"
        end

        for index3 in (round(Int,tbound1)+3):(round(Int,tbound2)+2)
            nstatus[G[index1],index3] = "I"
        end

        for index4 in (round(Int,tbound2)+3):(round(Int,endtime)+2)
            nstatus[G[index1],index4] = "R"
        end
    end
    return nstatus
end

# Function to count the occurrence of unique elements in vector sbm.nodes_status,
# in this case, the function counts the incidence of seir nodes
function func_countelements(v)
    u = unique(v)
    d = Dict([(i,count(x->x==i,v)) for i in u])
    D = convert(DataFrame, d)
    return D
end

# Solve the ODE system
#func_beta(beta,par_ode,tspan)

# Compute the parameters for disease properties
infectperiod = Gamma(infectperiod_shape,infectperiod_rate) #Infectious period of the disease
incubperiod = Gamma(incubperiod_shape,incubperiod_rate) #Incubation period of the disease

# Compute the parameters for the clusters and network
communitysize_avg_mat = [(communitysize_avg) for i=1:communitynum]
random_mat = rand(Uniform(-communitysize_range/2,communitysize_range/2), (communitynum, communitynum)) #Obtain a random community size value based on uniform distribution
random_mat = round.(Int, random_mat)
communitysize = communitysize_avg_mat + random_mat
studypopsize = sum(communitysize)
rates_mat = [rates_within + rates_between]
#Cluster = zeros(studypopsize) #Form an arry of the cluster with population size studypopsize

# Main algorithm
for timestep1 in 1:(round(Int,endtime))
    y = func_network(N) #Construct the who-infect-whom network
    G = func_uniqueS(y, nstatus, InfectSize, timestep1) #Make sure the infectees by nodes_names from y are susceptibles
    sbm = func_S2E(G, timestep1) #Spread the disease within the network (from susceptible to exposed)
    nstatus = func_E2I(sbm, G, timestep1) #Allow the disease to progress to infectious and recovered statuses

    #Count the incidence of SEIR
    for timestep2 in 1:(round(Int,endtime))
        D = func_countelements(nstatus[:,timestep2+1])
        D = insertcols!(D, 1, S=0, makeunique=true) #Set S as zero if column does not exist
        D = insertcols!(D, 1, E=0, makeunique=true) #Set E as zero if column does not exist
        D = insertcols!(D, 1, I=0, makeunique=true) #Set I as zero if column does not exist
        D = insertcols!(D, 1, R=0, makeunique=true) #Set R as zero if column does not exist
        D_s = D[1,:S] #Extract the S value from D
        D_e = D[1,:E] #Extract the E value from D
        D_i = D[1,:I] #Extract the E value from D
        D_r = D[1,:R] #Extract the R value from D

        #Put these SEIR incidence values into a DataFrame sbm_sol
        sbm_sol[timestep2+1,:S] = D_s #Put S value into the appropriate cell
        sbm_sol[timestep2+1,:E] = D_e #Put E value into the appropriate cell
        sbm_sol[timestep2+1,:I] = D_i #Put I value into the appropriate cell
        sbm_sol[timestep2+1,:R] = D_r #Put R value into the appropriate cell
        sbm_sol[timestep2+1,:N] = sbm_sol[timestep2+1,:S] + sbm_sol[timestep2+1,:E] + sbm_sol[timestep2+1,:I] + sbm_sol[timestep2+1,:R] #The total number of the population, for accurancy check
    end
    sbm_soln = sbm_sol[sbm_sol[:N].!=0,:] #Remove the 1st row, where N=0
    global G, sbm, nstatus, sbm_soln
end
#println(first(nstatus[11:20],10))
#println(last(sbm_soln,10))

# Write key results into files
CSV.write("./data/sbm_nodes.csv", sbm, writeheader=true)
CSV.write("./data/sbm_sol.csv", sbm_soln, writeheader=true)

#Plot graphs
plot([sbm_soln.S, sbm_soln.E, sbm_soln.I, sbm_soln.R], label=["S","E","I","R"], title="Stochastic Block Model of the Population", xlabel="Time (days)", ylabel="Number of Population")

print("Processing time:")
end #stop timeing the processing time
