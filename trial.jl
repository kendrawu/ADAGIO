# Individualized randomized clinical trial (iRCT)
# Author: Kendra Wu
# Date: 06 February 2020

# Begin timing the processing time
@timev begin #begin timing the processing time

# Introduce packages
using StatsBase # To use sample

# Call parmaeters and functions of the transmission dynamic model
include("sbm8.jl")

# Set parameter values
trial_communitynum = [1] # Community numbers that have been selected for trial
prop_in_trial = 0.8 # Proportion of people in the community who will be enrolled for trial
vac_efficacy = [0.9 0.9] # Vaccine efficacy of vaccine candidates
protection_threshold = 0.8 # Vaccine efficacy above this threshold implies protection

# Initializations
tstatus = fill(-1, N) # (-1=not in trial, 0=in control, 1=vaccine candidate number)
tstatus = hcat([1:1:N;], tstatus) # Put indexes on the 1st column

# Function to divide the population into groups: not in the trial, in the trial (control, vaccine candidates numbers)
function fn_trial_setup(par_disease, tstatus, communitynum_arr, trial_communitynum, nstatus, vac_efficacy, protection_threshold)

    # Inputs:
    # par_disease: User-defined parameters of the disease
    # tstatus: The trial statuses of each individual in the population
    # communitynum_arr: Holds the community number that each individual in the population belongs
    # trial_communitynum: The cluster numbers that are in the trial
    # nstatus: The health statuses of each individual in the population
    # vac_efficacy: Vaccine efficacy of vaccine candidate(s)
    # protection_threshold: A value to determine if a vaccinated person will be infected if exposed
    # timestep: Current time

    # Outputs:
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
        index1=1
        communitysize = communitysize_arr[index1] # Obtain size of community selected for trial
        enrolsize = round(Int,communitysize * prop_in_trial) # Size of enrollment in the community
        potential_nodes_in_trial_index = findall(x->x==index1, tstatus_fn[:,:,index1]) # Find everyone who are in this community

        if (size(potential_nodes_in_trial_index)[1])>0
            for index2 in 2:length(potential_nodes_in_trial_index)
                potential_nodes_in_trial_index_local = potential_nodes_in_trial_index
                push!(potential_nodes_in_trial::Array{Int16,1},potential_nodes_in_trial_index_local[CartesianIndex(index2,1)][1]) # Extract node names from potential_nodes_in_trial_index
            end

            # Select those who will be enrolled into the trial
            nodes_in_trial = sample(potential_nodes_in_trial, enrolsize, replace=false)

            # Assign trial statuses and update tstatus_fn
            for index3 in 1:length(nodes_in_trial)
                tstatus_fn[nodes_in_trial[index3],2] = 0 # In control group
            end
        end
    end

    tstatus_fn = tstatus_fn[:,1:2] # Only keep the first 2 columns
    return tstatus_fn
end

# Main algorithm
tstatus = fn_trial_setup(par_disease, tstatus, communitynum_arr, trial_communitynum, nstatus, vac_efficacy, protection_threshold)
println(tstatus)

print("Processing time:")
end # Stop timeing the processing time
