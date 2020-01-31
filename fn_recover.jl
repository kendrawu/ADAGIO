# Author: Kendra Wu
# Date: 31 January 2020

using SpecialFunctions
using Distributions
using StatsBase

# Function to move everyone in e forward by one timestep
function fn_recover(e, i, infectperiod, hosprate, enrol_nodes, timestep)

    # Inputs:
    # e: node names of exposed individuals
    # i: node names of infectious individuals
    # infectperiod: Infectious period distribution parameters
    # hosprate: Hospitalization rate distribution parameters
    # enrol_nodes: node names of those enrolled
    # timestep: Current time

    # Outputs:
    # e: node names of exposed individuals
    # i: node names of infectious individuals

    for index in 1:(size(e,2))
        i_enrolled = findall(x->x==e[index], enrol_nodes) # check if e[i] is in enrol_nodes
        if size(i_enrolled)[1]!=0
            node_name = e[index] # Obtain node_name that needs to progress from e
            Te[index] = Te[index]+1 # Obtain time of e and move e one timestep, assume e is not infectious
            Ti[index] = Ti[index]+round.(rand(Gamma(infectperiod[1]/infectperiod[2]),1)[1])+2 # The same node move from e to i
        end
    end

    hospad_prop = size(i,1) * hosprate # A proportion of randomly chosen i are admitted into hospital
    h = sample(1:(size(i,1)), round(Int, hospad_prop), replace=false) # Obtain a random sample of hospad_prop

    return e, i
end
