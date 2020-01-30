# Author: Kendra Wu
# Date: 31 January 2020

# Function to move everyone in e forward by one timestep
function fn_recover(e, i, infectperiod, hosprate, enrol_nodes, enroll_day, timestep)

    # Inputs:
    # e: node names of exposed individuals (vector of integers)
    # i: node names of infectious individuals (vector of integers)
    # infectperiod: Infectious period distribution parameters (two double parmeters)
    # hosprate: Hospitalization rate distribution parameters (two double parameters)
    # enrol_nodes: node names of those enrolled (vector of integers)
    # enroll_day: Day of enrollment for the node names of those enrolled
    # timestep: Current time

    # Outputs:
    # e: node names of exposed individuals
    # i: node names of infectious individuals

    for i in 1:(size(e,1))
        i_enrolled = findall(x->x==e[i], enrol_nodes) # check if e[i] is in enrol_nodes
        if size(i_enrolled)[1]!=0
            node_name = e[i] # Obtain node_name that needs to progress from e
            e[i, timestep] = e[i, timestep+1] # Move e one timestep, assume e is not infectious
            i[i, timestep+2] = i[i, timestep+infectperiod+2] # The same node move from e to i
        end
    end

    hospad_prop = i * hosprate
    h = sample(1:(size(i,1)), hospad_prop, replace=false) # A proportion of randomly chosen i are admitted into hospital

    return e, i
end

(e, i) = fn_recover(e, i, infectperiod, hosprate, enrol_nodes, enroll_day, timestep)
