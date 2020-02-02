# Author: Kendra Wu
# Date: 31 January 2020

using SpecialFunctions
using Distributions
using StatsBase

N = 20
contact_map = [1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               2 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               3 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               4 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               5 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1;
               6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
               20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
i = [5]
s = [1 3 8 9 10 11 12 13 14 15]
e = [6]
v = [2 4 7 16 17 18 19 20]
neighbour_list = [3 9 10 11 12 14 15 17 20]
highrisk_list = [4 8 13 16 18 19]
vac_efficacy = [0 0.7 0 0.7 0 0 0.8 0 0 0 0 0 0 0 0 0.8 0.8 0.9 0.9 0.9]
protection_threshold = 0.8
infect_prob_s = 0.9
infect_prob_v = 0.9
incubperiod = [3.11 0.32]
incubperiod_info = [1 0 0.0; 2 0 0.0; 3 0 0.0; 4 0 0.0; 5 0 0.0;
                    6 1 2.5; 7 0 0.0; 8 0 0.0; 9 0 0.0; 10 0 0.0;
                    11 0 0.0;12 0 0.0;13 0 0.0;14 0 0.0;15 0 0.0;
                    16 0 0.0; 17 0 0.0; 18 0 0.0; 19 0 0.0; 20 0 0.0]
timestep = 2
neighbour_scalar_s = 0.35
neighbour_scalar_v = 0.35
highrisk_scalar_s = 2.5
highrisk_scalar_v = 2.5

# Function to move nodes s and v that are contacts of node i to nodes e
function fn_spread(s, i, v, infect_prob, incubperiod, incubperiod_info, vac_efficacy, protection_threshold, neighbour_list, highrisk_list, neighbour_scalar_s, neighbour_scalar_v, highrisk_scalar_s, highrisk_scalar_v, contact_map, timestep)

    # Inputs:
    # s: Node names of susceptible individuals
    # i: Node names of infectious individuals
    # v: Node names of vaccinated individuals
    # infect_prob: Probability of infecting a normal contact in one day among susceptible and vaccinated
    # incubperiod: Array that holds parameters of incubation period
    # incubperiod_info: Matrix that holds node names, time reference, and end time of incubation period
    # vac_efficacy: Vaccine efficacy of vaccine candidate(s)
    # protection_threshold: A value to determine if a vaccinated person will be infected if exposed
    # neighbour_list: A list of neighour(s) of infectors
    # highrisk_list: A list of high-risk individuals
    # neighbour_scalar: Scalar of neighhours among susceptible and vaccinated
    # high_risk_scalar: Scalar of high-risk individuals among susceptible and vaccinated
    # contact_map: A list of contacts of infectors
    # timestep: Current time

    # Outputs:
    # e: Node names of exposed individuals
    # exposed_days: Number of days nodes e are exposed
    # incubperiod_info: Matrix that holds node names, time reference, and end time of incubation period

    # Assign type to variables
    contacts_list = Int[]
    susceptible_potential_infectee = Int[]
    susceptible_potential_infectee_neighbour = Int[]
    susceptible_potential_infectee_highrisk = Int[]
    vaccindated_potential_infectee = Int[]
    vaccindated_potential_infectee_neighbour = Int[]
    vaccindated_potential_infectee_highrisk = Int[]
    exposed_days = Float16[]
    global e
    global exposed_days

    # Determine the contacts of infectors
    for index1 in 1:(size(i,1))
        contact_map_i = contact_map[i[index1],:]
        sizeend = size(contact_map_i,1)
        node_location = contact_map_i[2:sizeend,:]
        for index2 in 1:(size(node_location,1))
            if node_location[index2] !=0
                push!(contacts_list,index2)
            end
        end

        # Determine the health status of normal contacts, neighhours, and high-risk individuals of i[index1]
        susceptible_potential_infectee = intersect(contacts_list,s) # The overall susceptible list
        susceptible_potential_infectee_neighbour = intersect(susceptible_potential_infectee, neighbour_list)
        susceptible_potential_infectee_highrisk = intersect(susceptible_potential_infectee, highrisk_list)
        susceptible_potential_infectee_normal = setdiff(setdiff(susceptible_potential_infectee, susceptible_potential_infectee_neighbour), highrisk_list)
        vaccinated_potential_infectee = intersect(contacts_list,v) # The overall vaccinated list
        vaccinated_potential_infectee_neighbour = intersect(vaccindated_potential_infectee, neighbour_list)
        vaccinated_potential_infectee_highrisk = intersect(vaccindated_potential_infectee, highrisk_list)
        vaccinated_potential_infectee_normal = setdiff(setdiff(vaccindated_potential_infectee, vaccindated_potential_infectee_neighbour), highrisk_list)

        # Determine size of exposure among the susceptible individuals
        susceptible_potential_exposed_size = size(susceptible_potential_infectee,1)
        susceptible_exposed_size_normal = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob_s), 1)[1]) # for normal contact
        susceptible_exposed_size_neighbour = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob * neighbour_scalar_s), 1)[1]) # for neighbours
        susceptible_exposed_size_highrisk = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob * high_risk_scalar_s), 1)[1]) # for high risk groups

        # Determine who will be exposed among the susceptible who are normal contacts, neighbours, and high-risk groups
        susceptible_exposed_size_normal_min = min(size(susceptible_potential_infectee_normal,1), susceptible_exposed_size_normal) # Make sure to have enough number of people to draw from among the normal contacts
        susceptible_exposed_size_neighbour_min = min(size(susceptible_potential_infectee_neighbour,1), susceptible_exposed_size_neighbour) # Make sure to have enough number of people to draw from among the neighbours
        susceptible_exposed_size_highrisk_min = min(size(susceptible_potential_infectee_highrisk,1), susceptible_exposed_size_highrisk) # Make sure to have enough number of people to draw from among the high risk groups
        susceptible_infectee_normal = sample(susceptible_potential_infectee_normal, susceptible_exposed_size_normal_min, replace=false) # For normal contacts
        susceptible_infectee_neighbour = sample(susceptible_potential_infectee_neighbour, susceptible_exposed_size_neighbour_min, replace=false) # For neighbours
        susceptible_infectee_highrisk = sample(susceptible_potential_infectee_highrisk, susceptible_exposed_size_highrisk_min, replace=false) # For high risk groups
        susceptible_infectee = union(susceptible_infectee_normal, susceptible_infectee_neighbour, susceptible_infectee_highrisk) # Put the lists together
        susceptible_infectee = sort(susceptible_infectee)

        # Determine size of exposure among the vaccinated individuals
        vaccinated_potential_exposed_size = size(vaccindated_potential_infectee,1)
        vaccinated_exposed_size_normal = round(Int, vaccinated_potential_exposed_size * infect_prob_v) # for normal contact
        vaccinated_exposed_size_neighbour = round(Int, vaccinated_potential_exposed_size * infect_prob * neighbour_scalar_v) # for neighbours
        vaccinated_exposed_size_highrisk = round(Int, vaccinated_potential_exposed_size * infect_prob * high_risk_scalar_v) # for high risk groups

        # Determine who will be exposed among the vaccinated who are normal contacts, neighbours, and high-risk groups
        vaccinated_exposed_size_normal_min = min(size(vaccinated_potential_infectee_normal,1), vaccinated_exposed_size_normal) # Make sure to have enough number of people to draw from among the normal contacts
        vaccinated_exposed_size_neighbour_min = min(size(vaccinated_potential_infectee_neighbour,1), vaccinated_exposed_size_neighbour) # Make sure to have enough number of people to draw from among the neighbours
        vaccinated_exposed_size_highrisk_min = min(size(vaccinated_potential_infectee_highrisk,1), vaccinated_exposed_size_highrisk) # Make sure to have enough number of people to draw from among the high risk groups
        vaccinated_infectee_normal = sample(vaccinated_potential_infectee_normal, vaccinated_exposed_size_normal, replace=false) # For normal contacts
        vaccinated_infectee_neighbour = sample(vaccinated_potential_infectee_neighbour, vaccinated_exposed_size_neighbour, replace=false) # For neighbours
        vaccinated_infectee_highrisk = sample(vaccinated_potential_infectee_highrisk, vaccinated_exposed_size_highrisk, replace=false) # For high risk groups
        vaccinated_infectee = union(vaccinated_infectee_normal, vaccinated_infectee_neighbour, vaccinated_infectee_highrisk) # Put the lists together
        vaccinated_infectee = sort(vaccinated_infectee)

        # Move s that is contact of i to e
        if size(susceptible_infectee,1)>0
            e = vcat(e, susceptible_infectee)

            for index3 in 1:(size(susceptible_infectee,1))

                # Add info of those exposed onto exposed_days
                exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                exposed_days[1,1] = susceptible_infectee[index3]
                exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1])

                # Add info of those exposed onto incubperiod_info
                incubperiod_info[susceptible_infectee[index3],2] = timestep
                incubperiod_info[susceptible_infectee[index3],3] = timestep + exposed_days[2,1]
            end
        end

        # Move v that has protection below protection_threshold and is a contact of i to e
        if size(vaccindated_infectee,1)>0
            #println(vaccindated_potential_infectee)
            for index4 in 1:(size(vaccindated_infectee,1))
                #println("vaccine efficacy = ", vac_efficacy[vaccindated_potential_infectee[index3]])
                if vac_efficacy[vaccindated_infectee[index4]] < protection_threshold
                    push!(e::Array{Int,1},vaccindated_infectee[index4])

                    # Add info of those exposed onto exposed_days
                    global exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                    exposed_days[1,1] = vaccindated_infectee[index4]
                    exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1])

                    # Add info of those exposed onto incubperiod_info
                    incubperiod_info[vaccindated_infectee[index4],2] = timestep
                    incubperiod_info[vaccindated_infectee[index4],3] = timestep + exposed_days[2,1]
                end
            end
        end

    end

    e = sort(e) # sort e in order of node names
    exposed_days = reshape(exposed_days, 2, :)' # Convert exposed_days from array to matrix
    
    return e, exposed_days, incubperiod_info
end

(e, exposed_days, incubperiod_info) = fn_spread(s, i, v, infect_prob, incubperiod, incubperiod_info, vac_efficacy, protection_threshold, neighbour_list, highrisk_list, neighbour_scalar_s, neighbour_scalar_v, highrisk_scalar_s, highrisk_scalar_v, contact_map, timestep)
