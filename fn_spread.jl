# Author: Kendra Wu
# Date: 02 February 2020

using SpecialFunctions
using Distributions
using StatsBase

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

    # Determine the contacts of infector(s) i
    for index1 in 1:(size(i,1))
        contact_map_i = contact_map[i[index1],:] # Extract the contact map out for i[index1]
        sizeend = size(contact_map_i,1) # Find size of contacts of i[index1]
        node_location = contact_map_i[2:sizeend,:] # Remove i[index1] from column 1
        for index2 in 1:(size(node_location,1))
            if node_location[index2] !=0
                push!(contacts_list,index2) # Determine who are the contacts of i[index1] by removing the zeros
            end
        end

        # Determine the health status of normal contacts, neighhours, and high-risk individuals of i[index1]
        susceptible_potential_infectee = intersect(contacts_list,s) # The overall susceptible list
        susceptible_potential_infectee_neighbour = intersect(susceptible_potential_infectee, neighbour_list) # For the neighbours
        susceptible_potential_infectee_highrisk = intersect(susceptible_potential_infectee, highrisk_list) # For high-risk groups
        susceptible_potential_infectee_normal = setdiff(setdiff(susceptible_potential_infectee, susceptible_potential_infectee_neighbour), highrisk_list) # For normal contacts
        vaccinated_potential_infectee = intersect(contacts_list,v) # The overall vaccinated list
        vaccinated_potential_infectee_neighbour = intersect(vaccinated_potential_infectee, neighbour_list) # For the neighbours
        vaccinated_potential_infectee_highrisk = intersect(vaccinated_potential_infectee, highrisk_list) # For high-risk groups
        vaccinated_potential_infectee_normal = setdiff(setdiff(vaccinated_potential_infectee, vaccindated_potential_infectee_neighbour), highrisk_list) # For normal contacts

        # Determine size of exposure among the susceptible individuals
        susceptible_potential_exposed_size = size(susceptible_potential_infectee,1) # The overall size of the exposed among the susceptible
        susceptible_exposed_size_normal = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob_s), 1)[1]) # Normal contacts
        susceptible_exposed_size_neighbour = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob * neighbour_scalar_s), 1)[1]) # The neighbours
        susceptible_exposed_size_highrisk = round(Int, rand(Poisson(susceptible_potential_exposed_size * infect_prob * high_risk_scalar_s), 1)[1]) # The high risks

        # Determine who will be exposed among the susceptible who are normal contacts, neighbours, and high-risk groups
        susceptible_exposed_size_normal_min = min(size(susceptible_potential_infectee_normal,1), susceptible_exposed_size_normal) # Make sure to have enough number of people to draw from among the normal contacts
        susceptible_exposed_size_neighbour_min = min(size(susceptible_potential_infectee_neighbour,1), susceptible_exposed_size_neighbour) # Make sure to have enough number of people to draw from among the neighbours
        susceptible_exposed_size_highrisk_min = min(size(susceptible_potential_infectee_highrisk,1), susceptible_exposed_size_highrisk) # Make sure to have enough number of people to draw from among the high risk groups
        susceptible_infectee_normal = sample(susceptible_potential_infectee_normal, susceptible_exposed_size_normal_min, replace=false) # Normal contacts
        susceptible_infectee_neighbour = sample(susceptible_potential_infectee_neighbour, susceptible_exposed_size_neighbour_min, replace=false) # The neighbours
        susceptible_infectee_highrisk = sample(susceptible_potential_infectee_highrisk, susceptible_exposed_size_highrisk_min, replace=false) # The high risks
        susceptible_infectee = union(susceptible_infectee_normal, susceptible_infectee_neighbour, susceptible_infectee_highrisk) # Put the lists together
        susceptible_infectee = sort(susceptible_infectee) # Sort the list

        # Determine size of exposure among the vaccinated individuals
        vaccinated_potential_exposed_size = size(vaccindated_potential_infectee,1) # The overall size of the exposed among the susceptible
        vaccinated_exposed_size_normal = round(Int, vaccinated_potential_exposed_size * infect_prob_v) # Normal contacts
        vaccinated_exposed_size_neighbour = round(Int, vaccinated_potential_exposed_size * infect_prob * neighbour_scalar_v) # The neighbours
        vaccinated_exposed_size_highrisk = round(Int, vaccinated_potential_exposed_size * infect_prob * high_risk_scalar_v) # The high risks

        # Determine who will be exposed among the vaccinated who are normal contacts, neighbours, and high-risk groups
        vaccinated_exposed_size_normal_min = min(size(vaccinated_potential_infectee_normal,1), vaccinated_exposed_size_normal) # Make sure to have enough number of people to draw from among the normal contacts
        vaccinated_exposed_size_neighbour_min = min(size(vaccinated_potential_infectee_neighbour,1), vaccinated_exposed_size_neighbour) # Make sure to have enough number of people to draw from among the neighbours
        vaccinated_exposed_size_highrisk_min = min(size(vaccinated_potential_infectee_highrisk,1), vaccinated_exposed_size_highrisk) # Make sure to have enough number of people to draw from among the high risk groups
        vaccinated_infectee_normal = sample(vaccinated_potential_infectee_normal, vaccinated_exposed_size_normal, replace=false) # The normal contacts
        vaccinated_infectee_neighbour = sample(vaccinated_potential_infectee_neighbour, vaccinated_exposed_size_neighbour, replace=false) # The neighbours
        vaccinated_infectee_highrisk = sample(vaccinated_potential_infectee_highrisk, vaccinated_exposed_size_highrisk, replace=false) # The high risk groups
        vaccinated_infectee = union(vaccinated_infectee_normal, vaccinated_infectee_neighbour, vaccinated_infectee_highrisk) # Put the lists together
        vaccinated_infectee = sort(vaccinated_infectee) # Sort the list

        # Move s that is contact of i to e
        if size(susceptible_infectee,1)>0
            e = vcat(e, susceptible_infectee) # Put the list of susceptible infectee(s) into the original e

            for index3 in 1:(size(susceptible_infectee,1))

                # Add info of those exposed onto exposed_days
                exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                exposed_days[1,1] = susceptible_infectee[index3] # The node name
                exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1]) # Incubation period

                # Add info of those exposed onto incubperiod_info
                incubperiod_info[susceptible_infectee[index3],2] = timestep # Current time
                incubperiod_info[susceptible_infectee[index3],3] = timestep + exposed_days[2,1] # End time of incubation period
            end
        end

        # Move v that has protection below protection_threshold and is a contact of i to e
        if size(vaccindated_infectee,1)>0
            for index4 in 1:(size(vaccindated_infectee,1))
                if vac_efficacy[vaccindated_infectee[index4]] < protection_threshold # Check if vaccine efficacy of the vaccinated is below threshold
                    push!(e::Array{Int,1},vaccindated_infectee[index4]) # Put node name into e

                    # Add info of those exposed onto exposed_days
                    global exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                    exposed_days[1,1] = vaccindated_infectee[index4] # The node name
                    exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1]) # Incubation period

                    # Add info of those exposed onto incubperiod_info
                    incubperiod_info[vaccindated_infectee[index4],2] = timestep # Current time
                    incubperiod_info[vaccindated_infectee[index4],3] = timestep + exposed_days[2,1] # End time of incubation period
                end
            end
        end

    end

    e = sort(e) # sort e in order of node names
    exposed_days = reshape(exposed_days, 2, :)' # Convert exposed_days from array to matrix

    return e, exposed_days, incubperiod_info
end
