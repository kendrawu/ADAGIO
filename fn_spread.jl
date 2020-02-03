# Author: Kendra Wu
# Date: 02 February 2020

using SpecialFunctions
using Distributions
using StatsBase

# Function to move nodes s and v that are contacts of node i to nodes e
function fn_spread(s, i, v, infect_prob, incubperiod, incubperiod_info, vac_efficacy, protection_threshold, neighbour_list, highrisk_list, neighbour_scalar, highrisk_scalar, contact_map, timestep)

    # Inputs:
    # s: Node names of susceptible individuals
    # i: Node name of infectious individual
    # v: Node names of vaccinated individuals
    # infect_prob: Probability of infecting a normal contact in one day among susceptible and vaccinated
    # incubperiod: Array that holds parameters of incubation period
    # incubperiod_info: Matrix that holds node names, time reference, and end time of incubation period
    # vac_efficacy: Vaccine efficacy of vaccine candidate(s)
    # protection_threshold: A value to determine if a vaccinated person will be infected if exposed
    # neighbour_list: A list of neighour(s) of infectors
    # highrisk_list: A list of high-risk individuals
    # neighbour_scalar: Scalar of neighhours among susceptible and vaccinated
    # highrisk_scalar: Scalar of high-risk individuals among susceptible and vaccinated
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
    vaccinated_potential_infectee = Int[]
    vaccinated_potential_infectee_neighbour = Int[]
    vaccinated_potential_infectee_highrisk = Int[]
    exposed_days = Float16[]
    global e
    global exposed_days

    # Initialization
    v_new = v

    # Determine the contacts of one infector
    for index1 in 1:(size(i,2))
        contact_map_i = contact_map[i[index1],:] # Extract the contact map out for i
        sizeend = size(contact_map_i,1) # Find size of contacts of i
        node_location = contact_map_i[2:sizeend] # Remove i from column 1
        for index2 in 1:(size(node_location,1))
            if node_location[index2] !=0
                push!(contacts_list,index2) # Determine who are the contacts of i[index1] by removing the zeros
            end
        end

        # Determine the health status of normal contacts, neighhours, and high-risk individuals of i[index1]
        susceptible_contacts = intersect(contacts_list, s) # The susceptible contacts
        vaccinated_contacts = intersect(contacts_list, v_new) # The vaccinated contacts
        icontacts = union(susceptible_contacts, vaccinated_contacts) # The potential infectees of i[index1]: the susceptible and vaccinated contacts of i[index1]
        icontacts_neighbour = intersect(icontacts, neighbour_list) # The neighbours
        icontacts_highrisk = intersect(icontacts, highrisk_list) # The highrisks
        icontacts_normal = setdiff(setdiff(icontacts, icontacts_neighbour), icontacts_highrisk) # Normal contacts

        # Determine size of exposure
        exposed_size_normal = size(i,2) * round(Int, rand(Poisson(infect_prob), 1)[1]) # Normal contacts
        exposed_size_neighbour = size(i,2) * round(Int, rand(Poisson(infect_prob * neighbour_scalar), 1)[1]) # The neighbours
        exposed_size_highrisk = size(i,2) * round(Int, rand(Poisson(infect_prob * highrisk_scalar), 1)[1]) # The high risks

        # Take into account susceptible depletion
        exposed_size_normal_min = min(size(icontacts_normal,1), exposed_size_normal)
        exposed_size_neighbour_min = min(size(icontacts_neighbour,1), exposed_size_neighbour)
        exposed_size_highrisk_min = min(size(icontacts_highrisk,1), exposed_size_highrisk)

        # Determine who will be exposed among normal contacts, neighbours, and high-risk groups by sampling
        exposed_normal = sample(icontacts_normal, exposed_size_normal_min, replace=false) # The normal contacts
        exposed_neighbour = sample(icontacts_neighbour, exposed_size_neighbour_min, replace=false) # The neighbours
        exposed_highrisk = sample(icontacts_highrisk, exposed_size_highrisk_min, replace=false) # The high risk groups
        exposed = union(exposed_normal, exposed_neighbour, exposed_highrisk) # Put the lists together

        # Separate the potential exposed who are currently susceptible and vaccinated
        susceptible_infectee = intersect(exposed, s) # The potential exposed that are currently susceptible
        vaccinated_infectee = intersect(exposed, v_new) # The potential exposed that are currently vaccinated
        susceptible_infectee = sort(susceptible_infectee)
        vaccinated_infectee = sort(vaccinated_infectee)

        # Move s that is contact of i to e
        if size(susceptible_infectee,1)>0
            e = vcat(e, susceptible_infectee) # Move s that is contact of i[index1] to e
            s = setdiff(s, e) # update s

            for index3 in 1:(size(susceptible_infectee,1))

                # Add info of those exposed onto exposed_days
                exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                exposed_days[1,1] = susceptible_infectee[index3] # Node name
                exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1]) # Incubation period

                # Add info of those exposed onto incubperiod_info
                incubperiod_info[susceptible_infectee[index3],2] = timestep # Current time
                incubperiod_info[susceptible_infectee[index3],3] = timestep + exposed_days[2,1] # End time of incubation period
            end
        end

        # Move v that has protection below protection_threshold and is a contact of i to e
        if size(vaccinated_infectee,1)>0
            for index4 in 1:(size(vaccinated_infectee,1))
                if vac_efficacy[vaccinated_infectee[index4]] < protection_threshold # Check if vaccine efficacy is below threshold
                    push!(e::Array{Int,1},vaccinated_infectee[index4]) # If so, add node name to e
                    v_new = setdiff(v_new, vaccinated_infectee[index4]) # remove e from v_new

                    # Add info of those exposed onto exposed_days
                    exposed_days = [[0 0]'; exposed_days] # Add a new row to exposed_days
                    exposed_days[1,1] = vaccinated_infectee[index4] # Node name
                    exposed_days[2,1] = round.(rand(Gamma(incubperiod[1]/incubperiod[2]),1)[1]) # Incubation period

                    # Add info of those exposed onto incubperiod_info
                    incubperiod_info[vaccinated_infectee[index4],2] = timestep # Current time
                    incubperiod_info[vaccinated_infectee[index4],3] = timestep + exposed_days[2,1] # End time of incubation period
                end
            end
        end
    end

    e = sort(e) # sort e in order of node names
    exposed_days = reshape(exposed_days, 2, :)' # Convert exposed_days from array to matrix

    return e, exposed_days, incubperiod_info
end
