# ADAGIO
A repo holds codes for the project ADAGIO (Adaptive Designs And Genomics In Outbreaks)

Transmission dynamic model: 

sbm.jl is a Julia code for simulating a transmission disease model that can be described in two parts: 1) a stochastic susceptible-exposed-infectious-removed-vaccinated (SEIRV) model that describes the epidemiological dynamics of the disease within a population; and 2) a network model that describes the spatial dynamics. The network model has three levels of model structure of clusters of varying sizes: i) small clusters, which may represent households or hospital wards; ii) communities of small clusters; and iii) a region of these communities. The rest of the document will assume these small clusters represent households, and each node within each household represent an individual. 

A potential infection occurs within the population if and only if there is a contact edge between two nodes, and there may be a contact edge between two nodes whether or not the nodes are from the same household or community. Rather, whether or not there is a contact edge between two nodes are determined by stochastic block network model (SBM), in which the probability of having contacts between individuals within the same cluster may be different from those that are between clusters. Subsequently, if there is a contact edge between two nodes, and one of them is an infectious (I) individual while the other is a susceptible (S), another stochastic block network model will determine whether there is a transmission edge between these two nodes. 

Infections are also introduced into the population randomly, that is, case importations occur to random clusters at different time points, and the disease importation rate, m, is defined as the number of cases per year arising completely independently from the population being studied. The per-timestep probability of infection for an individual in the population is proportional to the weighted sum of infectious cases in each cluster. Transmission through import case does not require a contact edge between two nodes.

DOI: 10.6084/m9.figshare.11777832
