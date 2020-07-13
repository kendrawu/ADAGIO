# ADAGIO
A repo holds codes for the project ADAGIO (Adaptive Designs And Genomics In Outbreaks). Note: This is working repository, so code and data are likely to change over time.

Transmission dynamic model: 

sbm.jl (or sbm_fn.jl) is Julia code for simulating a disease transmission using a model that can be described in two parts: 1) a stochastic susceptible-exposed-infectious-removed-vaccinated (SEIRV) model that describes the epidemiological dynamics of the disease within a population; and 2) a network model that describes the spatial dynamics. The network model has three levels of model structure of clusters of varying sizes: i) small clusters, which may represent households or hospital wards; ii) communities of small clusters; and iii) a region of these communities. The rest of the document will assume these small clusters represent households, and each node within each household represent an individual. 

A potential infection occurs within the population if and only if there is a contact edge between two nodes, and there may be a contact edge between two nodes whether or not the nodes are from the same household or community. Rather, whether or not there is a contact edge between two nodes are determined by stochastic block network model (SBM), in which the probability of having contacts between individuals within the same cluster may be different from those that are between clusters. Subsequently, if there is a contact edge between two nodes, and one of them is an infectious (I) individual while the other is a susceptible (S), another stochastic block network model will determine whether there is a transmission edge between these two nodes. 

Infections are also introduced into the population randomly, that is, case importations occur to random clusters at different time points, and the disease importation rate, m, is defined as the number of cases per year arising completely independently from the population being studied. The per-timestep probability of infection for an individual in the population is proportional to the weighted sum of infectious cases in each cluster. Transmission through import case does not require a contact edge between two nodes.

This transmission model is written with an intend that it is versatile enough so that it can be easily modified to describe a specific pathogen in a more specific setting if required. Or else, it can be used for comparing different pathogen(s) and/or different setting(s) using various parameter inputs. In terms of pathogens, this model can be applied to the following three scenarios: a) A pathogen that spreads sustainably between people in the community (e.g., influenza); b) A pathogen that spreads between people largely in hospitals or classrooms (e.g., SARS-CoV); and c) A pathogen that spreads in animal populations, but that can spillover and spread between humans unsustainably or otherwise (e.g., COVID-19, avian influenza, monkeypox).

Clinical trial design:

Using COVID-19 data, preliminary results indicates a cRCT design is more effective in lowering the number of infected cases. For greater details, and the clinical vaccine trial design (WIP), please refer to https://www.overleaf.com/read/zjgfvwqndhfz. NB. Not yet peer reviewed.

DOI: https://doi.org/10.6084/m9.figshare.11777832

Acknowledgment: This study/ project is funded by the National Institute for Health Research (NIHR) (Grant Reference Number PR-OD-1017-20006). The views expressed are those of the author(s) and not necessarily those of the NIHR or the Department of Health and Social Care.
