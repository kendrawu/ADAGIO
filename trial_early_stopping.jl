# Introduce packages
using SpecialFunctions # For generating gamma distribution
using LinearAlgebra
using Distributions
using StatsBase # To use sample
using Statistics
using Optim # To use optimize
using Roots
using CSV
using DataFrames
using DelimitedFiles
using Plots

# Set parameter values
N = 1000 # Total population size
begintime = 1.0
endtime = 200.0 # Duration of simulation (timestep is in a unit of days)
S0 = N # Number of suspectible people in the population at day 0
V0 = 0 # Number of vaccinated individuals in the population at day 0
casenum0 = 5 # Number of infectious individuals introduced into the community to begin the outbreak
immunenum0 = 0 # Number of people who are immune to the disease at the beginning of the outbreak
import_lambda = 0 # Number of occurrences variable for imported cases timeseries, assumed to follow Poisson Distribution
lambda0 = 0.0003 # Per-time-step probability of infection for a susceptible nodes from an infectious neighbour

## The households
# hhnum: Number of households in the network
# hhsize_avg: Average size of one household
# hhsize_range: Range of household sizes (sizes are uniformly distributed)
par_hh = DataFrame(hhnum=240, hhsize_avg=4, hhsize_range=2)
#par_hh = DataFrame(hhnum=1, hhsize_avg=500, hhsize_range=0)

## The clustered network
# communitynum: Number of communities in the clustered network
# communitysize_avg: Average size of one community
# communitysize_range: Range of community sizes (sizes are uniformly distributed)
par_community = DataFrame(communitynum=3, communitysize_avg=300, communitysize_range=10)
#par_community = DataFrame(communitynum=1, communitysize_avg=500, communitysize_range=0)

## Contact and transmission probabilities between nodes
# cprob_hhwithin_cwithin: Probability of contacts of an edge between two nodes in the same household and the same community
# cprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# cprob_hhbetween_cbetween: Probability of transmission of an edge between two nodes in different communities
# tprob_hhwithin_cwithin: Probability of transmission of an edge between two nodes in the same household and the same community
# tprob_hhbetween_cwithin: Probability of contacts of an edge between two nodes in different households but the same community
# tprob_hhbetween_cbetween: Probability of contacts of an edge between two nodes in different communities
par_prob = DataFrame(cprob_hhwithin_cwithin=1, cprob_hhbetween_cwithin=1, cprob_hhbetween_cbetween=1, tprob_hhwithin_cwithin=lambda0, tprob_hhbetween_cwithin=lambda0, tprob_hhbetween_cbetween=lambda0)

## Disease properties
# Use Ebola-like parameters (from Hitchings (2018)) - Gamma-distributed
#par_disease = DataFrame(incubperiod_shape=3.11, incubperiod_rate=0.32, infectperiod_shape=1.13, infectperiod_rate=0.226)
# Use COVID-19-like parameters (from Li et al (2020)) - Gamma-distributed
#par_disease = DataFrame(incubperiod_shape=2.11, incubperiod_rate=0.4, infectperiod_shape=3.0, infectperiod_rate=0.35)
par_disease = DataFrame(incubperiod_shape=3.45, incubperiod_rate=0.66, infectperiod_shape=5.0, infectperiod_rate=0.8)

nsim = 30 # Number of simulations
trial_begintime = 10.0
trial_endtime = endtime
prop_in_trial = 0.6 # Proportion of the population/ cluster will be enrolled in the trial
prop_in_highrisk = 0.2 # Proportion of the population/ in the cluster are at high risk
trial_communitynum = [1]
treatment_gp = 1 # it needs to be an integer, not an array
vac_efficacy = [0.6]
protection_threshold = 0.5
stage = 2 # Number of interim analyses will be done in the trial
allocation_ratio = [0.49 0.52]
#allocation_ratio = [0.9/(0.9+1) 1/(1+0.9)]
#allocation_ratio = [0.78/(0.0018+0.78) 0.0018/(0.0018+0.78)]
alpha = 0.9 # Desired type I error
power = 0.8
#delta = 0.5 # The clinical relevant different to power for
#sigma = 1 # The assumed value for the SD of patient responses
#two_sided = false # A logical variable indicating whether to use a two-sided null hypothesis
p_es = 0.65 # Effect size
p0_te = 0.55 # Treatment effect
e = 1.64 # Efficacy stopping boundary
f = 1.80 # Futility stopping boundary
gamma_infectperiod_duration = rand(Gamma(par_disease[1,:infectperiod_shape],1/par_disease[1,:infectperiod_rate]),1000)
gamma_infectperiod_maxduration = maximum(gamma_infectperiod_duration)
gamma_incubperiod_duration = rand(Gamma(par_disease[1,:incubperiod_shape],1/par_disease[1,:incubperiod_rate]),1000)
gamma_incubperiod_maxduration = maximum(gamma_incubperiod_duration)
method = "iRCT_non_adaptive"

# Initialization
samplesize_truecases = zeros(nsim) # Sample size based on true cases
samplesize_est = zeros(nsim) # Sample size based on estimation
Z = zeros(stage, length(allocation_ratio)) # Wald test statistics for efficacy and futility stopping boundaries

include("sbm_fn.jl") # load all the required functions
include("trial_fn.jl")

if method=="iRCT_MLE"
    MLE = fn_adapt_freq_MLE(par0, method, n_control, n_treatment, n_infectious_control, n_infectious_treatment)
    println("MLE ", MLE)
end

if method=="iRCT_Bayes"
    posterior = fn_adapt_Bayes(nsim)
    println("Posterior: ", posterior)
end

if method=="iRCT_non_adaptive"
    (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
    (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_iRCT(N, tstatus, prop_in_trial, allocation_ratio, vac_efficacy, protection_threshold)
    (nstatus1, tstatus1, soln1, T1, R01) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

    # Determine operation characteristics
    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true1) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize1 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE1 = fn_TTE(nstatus1, tstatus1, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    if nsim==1
        soln_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
        soln_mat[:,1] = soln1[:,1]
        soln_mat[:,2] = soln1[:,2]
        soln_mat[:,3] = soln1[:,3]
        soln_mat[:,4] = soln1[:,4]
        soln_mat[:,5] = soln1[:,5]
        nstatus_mat = nstatus1
        tstatus_mat = tstatus1
        VE_true_mat = VE_true1
        samplesize_mat = samplesize1
        TTE_mat = TTE1
        T = T1
        R0 = R01
    end

    if nsim>=2
        (soln_mat, nstatus_mat, tstatus_mat, VE_true_mat, samplesize_mat, TTE_mat, communitysize_arr, communitynum_arr, T, R0) = fn_iternation_iRCT_non_adpative(nsim, soln1, tstatus1, VE_true1, samplesize1, N, par_hh, par_community, par_prob, par_disease, prop_in_trial, import_lambda, casenum0, immunenum0, allocation_ratio, vac_efficacy, protection_threshold, treatment_gp, gamma_infectperiod_maxduration, trial_begintime, trial_endtime, endtime)
    end
#else
#    throw(ArgumentError("Adaptive method unknown."))
end

if method=="cRCT_non_adaptive"
    (nstatus, tstatus, sbm_sol, hhsize_arr, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc) = fn_pretransmission(N, par_hh, par_community, par_prob, par_disease, import_lambda, casenum0, immunenum0, endtime)
    (nstatus, tstatus, sbm_sol, T, R0) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, begintime+1, trial_begintime, endtime)
    (tstatus, nodes_in_control, nodes_in_treatment) = fn_trialsetup_cRCT(N, par_disease, tstatus, communitysize_arr, communitynum_arr, trial_communitynum, nstatus, allocation_ratio, vac_efficacy, protection_threshold)
    (nstatus1, tstatus1, soln1, T1, R01) = fn_transmodel(nstatus, tstatus, sbm_sol, par_hh, par_community, par_prob, par_disease, hhnum_arr, communitysize_arr, communitynum_arr, importcasenum_timeseries, Gc, trial_begintime+1, endtime, endtime)

    # Determine operation characteristics
    timestep_fn = endtime
    (n_control, n_treatment, n_infectious_control, n_infectious_treatment, n_exposed_control, n_exposed_treatment, VE_true1) = fn_vaccine_efficacy(tstatus, nstatus, timestep_fn, treatment_gp)
    samplesize1 = fn_samplesize_truecases(n_control, n_treatment, n_infectious_control, n_infectious_treatment, treatment_gp, timestep_fn, alpha, power)
    TTE1 = fn_TTE(nstatus1, tstatus1, treatment_gp, trial_begintime, trial_endtime, gamma_infectperiod_maxduration)

    if nsim==1
        soln_mat = zeros(Int, round(Int,endtime), 5) # Same as soln1, except it is a matrix, not a DataFrame
        soln_mat[:,1] = soln1[:,1]
        soln_mat[:,2] = soln1[:,2]
        soln_mat[:,3] = soln1[:,3]
        soln_mat[:,4] = soln1[:,4]
        soln_mat[:,5] = soln1[:,5]
        nstatus_mat = nstatus1
        tstatus_mat = tstatus1
        VE_true_mat = VE_true1
        samplesize_mat = samplesize1
        TTE_mat = TTE1
        T = T1
        R0 = R01
    end

    if nsim>=2
        (soln_mat, nstatus_mat, tstatus_mat, VE_true_mat, samplesize_mat, TTE_mat, communitysize_arr, communitynum_arr, T, R0) = fn_iternation_iRCT_non_adpative(nsim, soln1, tstatus1, VE_true1, samplesize1, N, par_hh, par_community, par_prob, par_disease, prop_in_trial, import_lambda, casenum0, immunenum0, allocation_ratio, vac_efficacy, protection_threshold, treatment_gp, gamma_infectperiod_maxduration, trial_begintime, trial_endtime, endtime)
    end
#else
#    throw(ArgumentError("Adaptive method unknown."))
end


# Results
println("Average sample size (from true cases): ", mean(samplesize_mat))
println("Average number of infectious people in the trial: ", n_infectious_control+n_infectious_treatment)
println("Average number of exposed people in the trial: ", n_exposed_control+n_exposed_treatment)
println("Time-to-Event: ", mean(filter(isfinite, TTE_mat[:,2]))
println("Vaccine efficacy: ", mean(filter(isfinite, VE_true_mat))

Y = fn_divide(soln_mat, endtime, nsim, 3)
lowerCI = zeros(round(Int,endtime))
upperCI = zeros(round(Int,endtime))
medianCI = zeros(round(Int,endtime))
for index = 1:round(Int,endtime)
    lowerCI[index] = quantile!(Y[index,:], 0.05)
    upperCI[index] = quantile!(Y[index,:], 0.95)
    medianCI[index] = median(Y[index,:])
end

X = fn_translate_nstatus(N, nstatus_mat, tstatus_mat, endtime, nsim, round(Int,endtime+1))
#infectednum = zeros(Int, nsim)
#for index in 1: nsim
#    infectednum[index] = fn_casecounting(X, N, prop_in_trial, allocation_ratio, index)
#end
#n_infected = sum(infectednum)
#println("Number of infected on average: ", n_infected/nsim)

(plot1, plot2) = fn_plot(X,Y)
plot(plot1,plot2,layout=(2,1))
