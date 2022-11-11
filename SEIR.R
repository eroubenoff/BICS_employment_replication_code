library(tidyverse) 
library(reshape2)
library(deSolve)
library(tidycensus)
library(socialmixr)

setwd("~/BICS_employment_replication_code")
source(file.path('utils.R'))


# For u and rho, use Davies values for <18, 35-46, 65+
# unscaled_u = c(0.39, 0.83, 0.74) # susceptibility for each age class
# rho        = c(0.26, 0.36, 0.69) # clinical fraction for each age class

df <- load_bics(wave = 4)$national
pop <- load_pop()
vu_list <- load_vu()
C <- employment_contact_matrix(pop, wave=4, fix_unknowns=TRUE)
baseline_init <- create_SEI2RDV_starting_pop(pop, 
                                             mu = c(3e-5, 0.002, 0.002, 0.069),
                                             u = c(0.39, 0.83, 0.83, 0.74),
                                             rho = c(0.3, 0.5, 0.5, 0.7),
                                             latentPeriod = 3,
                                             infectiousPeriod = 5)

# Alternate specifications of contact matrix (for sensitivity analysis)
C_noreallocation <- employment_contact_matrix(pop, wave = 4, fix_unknowns = FALSE)
C_noschoolclosure <- employment_contact_matrix(pop,  wave = 4, fix_unknowns = TRUE, school_contacts = TRUE) 


########## ---------------------------------------------------------- ##########
# SEI2RDV_ode
# Function to calculate derivatives for lsoda
########## ---------------------------------------------------------- ##########

SEI2RDV_ode = function(time, state, parameters){
  
  
  state_mat <- unlistmat(state)
  
  ## Derivatives
  derivatives <- state_mat 
  derivatives[,] <- 0
  
  # Don't let count drop below 0
  state_mat <- pmax(state_mat, 0)
  
  N <- state_mat["N",] 
  
  with(c(as.list(parameters)), {

    # Vaccine efficacy (first shot and second shot and booster)
    ve <- c(Va = Va, Vb = Vb, Vboost = Vboost)

    lambda <- u*(contact%*%(state_mat["Ic",] + alpha*state_mat["Isc",])/N)
    
    # dS
    derivatives["S", ] = -lambda*state_mat["S",] 
    derivatives["Sx", ] = -lambda*state_mat["Sx",] 
    
    # dVa First shot
    derivatives["Va", ] = -lambda*(1-ve["Va"])*state_mat["Va",] 
    
    # dVb Second shot
    derivatives["Vb", ] =  -lambda*(1-ve["Vb"])*state_mat["Vb",] 
    derivatives["Vbx", ] = -lambda*(1-ve["Vb"])*state_mat["Vbx",] 
    
    # dVboost booster dose
    derivatives["Vboost", ] = -lambda*(1-ve["Vboost"])*state_mat["Vboost",]  
    
    # dE 
    derivatives["E", ] = -sigma*state_mat["E",]  + 
                          lambda*state_mat["S",] + 
                          lambda*state_mat["Sx",] +
                          (1-ve["Va"])*lambda*state_mat["Va",] + 
                          (1-ve["Vb"])*lambda*state_mat["Vb",] + 
                          (1-ve["Vb"])*lambda*state_mat["Vbx",] +
                          (1-ve["Vboost"])*lambda*state_mat["Vboost",]
    
    # dIc
    derivatives["Ic", ]  = -gamma*state_mat["Ic", ] + rho*sigma*state_mat["E",] 
    derivatives["Cc", ]  = rho*sigma*state_mat["E",] 
    
    # dIsc
    derivatives["Isc", ] = -gamma*state_mat["Isc", ] + (1-rho)*sigma*state_mat["E",]
    derivatives["Csc", ] = (1-rho)*sigma*state_mat["E",]

    # dR
    derivatives["R", ]  = (1-mu)*gamma*state_mat["Ic", ] + gamma*state_mat["Isc", ]
    
    # dD
    derivatives["D", ] = mu*gamma*state_mat["Ic", ]
    
    return(list(as.vector(t(derivatives))))
  })
}








########## ---------------------------------------------------------- ##########
# vaccination function
# moves people from S to Va, Va to Vb, or Vb to Vboost depending on rates
########## ---------------------------------------------------------- ##########


vaccination <- function(time, state, parameters){
  state_mat <- unlistmat(state)
  
  boosters <- parameters$boosters
  daily_vax <- parameters$daily_vax
  priority_queue <- parameters$priority_queue
  split_proportion <- parameters$split_proportion
  
  
  if (!boosters) { # If we're not in the booster scenario
    
    # Do first shots
    if (time > parameters$v0) {
      
      vra <- vr_func2(vr = daily_vax/2, vmax = state_mat["S", ], priority_queue = priority_queue,
                      split_proportion = split_proportion)
      
      # Move them from S to Va 
      state_mat["S",] = state_mat["S",] - vra
      state_mat["Va", ] = state_mat["Va",] + vra
    }
    
    
    # Do second shots
    if (time > parameters$v0 + parameters$vt) {
      vrb <- vr_func2(vr = daily_vax/2, vmax = state_mat["Va", ], priority_queue = priority_queue,
                      split_proportion = split_proportion)
      state_mat["Vc", ] = state_mat["Vc", ] + vrb
      
      # Move them from Va to Vb
      state_mat["Va",] = state_mat["Va",] - vrb
      state_mat["Vb", ] = state_mat["Vb",] + vrb
    }
    

    
  } else if(boosters) {
    
    # Of the people waiting to get boosters, distribute 
    vrboost <- vr_func2(vr = daily_vax, vmax = state_mat["Vb", ], priority_queue = priority_queue,
                        split_proportion = split_proportion)
    state_mat["Vc", ] = state_mat["Vc", ] + vrboost
    
    # Move them from Vb to Vboost
    state_mat["Vb", ] = state_mat["Vb", ] - vrboost
    state_mat["Vboost",] = state_mat["Vboost",] + vrboost
    
  }
  
  return(as.vector(t(state_mat)))
  
}






########## ---------------------------------------------------------- ##########
# sim
# runs the set of simulations with supplied parameters
########## ---------------------------------------------------------- ##########


#' Runs a set of simulations with the supplied parameters
#' 
#' @details Uses the parameters passed to construct a set of simulations that either
#' (1) vaccinate no one; (2) prioritize 65+; (3) prioritize workers; (4)
#' split between the two; (5) prioritize low-risk workers; or (6) prioritize 
#' children. Once vaccines are used by the priority group they are split between
#' the remaining groups.
#' 
#' Simulation begins on January 1st 2021 and begins tallying counts 30 days
#' later, at which point the vaccines could start having effect.
#'
#' @param R0 Theoretical R0; contact matrix is scaled accordingly
#' @param infectiousPeriod Duration of infectiousness
#' @param latentPeriod Time spent in E compartment
#' @param rho numeric vector; Proportion of cases in each compartment that are clinical
#' default values from Davies et al 2020
#' @param mu numeric vector; mortality rates
#' @param start_date Simulation start date in "DD-MM-YYYY" format (used for baseline pop)
#' @param daily_vax Number of doses, split between first and second
#' @param v0 Date at which to begin vaccination (integer days before/after start_date)
#' @param vu numeric vector; vaccine uptake rate
#' @param alpha relative transmissibility of asymptomatic vs symptomatic cases
#' @param show_plot Flag to return plots
#' @param init Starting pop; if not passed, will call create_SEI2RDV_starting_pop
#' @param run_all if False, only runs scenarios none, senior, high risk, split
#' @param wave Which wave to run
#' @param method passed to ode
#' @param contact_matrix contact matrix for simulation
#' 
#' @return list of three elements:
#'     out: summary table
#'     plots: list of trajectory plots
#'     trajectories: list of the trajectories
  
sim <- function(
  R0 = 2.5,
  infectiousPeriod = 5, 
  latentPeriod = 3,
  rho = c(0.26, 0.36, 0.36, 0.69),
  mu = c(0.00004, 0.0023, 0.0023, 0.08), 
  start_date = "01-01-2021",
  alpha = 0.5,
  
  v0 = 0, 
  vu = c(0.7, 0.7, 0.7,0.7),
  daily_vax = 2 * 10^6,
  split_proportion = 0.5, 
  
  # Booster scenario
  boosters = FALSE,
  vu_0 = c(0,0,0,0), # Proportion of the pop that starts out fully vaccinated (size of Vb compartment)
  
  # Vaccine efficacy after shots 1, 2, and booster
  Va = 0.8, 
  Vb = 0.9,
  Vboost = 1.0, 

  # Other options
  init = baseline_init,
  run_all =  TRUE,
  wave = 4,
  method = "lsoda",
  contact_matrix =  C, 
  stats_start_date = 1, # Start tallying x days after,
  printout = TRUE
) {
  # susceptibility for each age class
  u_unscaled = c(0.39, 0.83, 0.83, 0.74)

  
  u<- u_unscaled/scale_u_for_R0(wanted_R0 = R0,
                                    u = u_unscaled, 
                                    C = contact_matrix,
                                    gamma = 1/infectiousPeriod,
                                    rho = rho,
                                    alpha = alpha)
  
  
  
  groups <- c("child", "adult_notworking_inperson", "adult_working_inperson", "senior")
  
  
  # Vaccine hesitancy
  if (boosters) {
    # in the booster scenario, start out with a certain proportion of the 
    # susceptibles already fully vaccinated (Vb) 
    init["Vb",] <- round(vu_0*init["S", ])
    
    # Vbx is a compartment of double-vaccinated people denying a booster
    init["Vbx", ] <- round(init["Vb",] * (1-vu))
    init["Vb",] <- init["Vb",] - init["Vbx",]
    
    # Assume the rest are denying the vacccine 
    init["Sx",] <- init["S", ] - (init["Vb",] + init["Vbx",])
    init["S",] <- 0
    
  } else if(!boosters) {
    # Sx is a compartment of people who deny the vaccine 
    init["Sx",] <- init["S",] * (1-vu)
    init["S", ] <- init["S",] - init["Sx", ]
    init["Vbx", ] <- 0
  }
  
  params <- list(
    R0 = R0,
    u= u,
    contact = as.matrix(contact_matrix),
    sigma      = 1/latentPeriod,   
    gamma      = 1/infectiousPeriod,    
    alpha      = alpha, 
    rho        = rho, 
    mu         = mu, 
    v0         = v0,
    vr         = c(0, 0, 0, 0),
    vt         = 25,
    vu         = vu,
    boosters = boosters,
    start_date = start_date,
    daily_vax  = daily_vax,
    Va         = Va,
    Vb         = Vb,
    Vboost     = Vboost,
    method = method
  )
  


  init <- listmat(init)
  

  vaccination_times <- seq(1,364,1)
  params$split_proportion <- c("senior" = split_proportion, "adult_working_inperson"= (1-split_proportion))

  # If we give no one the vaccine
  params$priority_queue <- list(NA)
  trajectory_none<- as.data.frame(lsoda(y = init,
                                        times = seq(from = 0, to = 365, by = 1),
                                        parms = params,
                                        func = SEI2RDV_ode,
                                        events=list(func=vaccination, time=vaccination_times),
                                        atol = 0.01
                                        ))
  # If we give them all to seniors
  params$priority_queue <- list("senior")
  trajectory_senior <- as.data.frame(lsoda(y = init,
                                           times = seq(from = 0, to = 365, by = 1),
                                           parms = params, 
                                           func = SEI2RDV_ode, 
                                           events=list(func=vaccination, time=vaccination_times),
                                           atol = 0.01
  ))
  # If we give them all to high risk 
  params$priority_queue <- list("adult_working_inperson")
  trajectory_hr <- as.data.frame(lsoda(y = init,
                                       times = seq(from = 0, to = 365, by = 1),
                                       parms = params, 
                                       func = SEI2RDV_ode, 
                                       events=list(func=vaccination, time=vaccination_times),
                                       atol = 0.01
  ))
  # If we split them between high risk and seniors
  params$priority_queue <- list(c("senior", "adult_working_inperson"))
  trajectory_split <- as.data.frame(lsoda(y = init,
                                          times = seq(from = 0, to = 365, by = 1),
                                          parms = params, 
                                          func = SEI2RDV_ode, 
                                          events=list(func=vaccination, time=vaccination_times),
                                          atol = 0.01
  ))

  
  # Tiered SR: Seniors, then HC Workers, then split
  params$priority_queue <- list("senior", "adult_working_inperson")
  trajectory_tiered_sr <- as.data.frame(lsoda(y = init,
                                          times = seq(from = 0, to = 365, by = 1),
                                          parms = params, 
                                          func = SEI2RDV_ode, 
                                          events=list(func=vaccination, time=vaccination_times),
                                          atol = 0.01
  ))
  
  # Tiered HC: Seniors, then HC Workers, then split
  params$priority_queue <- list("adult_working_inperson", "senior")
  trajectory_tiered_hc <- as.data.frame(lsoda(y = init,
                                              times = seq(from = 0, to = 365, by = 1),
                                              parms = params, 
                                              func = SEI2RDV_ode, 
                                              events=list(func=vaccination, time=vaccination_times),
                                              atol = 0.01
  ))

  
  trajectories <- list(trajectory_none = trajectory_none,
                       trajectory_senior = trajectory_senior,
                       trajectory_hr = trajectory_hr,
                       trajectory_split = trajectory_split,
                       trajectory_tiered_sr = trajectory_tiered_sr,
                       trajectory_tiered_hc = trajectory_tiered_hc)
  
  
  out <- calculate_SEI2RDV_stats(trajectories, params, stats_start_date = stats_start_date, printout=printout)
  return(list(out = out, trajectories = trajectories))
  
}







# Run a baseline 
s <- sim()


# Run a booster simulation
s_booster <- sim(
  vu_0 = vu_list$vu_jan2022, # Initial fully-vaxxed pop
  R0 = 3.5,
  daily_vax = 1000000,
  
  # Vaccine Efficacy
  Va = 0.5,
  Vb = 0.6,
  Vboost = 0.9,
  boosters = TRUE
)



