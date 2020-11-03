#####################################
## Title: Run RegrImp simulation
## Author: Rose Sisk
## Date Created: 01/04/2020
## Date edited: 16/06/2020
#######################################
source("simFunctions.R") # main functions to run simulation

library(tidyverse)
set.seed(147)

args <- commandArgs(trailingOnly = T)
s <- as.numeric(args[1])## pull in task ID from bash script. Use this as an index for the params dataset.

grd4 = c(-0.69, 0, 0.5, 1.1)
grd2 = c(0.1, 1)


### set parameter combinations.
params <- crossing(n = 10000,
                   pi_Y = 0.1,
                   pi_R1 = c(0.1, 0.25, 0.5, 0.75),
                   sigma_X1 = 1,
                   sigma_X2 = 1,
                   alpha_0 = 0,
                   eta_0 = 0,
                   beta_0 = 0,
                   gamma_0 = 0,
                   beta_X1 = grd4,
                   beta_X2 = grd4,
                   gamma_X1 = grd2,
                   gamma_X2 = grd2,
                   gamma_X1X2 = c(0,1))


## Edit the intercepts so that the output corresponds to the required missing level/outcome prevalence. 
params <- params %>%
  mutate(beta_0 = log(pi_R1/(1-pi_R1)),
         gamma_0 = log(pi_Y/(1-pi_Y)))

params <- as.matrix(params)


#plan(multiprocess)

###stores results alongside parameter values in a nested dataset (retaining names of outputs)
#sims <- params %>%
 # mutate(results = future_pmap(list(n, pi_Y, pi_R1, sigma_X1, sigma_X2, 
  #                                  alpha_0, eta_0, beta_0, gamma_0, 
   #                                 beta_X1, beta_X2, gamma_X1, gamma_X2, nrun),
    #                           simulate_nrun))

### try some tests to see if/how this speeds things up.


#system.time({sims <- params[1:50,] %>%
 # mutate(results = future_pmap(list(n, pi_Y, pi_R1, sigma_X1, sigma_X2, 
  #                                  alpha_0, eta_0, beta_0, gamma_0, 
   #                                 beta_X1, beta_X2, gamma_X1, gamma_X2, nrun),
    #                           simulate_nrun))
#})

simulate_nrun(n = params[s, 1], 
              pi_Y = params[s, 2], 
              pi_R1 = params[s, 3], 
              sigma_X1 = params[s, 4], 
              sigma_X2 = params[s, 5], 
              alpha_0 = params[s, 6],  
              eta_0 = params[s, 7],
              beta_0 = params[s, 8],
              gamma_0 = params[s, 9],
              beta_X1 = params[s, 10],
              beta_X2 = params[s, 11],
              gamma_X1 = params[s, 12],
              gamma_X2 = params[s, 13],
              gamma_X1X2 = params[s, 14],
              nrun = 200)











