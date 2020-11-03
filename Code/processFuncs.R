#################################################################################
## Regression Imputation vs Multiple Imputation for prediction: a simulation
## Set of helper functions to process results (RDS files) from CSF run of simulation
#################################################################################

library(tidyverse)
library(purrr)
library(ggpubr)
library(wesanderson)
library(here)

## create a tibble to store results. Have one col to record each of the variable 
##parameters. Then one col per performance measure/setting.
grd4 = c(-0.69, 0, 0.5, 1.1)
grd2 = c(0.1, 1)

params <- crossing(pi_R1 = c(0.1, 0.25, 0.5, 0.75),
                   beta_X1 = grd4,
                   beta_X2 = grd4,
                   gamma_X1 = grd2,
                   gamma_X2 = grd2,
                   gamma_X1X2 = c(0,1))

## function that takes parameter values as inputs, then reads in and averages the corresponding results.
empSE <- function(x, nsim = 200) {
  mn <- mean(x)
  sqrt(
    (1/(nsim - 1))*sum((x-mn)^2)
  )/sqrt(nsim)
}

CI_emp <- function(x, upper = T){
  SE <- empSE(x)
  if (upper == T){
    CI <- mean(x) + 1.96*SE
  } else if (upper == F){
    CI <- mean(x) - 1.96*SE
  }
  CI
}

load_res <- function(pi_R1 = 0.1,
                     beta_X1 = -0.69,
                     beta_X2 = -0.69,
                     gamma_X1 = 0.1,
                     gamma_X2 = 0.1,
                     gamma_X1X2 = 0) {
  
  flnm <- paste0("resMat", "pi_R1", pi_R1, "beta_X1", beta_X1, "beta_X2", beta_X2, "gamma_X1", gamma_X1, "gamma_X2", gamma_X2, "gamma_X1X2", gamma_X1X2, ".RDS")
  dat <- readRDS(here("Outputs", flnm))
  dat <- as.tibble(dat)
  
  outp <- dat %>%
    summarise(across(everything(), 
                     list(mean = mean, SE = ~empSE(.x, 200), CI_L = ~CI_emp(.x, upper = F), CI_U = ~CI_emp(.x, upper = T))))
  
  outp
  
}

load_raw <- function(pi_R1 = 0.1,
                     beta_X1 = 1.1,
                     beta_X2 = -0.69,
                     gamma_X1 = 1,
                     gamma_X2 = 1,
                     gamma_X1X2 = 0) {
  
  flnm <- paste0("resMat", "pi_R1", pi_R1, "beta_X1", beta_X1, "beta_X2", beta_X2, "gamma_X1", gamma_X1, "gamma_X2", gamma_X2, "gamma_X1X2", gamma_X1X2, ".RDS")
  dat <- readRDS(here("Outputs", flnm))
  dat <- as_tibble(dat)
  dat
  
}

pivot_metric <- function(data, metric = "CS") {
  data %>%
    select(ends_with(metric)) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(-iteration,
                 names_to = "method",
                 values_to = metric) %>%
    mutate(Imp_method = ifelse(str_detect(method, "\\d") == TRUE, paste("MI", sep = ": ", parse_number(method)), 
                               ifelse(str_detect(method, "regr") == TRUE, "Regr. Imp.",
                                      ifelse(str_detect(method, "CC") == TRUE, "Complete Case", "No missing"))),
           Y_inc = ifelse(str_detect(method, "no_Y") == TRUE | str_detect(method, "_noY") == TRUE, "No Y in Dev or Val",
                          ifelse(str_detect(method,"YnoY") == TRUE, "Y in Dev NOT val",
                                 ifelse(str_detect(method,"_Y_") == TRUE, "Y in Dev AND Val", "No missing/Complete Case"))),
           Ind = ifelse(str_detect(method, "_ind_") == TRUE, "Missing Indicator", "No Missing Indicator"),
           Imp_method = factor(Imp_method, levels = c("Complete Case", "No missing", "MI: 5", "MI: 10", "MI: 20", "MI: 50", "Regr. Imp.")))
}

## For now focus on gamma_X1 and gamma_X2 = 1 as this is the stronger effect.
MCAR1 <- load_raw(pi_R1 = 0.5,
                  beta_X1 = 0,
                  beta_X2 = 0,
                  gamma_X1 = 1,
                  gamma_X2 = 1,
                  gamma_X1X2 = 0)

MAR1 <- load_raw(pi_R1 = 0.5,
                 beta_X1 = 0,
                 beta_X2 = 1.1,
                 gamma_X1 = 1,
                 gamma_X2 = 1,
                 gamma_X1X2 = 0)

MNAR1 <- load_raw(pi_R1 = 0.5,
                  beta_X1 = 1.1,
                  beta_X2 = 1.1,
                  gamma_X1 = 1,
                  gamma_X2 = 1,
                  gamma_X1X2 = 0)


## might be useful to reshape the data into long format, so we can use groups.
## so e.g. have all the C-statistics in the same column, with each row representing a "method", then we can pass a group argument to ggplot.

pivot_metric <- function(data, metric = "CS") {
  data %>%
    select(ends_with(metric) & !starts_with("CC")) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(-iteration,
                 names_to = "method",
                 values_to = metric) %>%
    mutate(Imp_method = ifelse(str_detect(method, "\\d") == TRUE, paste("MI", sep = ": ", parse_number(method)), 
                               ifelse(str_detect(method, "regr") == TRUE, "Regr. Imp.",
                                      ifelse(str_detect(method, "CC") == TRUE, "Complete Case", "No missing"))),
           Y_inc = ifelse(str_detect(method, "no_Y") == TRUE | str_detect(method, "_noY") == TRUE, "No Y in Dev or Val",
                          ifelse(str_detect(method,"YnoY") == TRUE, "Y in Dev NOT val",
                                 ifelse(str_detect(method,"_Y_") == TRUE, "Y in Dev AND Val", "No missing"))),
           Ind = ifelse(str_detect(method, "_ind_") == TRUE, "Missing Indicator", "No Missing Indicator"),
           Imp_method = factor(Imp_method, levels = c("No missing", "MI: 5", "MI: 10", "MI: 20", "MI: 50", "Regr. Imp.")))
}