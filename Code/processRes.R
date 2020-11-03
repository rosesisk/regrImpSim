#################################################
## Process and combine outputs from simulation
## Date created: 06/07/2020
## Author: Rose Sisk
#################################################
library(tidyverse)
library(purrr)
library(ggpubr)
library(wesanderson)
library(here)
library(pivottabler)

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



pivot_metrics <- function(data) {

C <- data %>%
  select(ends_with("C")) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(-iteration,
               names_to = "method",
               values_to = "C")

CS <- data %>%
  select(ends_with("CS")) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(-iteration,
               names_to = "method",
               values_to = "CS") %>%
  select(-method)

CITL <- data %>%
  select(ends_with("CITL")) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(-iteration,
               names_to = "method",
               values_to = "CITL") %>%
  select(-method)

BScore <- data %>%
  select(ends_with("BScore")) %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(-iteration,
               names_to = "method",
               values_to = "BScore") %>%
  select(-method)

C %>%
  cbind(CS, CITL, BScore) %>%
  select(1, method, C, CS, CITL, BScore) %>%
  mutate(Imp_method = ifelse(str_detect(method, "\\d") == TRUE, paste("MI", sep = ": ", parse_number(method)), 
                             ifelse(str_detect(method, "regr") == TRUE, "Regr. Imp.",
                                    ifelse(str_detect(method, "nomis") == TRUE, "No missing", "Complete case"))),
         Y_inc = ifelse(str_detect(method, "no_Y") == TRUE | str_detect(method, "_noY") == TRUE, "No Y in Dev or Val",
                        ifelse(str_detect(method,"YnoY") == TRUE, "Y in Dev NOT val",
                               ifelse(str_detect(method,"_Y_") == TRUE, "Y in Dev AND Val", "No imputation"))),
         Ind = ifelse(str_detect(method, "_ind_") == TRUE, "Missing Indicator", "No Missing Indicator"),
         Imp_method = factor(Imp_method, levels = c("Complete case","No missing", "MI: 5", "MI: 10", "MI: 20", "MI: 50", "Regr. Imp."))) %>%
  select(-method)

}


MCAR1_summ <- pivot_metrics(MCAR1)
MAR1_summ <- pivot_metrics(MAR1)
MNAR1_summ <- pivot_metrics(MNAR1)

MNAR1_summ %>%
  #filter(Imp_method != "Complete case") %>% ## maybe use this as CC metrics are not comparable - different sample size.
  ggplot(aes(x = Imp_method, y = C, fill = Ind)) +
  geom_boxplot() +
  rotate_x_text() +
  facet_grid(Y_inc~., scales = "free_y") +
  xlab("Imputation method") +
  ylab("C-statistic") +
  scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
  ggtitle("MNAR, 50% missing X1") +
  coord_flip() 


## now merge together and try t
## boxplot for each regression imputation method

MNAR1_C %>%
ggplot(aes(x = Imp_method, y = C, fill = Ind)) +
  geom_boxplot() +
  rotate_x_text() +
  facet_grid(Y_inc~., scales = "free_y") +
  xlab("Imputation method") +
  ylab("C-statistic") +
  scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
  ggtitle("MNAR, 50% missing X1") +
  coord_flip() 

MNAR1_BScore %>%
  ggplot(aes(x = Imp_method, y = BScore, fill = Ind)) +
  geom_boxplot() +
  rotate_x_text()+
  facet_grid(Y_inc~., scales = "free_y") +
  xlab("Imputation method") +
  ylab("Brier Score") +
  scale_fill_manual(values=wes_palette(n=2, name="Royal1"))+
  ggtitle("MNAR, 50% missing X1") +
  coord_flip() 

MCAR1_CITL %>%
  ggplot(aes(x = Imp_method, y = CITL, fill = Ind)) +
  geom_boxplot() +
  rotate_x_text()+
  facet_wrap(~Y_inc, ncol = 4, scales = "free_x") +
  xlab("Imputation method") +
  ylab("Calibration in the large") +
  scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
  ggtitle("MCAR, 50% missing X1")

CS %>%
  ggplot(aes(x = Imp_method, y = CS, fill = Ind)) +
  geom_boxplot() +
  rotate_x_text() +
  facet_wrap(~Y_inc, ncol = 4, scales = "free_x") +
  xlab("Imputation method") +
  ylab("Calibration slope") +
  scale_fill_manual(values=wes_palette(n=2, name="Royal1"))+
  ggtitle("MNAR X1 and X2, 50% missing X1")



#### Create a complex table of results, one for each parameter combination of interest: 
##first take MCAR1/MAR1/MNAR1. See how it looks in 1 table or a set of separate tables.
MCAR1_C %>%
  group_by(Y_inc, Ind, Imp_method) %>%
  summarise(mean_C = mean(C),
            SE_C = empSE(C))

pt <- PivotTable$new()
pt$addData(MCAR1_C)
pt$addColumnDataGroups("Y_inc")
pt$addColumnDataGroups("Ind")
pt$addRowDataGroups("Imp_method")
pt$defineCalculation(calculationName = "Mean", summariseExpression = "mean(C)")
pt$renderPivot()

qhpvt(MCAR1_C, "Imp_method", c("Y_inc", "Ind"), 
      c("Mean C"="mean(C)", "SE C"="empSE(C)"), 
      formats=list("%.3f", "%.4f"))

pvt <- 


######SAVED PLOTS IN FOLDER ARE USING t1 <- load_raw(pi_R1 = 0.5,
#beta_X1 = 0,
#beta_X2 = 1.1,
#gamma_X1 = 1,
#gamma_X2 = 1,
#gamma_X1X2 = 0)
#AND
#t1 <- load_raw(pi_R1 = 0.5,
    #           beta_X1 = 1.1,
     #          beta_X2 = 1.1,
        #       gamma_X1 = 1,
     #          gamma_X2 = 1,
        #       gamma_X1X2 = 0)



### INITIAL NOTES AND OBSERVATIONS
## YnoY method seems to be failing in all scenarios, especially when missingness is at least 50%.
## Seems odd because this is probably  "recommended" - but seemingly inconsistency between dev imp model and val imp
##models is causing issues in performance. 
## NEED to draw attention to pragmatic vs ?? validation as the differences will be huge. 
##Increasing the number of imputed datasets seems to have little to no impact on performance measures.
## clear optimism when using Y in the imputation model (especially for Regr Imp)


## when the predictor effects are very weak, we get HUGE variability in the performance measures. Even when we have MCAR. Not sure why this is.
## Massive variability in the estimates of the performance metrics. This might need looking into.
## Perhaps for now focus on situation where gammas = 1, then revisit/seek advice on the 0.1 scenario at a later date.

res_raw <- pmap_dfr(params, load_res) ## take parameter combinations  and pass these as arguments to load_res()
res <- cbind(params, res_raw)

## might be useful to store median, IQR, range for boxplots too. Or display these separately using boxplots. 


saveRDS(res, "Outputs/combinedResults.RDS")

res2 <- readRDS("Outputs/combinedResults.RDS")

###################################################################
###################################################################
### OLD CODE

pivot_metric <- function(data, metric = "CS") {
  data %>%
    select(ends_with(metric) & !starts_with("CC")) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(-iteration,
                 names_to = "method",
                 values_to = metric) %>%
    mutate(Imp_method = ifelse(str_detect(method, "\\d") == TRUE, paste("MI", sep = ": ", parse_number(method)), 
                               ifelse(str_detect(method, "regr") == TRUE, "Regr. Imp.",
                                      ifelse(str_detect(method, "nomis") == TRUE, "No missing", NA))),
           Y_inc = ifelse(str_detect(method, "no_Y") == TRUE | str_detect(method, "_noY") == TRUE, "No Y in Dev or Val",
                          ifelse(str_detect(method,"YnoY") == TRUE, "Y in Dev NOT val",
                                 ifelse(str_detect(method,"_Y_") == TRUE, "Y in Dev AND Val", "No missing"))),
           Ind = ifelse(str_detect(method, "_ind_") == TRUE, "Missing Indicator", "No Missing Indicator"),
           Imp_method = factor(Imp_method, levels = c("No missing", "MI: 5", "MI: 10", "MI: 20", "MI: 50", "Regr. Imp.")))
}


MAR1_CS <- pivot_metric(MAR1, "CS")
MAR1_CITL <- pivot_metric(MAR1, "CITL")
MAR1_BScore <- pivot_metric(MAR1, "BScore")
MAR1_C <- pivot_metric(MAR1, "C")

MNAR1_CS <- pivot_metric(MNAR1, "CS")
MNAR1_CITL <- pivot_metric(MNAR1, "CITL")
MNAR1_BScore <- pivot_metric(MNAR1, "BScore")
MNAR1_C <- pivot_metric(MNAR1, "C")
