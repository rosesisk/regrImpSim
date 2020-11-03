#################################################################
## Title: Simulation function(s) for Regression Imputation study
#################################################################
library(tidyverse)
library(mice)
library(broom)
library(pROC)

expit <- function(x) {1/(1+exp(-x))}

genDataSingle <- function(n = 10000,
                          pi_Y = 0.1,
                          pi_R1 = 0.5,
                          sigma_X1 = 1,
                          sigma_X2 = 1,
                          alpha_0 = 0,
                          eta_0 = 0,
                          beta_0 = 0,
                          gamma_0 = 0,
                          beta_X1 = 0,
                          beta_X2 = 0,
                          gamma_X1 = 1,
                          gamma_X2 = 1,
                          gamma_X1X2 = 1)

{
  df <- tibble(X1 = rnorm(n = n, mean = alpha_0, sd = sigma_X1),
               X2 = rnorm(n = n, mean = eta_0, sd = sigma_X2),
               R1 = rbinom(n = n, size = 1, prob = expit(beta_0 + beta_X1 * X1 + beta_X2 * X2)),
               Y = rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_X1 * X1 + gamma_X2 * X2 + gamma_X1X2 * X1*X2))
  )
  
  df
  
}

#df <- genDataSingle() ## create a dataset to test functions with default param values.


#### this function should input a simulated dataset, output the relevant performance metrics in a matrix
runSim <- function(df) {
  
  n <- nrow(df)
  
  df_comp <- df 
  
  df$X1[df$R1 == 1] <- NA
  

  #randomly split data into training and test sets
  rid_test <- sample(1:n, size = n/2, replace = F)
  
  df_test <- df[rid_test, ]
  df_train <- df[-rid_test, ]

  
  ## Complete case model - need to be careful with how this is compared to others since sample size is different.
cccoefs <- glm(Y ~ X1 + X2 + X1:X2, family = binomial(link = "logit"), data = df_train) %>%
  tidy() %>%
  pull(estimate)
names(cccoefs) <- c("beta_0", "beta_X1", "beta_X2", "beta_X1X2")

### now get coefficients for dataset with no missing items. 
df_comp_test <- df_comp[rid_test,]
df_comp_train <- df_comp[-rid_test,]

nomiscoefs <- glm(Y ~ X1 + X2 + X1:X2, family = binomial(link = "logit"), data = df_comp_train) %>%
  tidy() %>%
  pull(estimate)
names(nomiscoefs) <- c("beta_0", "beta_X1", "beta_X2", "beta_X1X2")

### build a function to run multiple imputation with/without Y in imp model. Save out parameters.
### first setup predictorMatrix so that Y is not used to impute X1.
dummyMice <- mice(df_train, m = 1, maxit = 0)
predMat <- dummyMice$predictorMatrix
predMat[1, 4] <- 0 ## do not allow Y in imputation model for X1.

### runs MI under relevant conditions. Outputs mids object.
miceY <- function(df, m) {
  
  ## imputes separately within subsets defined by Y, as per Tilling et al, 2016.
mice_Y0 <- df %>%
  filter(Y == 0) %>%
  mice(m = m, method = "norm", maxit = 1, print = F)

mice_Y1 <- df %>%
  filter(Y == 1) %>%
  mice(m = m, method = "norm", maxit = 1, print = F)

mice_Y <- rbind(mice_Y0, mice_Y1)

}

MIrun <- function(df, m = 10, Y = TRUE) {
  if (Y == TRUE){
    MI <- miceY(df, m)}
  
  else if (Y == FALSE)  {
    MI <- mice(df, m = m, method = "norm", predictorMatrix = predMat,  maxit = 1, print = F) }
  MI }

## recovers relevant parameter estimates from MI object.
MIparams <- function(MI, m, ind = FALSE){  
  if (ind == TRUE){
    
  mod <- with(MI, glm(Y ~ X1 + X2 + X1:X2 + R1, family = binomial(link = "logit")))
  params <- summary(pool(mod))$estimate
  names(params) <- c("beta_0", "beta_X1", "beta_X2", "beta_R1", "beta_X1X2")
  
  }
   
  else if (ind == FALSE) {
  mod <- with(MI, glm(Y ~ X1 + X2 + X1:X2, family = binomial(link = "logit")))
  params <- summary(pool(mod))$estimate
  names(params) <- c("beta_0", "beta_X1", "beta_X2", "beta_X1X2")
  }
  
  params
  
}

## packages up MIrun and MItrain for use in training data - return parameter estimates with one call. 
MItrain <- function(df, m = 10, Y = TRUE, ...) {
  MI <- MIrun(df, m = m, Y = Y)
  MIparams(MI, ...) }

#MItrain(df, ind = F) #testing

#### Now do the same for Regression Imputation, fixed and stochastic. The fixed version uses a deterministic prediction from a fixed model.
## Stochastic adds noise to this.

m <- c(5, 10, 20, 50)

outp_mat_Y_noind <- sapply(m, function(x) MItrain(df_train, m = x, Y = TRUE, ind = F))
outp_mat_Y_ind <- sapply(m, function(x) MItrain(df_train, m = x, Y = TRUE, ind = T))
#store parameters from each m for models including Y, with/without missing indicator
outp_mat_noY_noind <- sapply(m, function(x) MItrain(df_train, m = x, Y = FALSE, ind = F))
outp_mat_noY_ind <- sapply(m, function(x) MItrain(df_train, m = x, Y = FALSE, ind = T))
#as above for models without Y

colnames(outp_mat_Y_noind) <- paste("Imps", m)
colnames(outp_mat_Y_ind) <- paste("Imps", m)
colnames(outp_mat_noY_noind) <- paste("Imps", m)
colnames(outp_mat_noY_ind) <- paste("Imps", m)

  ### Initially look at imputing data separately in the training set, using corresponding number of imputations.
imputed_test_Y <- lapply(m, function(x) MIrun(df_test, m = x, Y = TRUE)) ## produces list of mids objects
imputed_test_noY <- lapply(m, function(x) MIrun(df_test, m = x, Y = FALSE)) 

### Apply each model derived in training set to corresponding m datasets in this list. 
### create a list of stacked datasets, each corresponding to an m, number of imputations. 
### write a function that takes m (imputed datasets), set of Betas. Outputs the predicted Y and linear predictor.


#### Function to take a set of multiply imputed datasets (mids object) and output the predicted Y and linear predictor, 
## based on the appropriate coefficients derived earlier.

test_preds <- function(imps, ind = T, Y = T) {
 
  nimps = imps$m
  index <- which(m == nimps)
  
  if (Y == TRUE & ind == T){
    betas <- outp_mat_Y_ind[,index]
  }
  else if(Y == TRUE & ind == F){
    betas <- outp_mat_Y_noind[,index]
  } 
  else if (Y == FALSE & ind == T) {
    betas <- outp_mat_noY_ind[,index]
  }
  else if (Y == FALSE & ind == F){
    betas <- outp_mat_noY_noind[,index]
  }
    ## append a column of 1s to this. 
  ## and select the relevant columns IN ORDER OF COEFFICIENTS. X1, X2, R1, X1:X2
  if (ind == T){
  des_mat <- complete(imps, action = "long", include = F) %>%
    mutate(intercept = rep(1, nimps * 5000),
           X1X2 = X1*X2) %>%
    select(intercept, X1, X2, R1, X1X2, Y, .imp)}
  
  else if (ind == F){
    des_mat <- complete(imps, action = "long", include = F) %>%
      mutate(intercept = rep(1, nimps * 5000),
             X1X2 = X1*X2) %>%
      select(intercept, X1, X2, X1X2, Y, .imp)}
  
  pred_LP <- as.matrix(des_mat[,1:(length(betas))]) %*% betas
  pred_Y <- expit(pred_LP) 
  
  outp <- tibble(pred_LP = pred_LP, pred_Y = pred_Y, Y = des_mat$Y, .imp = des_mat$.imp)
  outp
}


#### CHECK SOME OF THESE VALUES MANUALLY. Unsure if function is pulling out all correct estimates.
test_predictions_Y_noind <- lapply(1:length(m), 
                                   function(x) test_preds(imputed_test_Y[[x]], 
                                                          ind = F, 
                                                          Y = T))
test_predictions_Y_ind <- lapply(1:length(m), 
                                 function(x) test_preds(imputed_test_Y[[x]], 
                                                        ind = T, 
                                                        Y = T))
test_predictions_noY_ind <- lapply(1:length(m), 
                                   function(x) test_preds(imputed_test_noY[[x]], 
                                                          ind = T, 
                                                          Y = F))
test_predictions_noY_noind <- lapply(1:length(m), 
                                     function(x) test_preds(imputed_test_noY[[x]], 
                                                            ind = F, 
                                                            Y = F)) 
test_predictions_YnoY_noind <- lapply(1:length(m), 
                                      function(x) test_preds(imputed_test_noY[[x]], 
                                                             ind = F, 
                                                             Y = T)) 
test_predictions_YnoY_ind <- lapply(1:length(m), 
                                    function(x) test_preds(imputed_test_noY[[x]], 
                                                           ind = T, 
                                                           Y = T))

## combine into one list, then can apply the performance metrics all at once. Just be sure to retain order/name correctly within list.
test_predictions_mi <- do.call(c, list(test_predictions_Y_noind,
                                    test_predictions_Y_ind,
                                    test_predictions_noY_noind,
                                    test_predictions_noY_ind,
                                    test_predictions_YnoY_noind,
                                    test_predictions_YnoY_ind))


names(test_predictions_mi) <- c(paste(m, "Y_noind", sep = "_"),
                             paste(m, "Y_ind", sep = "_"), 
                             paste(m, "noY_noind", sep = "_"),
                             paste(m, "noY_ind", sep = "_"),
                             paste(m, "YnoY_noind", sep = "_"),
                             paste(m, "YnoY_ind", sep = "_"))

##### Now train models using regression imputation.
frimp_Y1 <- df_train %>%
  filter(Y==1) %>%
  mice(method = "norm.predict", m = 1, maxit = 1, print= F)

frimp_Y0 <- df_train %>%
  filter(Y==0) %>%
  mice(method = "norm.predict", m = 1, maxit = 1, print= F)

fixRegImpY <- rbind(frimp_Y1, frimp_Y0)
fixRegImpnoY <- mice(df_train, method = "norm.predict", m = 1, maxit = 1, predictorMatrix = predMat, print = F)

## extract coefficients
fixRegYnoIndCoefs <- with(fixRegImpY, glm(Y ~ X1 + X2 + X1:X2, family = binomial(link = "logit"))) %>%
  summary %>%
  select(estimate) %>%
  pull

fixRegnoYnoIndCoefs <- with(fixRegImpnoY, glm(Y ~ X1 + X2 + X1:X2, family = binomial(link = "logit"))) %>%
  summary %>%
  select(estimate) %>%
  pull

fixRegYIndCoefs <- with(fixRegImpY, glm(Y ~ X1 + X2 + R1 + X1:X2, family = binomial(link = "logit"))) %>%
  summary %>%
  select(estimate) %>%
  pull

fixRegnoYIndCoefs <- with(fixRegImpnoY, glm(Y ~ X1 + X2 + R1 + X1:X2, family = binomial(link = "logit"))) %>%
  summary %>%
  select(estimate) %>%
  pull

### run imputation in test data using regression imp.
fixRegImpYTestImp <- mice(df_test, method = "norm.predict", m = 1, maxit = 1, print = F)
fixRegImpnoYTestImp <- mice(df_test, method = "norm.predict", m = 1, maxit = 1, predictorMatrix = predMat, print = F)

### Use coefficients derived in training data to append predicted Y and LP to test data.
### Need to do this with/without missing indicator for each combo of YY, noYnoY, YnoY. 
test_preds_regimp <- function(Ydev = T, Ytest = T, ind = T) {
  if (Ydev == T & ind == T) {
    betas <- fixRegYIndCoefs
  } else if (Ydev == F & ind == T){
    betas <- fixRegnoYIndCoefs 
  } else if (Ydev == T & ind == F) {
    betas <- fixRegYnoIndCoefs
  } else if (Ydev == F & ind == F){
    betas <- fixRegnoYnoIndCoefs
  }
  
  if (Ytest == T){
    impd <- fixRegImpYTestImp
  } else if(Ytest == F){
    impd <- fixRegImpnoYTestImp
  }
  
  if (ind == T) {
  test_dat <- complete(impd) %>%
    mutate(intercept = rep(1, 5000),
           X1X2 = X1*X2) %>%
    select(intercept, X1, X2, R1, X1X2, Y)
  }
  else if (ind == F){
    test_dat <- complete(impd) %>%
      mutate(intercept = rep(1, 5000),
             X1X2 = X1*X2) %>%
      select(intercept, X1, X2, X1X2, Y)
  }
  
  pred_LP <- as.matrix(test_dat[,1:length(betas)]) %*% betas
  pred_Y <- expit(pred_LP)

  outp <- tibble(pred_LP = pred_LP, pred_Y = pred_Y, Y = test_dat$Y)
  outp
  
}

test_pred_fixreg_YnoInd <- test_preds_regimp(Ydev = T, Ytest = T, ind = F)
test_pred_fixreg_Yind <- test_preds_regimp(Ydev = T, Ytest = T, ind = T)
test_pred_fixreg_noYnoInd <- test_preds_regimp(Ydev = F, Ytest = F, ind = F)
test_pred_fixreg_noYInd <- test_preds_regimp(Ydev = F, Ytest = F, ind = T)
test_pred_fixreg_YnoYnoInd <- test_preds_regimp(Ydev = T, Ytest = F, ind = F)
test_pred_fixreg_YnoYInd <- test_preds_regimp(Ydev = T, Ytest = F, ind = T)

test_preds_fixreg <- list(test_pred_fixreg_YnoInd,
                                     test_pred_fixreg_Yind,
                                     test_pred_fixreg_noYnoInd,
                                     test_pred_fixreg_noYInd,
                                     test_pred_fixreg_YnoYnoInd,
                                     test_pred_fixreg_YnoYInd)

names(test_preds_fixreg) <- c("regr_Y_noind",
                              "regr_Y_ind",
                              "regr_noY_noind",
                              "regr_no_Y_ind",
                              "regr_YnoY_noind",
                              "regr_YnoY_ind")

################################### Calculate performance measures. ##################################
##### Complete Case Model - need to remember these are not comparable with imputation methods' measures
ccYTest <- df_test %>%
  mutate(pred_LP = cccoefs[1] + cccoefs[2] * X1 + cccoefs[3] * X2 + cccoefs[4]*X1*X2,
        pred_Y = expit(pred_LP))

CSCC <- glm(Y ~ pred_LP, family = "binomial", x = T, y = T, data = ccYTest) %>%
  tidy %>%
  filter(term == "pred_LP") %>%
  pull(estimate)

CITLCC <-  glm(Y ~ offset(pred_LP), family = "binomial", x= T, y = T, data = ccYTest) %>%
  tidy %>%
  pull(estimate)

CStatCC <- auc(Y ~ pred_Y, data = ccYTest) %>%
  as.numeric()

BrierCC <- ccYTest %>%
  summarise(BScore = mean( (pred_Y - Y)^2, na.rm = T)) %>%
  pull

###### no missing data (gold standard)
nomis_test_predictions <- df_comp_test %>%
  mutate(pred_LP = nomiscoefs[1] + nomiscoefs[2] * X1 + nomiscoefs[3] * X2 + nomiscoefs[4] * X1*X2,
         pred_Y = expit(pred_LP))

CS_nomis <- glm(Y ~ pred_LP, family = "binomial", x = T, y = T, data = nomis_test_predictions) %>%
  tidy %>%
  filter(term == "pred_LP") %>%
  pull(estimate)

CITL_nomis <- glm(Y ~ offset(pred_LP), family = "binomial", x= T, y = T, data = nomis_test_predictions) %>%
  tidy %>%
  pull(estimate)

CStat_nomis <- auc(Y ~ pred_Y, data = nomis_test_predictions) %>%
  as.numeric()

BrierStat_nomis <- nomis_test_predictions %>%
  summarise(BScore = mean( (pred_Y - Y)^2)) %>%
  pull
  

########### Regression Imputation 
CS_fixreg <- sapply(test_preds_fixreg, function(x) {
  glm(Y ~ pred_LP, family = "binomial", x= T, y = T, data = x) %>%
    tidy %>%
    filter(term == "pred_LP") %>%
    select(estimate) %>%
    as.numeric
}) 


CITL_fixreg <- sapply(test_preds_fixreg, function(x) { glm(Y ~ offset(pred_LP), family = "binomial", x= T, y = T, data = x) %>%
  tidy %>%
  select(estimate) %>%
  as.numeric })

C_fixreg <- sapply(test_preds_fixreg, function(x) {auc(Y ~ pred_Y, data = x) %>%
    as.numeric()})

 
Brier_fixreg <- sapply(test_preds_fixreg, function(x) {x %>%
    summarise(BScore = mean( (pred_Y - Y)^2) ) %>%
    pull}) 

## combine together.
res_fixreg <- t(rbind(CS_fixreg, CITL_fixreg, C_fixreg, Brier_fixreg))
### remember that for MI it's a little bit trickier - have to do it on each individual imputed dataset.

#input: list of stacked datasets, all combinations of m and Y.
#output: average performance measures for each combination of m and Y/no Y. Averaged across imputed datasets. 

calculateCS <- function(data) { 
  data %>%
  group_by(.imp) %>%
  do(modCS = glm(Y ~ pred_LP, family = "binomial", x= T, y = T, data = .)) %>%
  tidy(., modCS) %>%
  filter(term == 'pred_LP') %>%
  select(estimate) %>%
  ungroup() %>%
  summarise(mean(estimate)) %>%
    as.numeric() ## so that sapply outputs a vector rather than list of datasets
  }

calculateCITL <- function(data) {
  data %>%
    group_by(.imp) %>%
    do(modCITL = glm(Y ~ offset(pred_LP), family = "binomial", x= T, y = T, data = .)) %>%
    tidy(., modCITL) %>%
    select(estimate) %>%
    ungroup() %>%
    summarise(mean(estimate)) %>%
    as.numeric() 
}

calculateC <- function(data) {
  data %>%
    group_by(.imp) %>%
    do(c_stat = auc(Y ~ pred_Y, data = .)) %>%
    select(c_stat) %>%
    ungroup() %>%
    summarise(mean(as.numeric(c_stat))) %>%
    as.numeric()
}


resMat <- matrix(nrow = length(m) *2, ncol = 4, length(m) * 2) ##create a matrix to store out results (may not be necessary later)

## calculate performance measures for multiply imputed datasets.
CS_mi <- sapply(test_predictions_mi, calculateCS) ## this produces a vector, one for each value of m.
CITL_mi <- sapply(test_predictions_mi, calculateCITL)
C_mi <- sapply(test_predictions_mi, calculateC)
BScore_mi <- sapply(test_predictions_mi, function(x) {x %>%
                      summarise(BScore = mean( (pred_Y - Y)^2) ) %>%
                      pull}) # not working- all NAs

### Combine with regression imputation and CC.
### names of Imputed versions are above
resMatImps <- t(rbind(CS_mi, CITL_mi, C_mi, BScore_mi))


resCC <- c(CSCC, CITLCC, CStatCC, BrierCC)
res_nomis <- c(CS_nomis, CITL_nomis, CStat_nomis, BrierStat_nomis)

resMat <- rbind(resMatImps, res_fixreg, resCC, res_nomis)
colnames(resMat) <- c("CS", "CITL", "C", "BScore")

## reshape to a vector
resVec <- as.vector(resMat)
names(resVec) <- paste(rownames(resMat), rep(colnames(resMat), each = nrow(resMat)), sep = "_") ### CHECK THIS.
resVec

}

simulate_nrun <- function(n = 10000,
                          pi_Y = 0.1,
                          pi_R1 = 0.5,
                          sigma_X1 = 1,
                          sigma_X2 = 1,
                          alpha_0 = 0,
                          eta_0 = 0,
                          beta_0 = 0,
                          gamma_0 = 0,
                          beta_X1 = 0,
                          beta_X2 = 0,
                          gamma_X1 = 1,
                          gamma_X2 = 1,
                          gamma_X1X2 = 1,
                          nrun = 2) {
 
   resMat <- matrix(nrow = nrun, ncol = 128)
   
   m <- c(5, 10, 20, 50)
   Yind_ord <- c("Y_noind", 
                 "Y_ind", 
                 "noY_noind", 
                 "no_Y_ind", 
                 "YnoY_noind", 
                 "YnoY_ind")
   
   measures <- c("CS",
                 "CITL",
                 "C",
                 "BScore")
  
   one_meas <- c(  
   paste(rep(m, length(Yind_ord)), rep(Yind_ord, each = length(m)), sep = "_"),
   paste("regr", Yind_ord, sep = "_"),
   "CC",
   "nomis")

   colnames(resMat) <- paste(one_meas, rep(measures, each=length(one_meas)), sep = "_")
     
     
  for (i in 1:nrun){
   
  data <- genDataSingle(n = n, 
                        pi_Y = pi_Y, 
                        pi_R1 = pi_R1, 
                        sigma_X1 = sigma_X1, 
                        sigma_X2 = sigma_X2, 
                        alpha_0 = alpha_0, 
                        eta_0 = eta_0, 
                        beta_0 = beta_0, 
                        gamma_0 = gamma_0, 
                        beta_X1 = beta_X1, 
                        beta_X2 = beta_X2, 
                        gamma_X1 = gamma_X1, 
                        gamma_X2 = gamma_X2, 
                        gamma_X1X2 = gamma_X1X2)
  
  resMat[i,] <- runSim(data)
  
  }
   flname <- paste0("resMat", "pi_R1", pi_R1, "beta_X1", beta_X1, "beta_X2", beta_X2, "gamma_X1", gamma_X1, "gamma_X2", gamma_X2, "gamma_X1X2", gamma_X1X2, ".RDS")
   saveRDS(resMat, file = paste0("Outputs/", flname))
   
}

