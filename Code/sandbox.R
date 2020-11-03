library(tidyverse)
library(mice)
library(broom)
library(pROC)
library(DescTools)


expit <- function(x) {1/(1+exp(-x))}

genDataSingle <- function(n = 10000,
                          pi_Y = 0.1,
                          pi_R1 = 0.5,
                          sigma_X1 = 1,
                          sigma_X2 = 1,
                          alpha_0 = 0,
                          beta_0 = 0,
                          eta_0 = 0,
                          gamma_0 = 0,
                          gamma_X1 = 1,
                          gamma_X2 = 1,
                          gamma_X1X2 = 1,
                          beta_X1 = 0,
                          beta_X2 = 0)

{
  df <- tibble(X1 = rnorm(n = n, mean = alpha_0, sd = sigma_X1),
               X2 = rnorm(n = n, mean = eta_0, sd = sigma_X2),
               R1 = rbinom(n = n, size = 1, prob = expit(beta_0 + beta_X1 * X1 + beta_X2 * X2)),
               Y = rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_X1 * X1 + gamma_X2 * X2 + gamma_X1X2 * X1*X2))
  )
  
  df
  
}

set.seed(1247585)
df <- genDataSingle()

n <- nrow(df)

df_comp <- df 

df$X1[df$R1 == 1] <- NA
df <- df %>% mutate(X1X2 = X1 * X2)

#randomly split data into training and test sets
rid_test <- sample(1:n, size = n/2, replace = F)

df_test <- df[rid_test, ]
df_train <- df[-rid_test, ]

m3 <- mice(df_train, m = 3, method = "norm", maxit = 1, seed = 187)
m5 <- mice(df_train, m = 5, method = "norm", maxit = 1, seed = 187)

### produces the same set of imputed datasets when m is the same, but once m is changed then this is not the case.
### none of the imputed datasets are the same in sets m3 and m5. So doing mice(m = 100), and taking the first 3 datasets is not the same as running mice(m=3).



