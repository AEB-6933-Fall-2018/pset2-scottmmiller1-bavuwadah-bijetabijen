# Problem Set 2: Cluster-robust Standard Errors and Type I Error
# AEB 6933: Advanced Econometrics
# Fall 2018 
# Travis McArthur 
# Due: Tuesday, December 11 11:59pm

# Make sure to install these packages:
# install.packages("mvtnorm")
# install.packages("lfe")
# install.packages("ICC")
library("mvtnorm")
library("lfe")
library("ICC")



# gen_cluster function is adopted from:
# https://yukiyanai.github.io/teaching/rm1/contents/R/clustered-data-analysis.html


gen_cluster <- function(param = c(.1, .5), N_g = 100, n_cluster = 50, rho = .5) {
  # Function to generate clustered data
  # Required package: mvtnorm
  
  # individual level
  Sigma_i <- matrix(c(1, 0, 0, 1 - rho), ncol = 2)
  values_i <- rmvnorm(n = N_g * n_cluster, sigma = Sigma_i)
  
  # cluster level
  cluster_name <- rep(0:(n_cluster-1), each = N_g)
  # Subtle change from 1:n_cluster to 0:(n_cluster-1) to
  # make this work: as.factor(trunc(clustered.df$cluster/5))
  Sigma_cl <- matrix(c(1, 0, 0, rho), ncol = 2)
  values_cl <- rmvnorm(n = n_cluster, sigma = Sigma_cl)
  
  # predictor var consists of individual- and cluster-level components
  x <- values_i[ , 1] + rep(values_cl[ , 1], each = N_g)
  
  # error consists of individual- and cluster-level components
  error <- values_i[ , 2] + rep(values_cl[ , 2], each = N_g)
  
  # data generating process
  y <- param[1] + param[2]*x + error
  
  df <- data.frame(x, y, cluster = cluster_name)
  return(df)
}


# QUESTION 2

# Using a Monte Carlo simulation with 1000 simulated draws from
# the data generating process below, do the following:

# a) Estimate an OLS regression of y on x and obtain the 
# p-values without clustering correction from lm(). Report 
# the proportion of p-values that are less than 0.05.
# Given that the null hypothesis of beta_1 = 0 is true, 
# is the actual type I error with the uncorrected p-values 
# too high, i.e. is there overrejection?

# b) do the same as in (a), but with the correction for
# clustering, with felm()

# Monte Carlo skeleton

set.seed(100)
n.sims <- 1000

# vectors to store p-values
p.vals.ols <- matrix(NA,n.sims,2)
p.vals.cl <- matrix(NA,n.sims,2)

for ( i in seq_len(n.sims)) {
  clustered.df.2 <- gen_cluster(param = c(1, 0), N_g = 100, n_cluster = 50, rho = .5)
  clustered.df.2$cluster.factor <- as.factor(clustered.df.2$cluster)
  
  # Do something here to obtain the p-values
  # Hint: use str() on summary(lm(y ~ x, data = clustered.df.2)) to get the p-values directly
  
  # OLS
  ols <- summary(lm(y ~ x, data = clustered.df.2)) # estimation
  p.vals.ols[i,1] <- ols$coefficients[2,4] # capture p-values
  ifelse( p.vals.ols[i,1] < .05 , p.vals.ols[i,2] <- 1 , p.vals.ols[i,2] <- 0 ) # gen row 2 = 1 if p-val < .05
  
  # Cluster correction
  cl <- summary(felm(y ~ x | clustered.df.2$cluster.factor , data = clustered.df.2)) # estimation
  p.vals.cl[i,1] <- cl$coefficients[1,4] # capture p-values
  ifelse( p.vals.cl[i,1] < .05 , p.vals.cl[i,2] <- 1 , p.vals.cl[i,2] <- 0 ) # gen row 2 = 1 if p-val < .05
  
}

# type I error
mean(p.vals.ols[,2]) # ols = .7
mean(p.vals.cl[,2]) # cluster correction = .05

