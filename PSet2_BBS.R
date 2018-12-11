rm(list = ls())

# install.packages("mvtnorm")
# install.packages("lfe")
# install.packages("ICC")

library("mvtnorm")
library("lfe")
library("ICC")


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


# 3.

# a)

set.seed(100)

clustered.df.3 <- gen_cluster(param = c(1, 0), N_g = 500, n_cluster = 500, rho = .5)
# True DGP is y = 1 + 0 * x + e, so beta_1 = 0
clustered.df.3$cluster.factor <- as.factor(clustered.df.3$cluster)

reg1 <- lm (y ~ x, data = clustered.df.3)

var1 <- summary(reg1)$coefficients[2,2]^2

# Without correction for clustering, the variance of the estimate of beta_1 is 1.97053e-06.

reg2 <- felm(y ~ x|cluster.factor, data=clustered.df.3)

var2 <- summary(reg2)$coefficients[1,2]^2

# With correction for clustering, the variance of the estimate of beta_1 is 1.996816e-06.


ratio <- var2/var1

# The ratio of these variances (corrected/uncorrected) is 1.01334.


# b)

e <- resid(reg1)

rho_xk <- ICCest (cluster.factor, x, data = clustered.df.3, alpha = 0.05, CI.type = "S")

rho_u <- ICCest (cluster.factor, e, data = clustered.df.3, alpha = 0.05, CI.type = "S")

N_g <- 500

vif <- 1 + rho_xk$ICC*rho_u$ICC*(N_g - 1)

# Here, the variance inflation factor is 135.2376.


# c)

# Higher clustering level:

clustered.df.3$cluster.higher.factor <- as.factor(trunc(clustered.df.3$cluster/5))

reg3 <- lm (y ~ x, data = clustered.df.3)

var3 <- summary(reg1)$coefficients[2,2]^2

# Without correction for clustering, the variance of the estimate of beta_1 is 1.97053e-06.

reg4 <- felm(y ~ x|cluster.higher.factor, data=clustered.df.3)

var4 <- summary(reg4)$coefficients[1,2]^2

# With correction for clustering, the variance of the estimate of beta_1 is 2.0141e-06.


ratio_c <- var4/var3

# The ratio of these variances (corrected/uncorrected) is 1.022111.


# d)

e_d <- resid(reg3)

rho_xk_d <- ICCest (cluster.higher.factor, x, data = clustered.df.3, alpha = 0.05, CI.type = "S")

rho_u_d <- ICCest (cluster.higher.factor, e_d, data = clustered.df.3, alpha = 0.05, CI.type = "S")

N_g <- 500

vif_d <- 1 + rho_xk_d$ICC*rho_u_d$ICC*(N_g - 1)

# Here, the variance inflation factor is 5.816075.




