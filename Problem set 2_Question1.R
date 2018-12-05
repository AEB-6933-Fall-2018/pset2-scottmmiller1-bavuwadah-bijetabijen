# Make sure to install these packages:
install.packages("mvtnorm")
install.packages("lfe")
install.packages("ICC")
library("mvtnorm")
library("lfe")
library("ICC")
library('plm')
library("lmtest")


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



# QUESTION 1

# Using the following data generating process that creates the
# clustered.df.1 dataframe, do the following:

clustered.df.1 <- gen_cluster()

# Data generating process (DGP):

set.seed(100)

clustered.df.1 <- gen_cluster(param = c(1, 0), N_g = 100, n_cluster = 50, rho = .5)
# True DGP is y = 1 + 0 * x + e, so beta_1 = 0
clustered.df.1$cluster.factor <- as.factor(clustered.df.1$cluster)

# a) Report the standard error that results from an OLS regression
# of y on x , i.e. y = beta_0 + beta_1 * x + e,  without any 
# correction for clustering or heteroskedasitity (use lm() )

reg <- lm( y ~ x, data = clustered.df.1 )
summary( reg )



# b) Report the intra-class correlation for clustered.df.1$x that is 
# estimated via ICCest()

Reg.ICC <- ICCest( cluster.factor, clustered.df.1$x, data = clustered.df.1, CI.type = "S" )
Reg.ICC

# c) Report the intra-class correlation of the residual 
# of the regression y = beta_0 + beta_1 * x + e

resid <- reg$residuals

Resid.ICC <- ICCest( cluster.factor, resid, data = clustered.df.1, CI.type = "S" )
Resid.ICC

# d) Report the standard error that results from an OLS regression
# of y on x _with_ correction for clustering. Use felm(). Don't
# use fixed effects.

reg2 <- felm(y ~ x + G(clustered.df.1$cluster.factor), exactDOF = T, data = clustered.df.1)
alpha <- getfe(reg2)
summary(reg2)

