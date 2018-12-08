# #############################################
# Question 4
# #############################################

# In this question you will explore Type I error
# via simple (individual observation-level) bootstrapping
# and block bootstrapping

# a) Estimate the Type I error rate with alpha = 0.05
# (i.e. p-value = 0.05) when using simple  bootstrapping.
# One approach to estimate the Type I error rate
# is to construct 95% confidence intervals using
# quantile() on the boostrap estimates and calculating
# the proportion of confidnece intervals that do not
# include the true beta_1, which is 0. You will need to
# create a loop within the main loop I have constructed below,
# i.e. create a nested loop. Computation may take some minutes.

set.seed(100)
n.sims <- 100 
# Number of simulations
n.boot <- 100 
# Bootstrap iterations to 100

beta1 <- c() # Store Beta_1_hat
ci    <- matrix(0,n.sims,2) # Store Confidence Intervals
True  <- matrix(0,n.sims,1) # =1 if outside of CI; =0 if inside of CI

for ( i in seq_len(n.sims)) {
  clustered.df.4 <- gen_cluster( param = c(1, 0), N_g = 100, n_cluster = 50, rho = .5)
  for(k in seq_len( n.boot ) ){
    bsample <- sample( nrow( clustered.df.4 ),size = 5000, replace = TRUE)
    bdta    <- clustered.df.4[bsample,]
    reg4    <- lm(y ~ x, data = bdta)
    beta1[k]   <- c( summary( reg4 )$coefficients[2, 1] )
  }
  low.ci  <- quantile(beta1, probs = .025) # 0.025
  high.ci <- quantile(beta1, probs = .975) # 0.975
  ci[i,]  <- cbind(low.ci,high.ci)
  True[i] <- ifelse(low.ci > 0 | high.ci < 0, 1 , 0)
  cat(i, base::date(), "\n")
}
# Compute the proportion of "1": b_1_hat outside of CI
prop.True <- (sum(True) / n.sims) * 100
cat('Proportion of b_1_hat outside of CI:',prop.True,'%',"\n") # 80%

# b) Do the same as in (a), but with block boostrapping,
# where blocks are the clusters. Note: Constructing
# block bootstrap samples via elegant code may be
# difficult. One approach to accomplish the task is
# with concise syntax is to use split(),
# do.call(), rbind(), and clever list indexing.

set.seed(100)
n.sims    <- 100
n.boot    <- 100
n.blocks  <- 25 # number of blocks to select on the bootstrapping procedure

ci.b        <- matrix(0, n.sims, 2)
beta1.block <- c()
True.block  <- matrix(0, n.sims, 1)

for ( i in seq_len(n.sims)) {
  clustered.df.4 <- gen_cluster(param = c(1, 0), N_g = 100, n_cluster = 50, rho = .5)
  clustered.df.4$cluster.factor <- as.factor(clustered.df.4$cluster)
  for(k in seq_len(n.boot)){
    bsample.block <- subset(clustered.df.4, cluster %in% c(sample(1:50,size=n.blocks)))
    reg4.block    <- lm(y ~ x, data=bsample.block)
    beta1.block[k]   <- c(summary(reg4.block)$coefficients[2, 1])
  }
  low.ci.b      <- quantile(beta1.block, probs = .025) # 0.025
  high.ci.b     <- quantile(beta1.block, probs = .975) # 0.975
  ci.b[i,]      <- cbind(low.ci.b, high.ci.b)
  True.block[i] <- ifelse(low.ci.b > 0 | high.ci.b < 0, 1, 0)
  cat(i, base::date(), "\n")
}
# Compute the proportion of "1": b_1_hat outside of CI
prop.True.block <- (sum(True.block) / n.sims) * 100
cat('Proportion of b_1_hat outside of CI:',prop.True.block,'%',"\n") # 5%
