##lets bring some data in
library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)


#################################
#WP Histo 1-45 Samples
#################################

WPH_measures <- read_csv("data/WP_sonora.csv")

#step is VERY important or you get a bunch of NAs

WPH_measures$value <- as.numeric(WPH_measures$value)
WPH_measures$measurement <- as.factor(WPH_measures$measurement)


##get one value each
df1Tall <- WPH_measures %>% group_by(ID, species, worded_health, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

#new df
df4 <- df1Tall %>% pivot_wider(id_cols = c(ID, species, worded_health), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$ID <- paste(df4$species, df4$ID, df4$worded_health, sep = "_")

View(df4)


# Group condition
##HD are included as DD because of previous research (they behave like DD)
## lets also remove any possible NAs before we start 
df4$state <- ifelse(df4$worded_health %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)
df5 <- df4[,c(1,2,3,4, 11, 5, 22, 23, 26, 36)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so they're centered on zero and won't have huge outliers
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##complete.cases??


# Exopercent = x1
df6$x1 <- (df6$mean_prop_exo - mean(df6$mean_prop_exo))/sd(df6$mean_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$mean_gastro - mean(df6$mean_gastro))/sd(df6$mean_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$max_avg_symb - mean(df6$max_avg_symb))/sd(df6$max_avg_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


#values from -1 to 1, 0 means no linear relationship
cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #0.60
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #0.18
cor(df6$mean_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.03
cor(df6$mean_gastro, df6$max_symb_vac, use='complete.obs') #-0.25
cor(df6$mean_gastro, df6$mean_degr_symb, use='complete.obs') #0.78
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #-0.46


##nothing is correlated!!!

# Enumerate species for indexing as random effect
df6$sp <- as.numeric(as.factor(df6$species))

# Number of species
J <- length(unique(df6$species))

# Number of parameters
##this is one more than params because of intercept
K <- 6

# Create identity matrix for Wishart dist'n
# (Number of parameters to estimate (K))
W <- diag(K)


# Load data
data <- list(y = df6$state, 
             group = df6$sp,
             n = dim(df6)[1], 
             J = J,
             x1 = df6$x1,
             x2 = df6$x2,
             x3 = df6$x3, 
             x4 = df6$x4,
             x5 = df6$x5,
             K = K, 
             W = W)


# Define the model in the BUGS language and write a text file
##
sink("model.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1 + exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i] + beta4[group[i]] * x4[i]	+ beta5[group[i]] * x5[i]# linear predictor    
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2] # Exo
    beta2[j] <- BB[j,3] # Gastro
    beta3[j] <- BB[j,4] # Zoox
    beta4[j] <- BB[j,5] # symb size
    beta5[j] <- BB[j,6] # degraded symb
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivariate normal
    
    BB.hat[j,1] <- mu.alpha
	  BB.hat[j,2] <- mu.beta1
	  BB.hat[j,3] <- mu.beta2
	  BB.hat[j,4] <- mu.beta3
	  BB.hat[j,5] <- mu.beta4
	  BB.hat[j,6] <- mu.beta5
	  # Use below for 2nd-level of model
    #BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    #BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    mu.beta4 ~ dnorm(0, 0.0001)
    mu.beta5 ~ dnorm(0, 0.0001)
    #gamma0b ~ dnorm(0, 0.0001)
    #gamma0b2 ~ dnorm(0, 0.0001)
    #gamma0b3 ~ dnorm(0, 0.0001)

    #gamma1b ~ dnorm(0, 0.0001)
    #gamma1b2 ~ dnorm(0, 0.0001)
    #gamma1b3 ~ dnorm(0, 0.0001)
 
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()


# Initial values
inits <- function(){
  list(mu.alpha = rnorm(1), 
       mu.beta1 = rnorm(1), 
       mu.beta2 = rnorm(1), 
       mu.beta3 = rnorm(1), 
       mu.beta4 = rnorm(1), 
       mu.beta5 = rnorm(1),
       BB = matrix(rnorm(J*K),nrow=J,ncol=K),
       Tau.B = rwish(K+1,diag(K))
       #gamma0b = rnorm(1), 
       #gamma1b=rnorm(1),
       #gamma0b2=rnorm(1),
       #gamma0b3=rnorm(1),
       #gamma1b2=rnorm(1),
       #gamma1b3=rnorm(1) 
  )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta1","mu.beta2","mu.beta3", "mu.beta4",
                "mu.beta5", "BB", "Sigma.B"
                #"gamma0b","gamma1b","gamma0b2",
                #"gamma0b3","gamma1b2","gamma1b3"
)


# MCMC settings
ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

data$y <- as.numeric(data$y)

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 300 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.15 minutes at time 2026-01-13 13:34:30.977053.
# 
# mean      sd     2.5%     50%    97.5% overlap0     f  Rhat n.eff
# mu.alpha       -7.131   5.275  -19.444  -6.208    0.542     TRUE 0.967 1.550     8
# mu.beta1        7.165   5.151    0.070   6.003   18.451    FALSE 0.977 1.791     6
# mu.beta2       -7.766   5.161  -19.855  -6.674   -0.759    FALSE 0.987 1.717     6
# mu.beta3       -8.108   8.018  -27.044  -6.319    3.187     TRUE 0.936 1.388    10
# mu.beta4       -5.384   9.404  -26.301  -3.978   10.302     TRUE 0.800 1.165    33
# mu.beta5        0.519   5.820  -14.179   0.973    9.528     TRUE 0.639 1.164   257
# BB[1,1]        -9.483   5.250  -19.926  -8.829   -1.427    FALSE 1.000 2.411     4
# BB[2,1]        -4.820   4.640  -18.109  -4.162    2.971     TRUE 0.907 1.319    15
# BB[1,2]         7.533   5.245   -0.099   6.780   19.407     TRUE 0.973 1.759     6
# BB[2,2]         6.789   5.068    0.452   5.192   18.593    FALSE 0.988 1.839     6
# BB[1,3]        -7.189   4.985  -20.999  -5.929   -0.910    FALSE 0.996 1.955     5
# BB[2,3]        -8.344   5.515  -22.242  -6.907   -0.976    FALSE 0.991 1.529     8
# BB[1,4]        -2.347   3.183   -8.512  -2.271    4.325     TRUE 0.783 1.164    20
# BB[2,4]       -13.811   8.692  -32.581 -11.503   -2.358    FALSE 0.998 2.232     5
# BB[1,5]       -13.846   8.447  -31.616 -13.264   -0.342    FALSE 0.983 2.251     5
# BB[2,5]         2.972   3.058   -1.035   2.343   10.372     TRUE 0.898 1.509     8
# BB[1,6]         2.283   2.649   -0.912   1.656   10.797     TRUE 0.869 1.384    11
# BB[2,6]        -1.234   8.217  -27.159   0.077   11.697     TRUE 0.492 1.225    35
# Sigma.B[1,1]   20.318 122.744    0.198   3.947  128.878    FALSE 1.000 1.255   142
# Sigma.B[2,1]   -1.162  42.914  -31.613  -0.128   23.414     TRUE 0.549 1.124 10446
# Sigma.B[3,1]   -0.490  92.660  -42.299  -0.045   42.562     TRUE 0.517 1.263  3539
# Sigma.B[4,1]  -28.297 240.238 -183.726  -6.553    8.868     TRUE 0.825 1.269   225
# Sigma.B[5,1]   44.578 273.481   -9.432  11.294  286.596     TRUE 0.845 1.239   139
# Sigma.B[6,1]  -17.648 135.480 -146.657  -0.740   16.547     TRUE 0.650 1.271   105
# Sigma.B[1,2]   -1.162  42.914  -31.613  -0.128   23.414     TRUE 0.549 1.124 10446
# Sigma.B[2,2]    7.420  33.989    0.174   1.828   45.726    FALSE 1.000 1.001  1530
# Sigma.B[3,2]   -0.413  35.792  -24.303  -0.008   21.491     TRUE 0.504 1.203  2312
# Sigma.B[4,2]    0.596 103.800  -63.585   0.471   64.947     TRUE 0.574 1.164  1668
# Sigma.B[5,2]   -7.798 128.076 -121.250  -0.898   70.185     TRUE 0.589 1.130  1874
# Sigma.B[6,2]   -1.578  47.106  -37.491  -0.069   28.548     TRUE 0.529 1.016  4213
# Sigma.B[1,3]   -0.490  92.660  -42.299  -0.045   42.562     TRUE 0.517 1.263  3539
# Sigma.B[2,3]   -0.413  35.792  -24.303  -0.008   21.491     TRUE 0.504 1.203  2312
# Sigma.B[3,3]   14.307  96.848    0.186   2.614   94.418    FALSE 1.000 1.204   258
# Sigma.B[4,3]   10.970 220.102  -70.312   0.491  139.896     TRUE 0.573 1.259  3678
# Sigma.B[5,3]  -11.061 245.230 -178.322  -0.801  112.541     TRUE 0.577 1.252  6786
# Sigma.B[6,3]   -0.663  64.497  -53.700  -0.067   48.966     TRUE 0.526 1.196  9733
# Sigma.B[1,4]  -28.297 240.238 -183.726  -6.553    8.868     TRUE 0.825 1.269   225
# Sigma.B[2,4]    0.596 103.800  -63.585   0.471   64.947     TRUE 0.574 1.164  1668
# Sigma.B[3,4]   10.970 220.102  -70.312   0.491  139.896     TRUE 0.573 1.259  3678
# Sigma.B[4,4]   90.519 613.881    0.308  22.036  534.035    FALSE 1.000 1.263   161
# Sigma.B[5,4] -117.377 677.134 -688.552 -32.205    0.326     TRUE 0.953 1.262   107
# Sigma.B[6,4]   29.535 284.313  -54.959   1.499  270.533     TRUE 0.626 1.257   167
# Sigma.B[1,5]   44.578 273.481   -9.432  11.294  286.596     TRUE 0.845 1.239   139
# Sigma.B[2,5]   -7.798 128.076 -121.250  -0.898   70.185     TRUE 0.589 1.130  1874
# Sigma.B[3,5]  -11.061 245.230 -178.322  -0.801  112.541     TRUE 0.577 1.252  6786
# Sigma.B[4,5] -117.377 677.134 -688.552 -32.205    0.326     TRUE 0.953 1.262   107
# Sigma.B[5,5]  175.827 855.077    0.443  53.218 1027.016    FALSE 1.000 1.236    82
# Sigma.B[6,5]  -34.204 278.937 -337.925  -2.361   90.945     TRUE 0.629 1.230   109
# Sigma.B[1,6]  -17.648 135.480 -146.657  -0.740   16.547     TRUE 0.650 1.271   105
# Sigma.B[2,6]   -1.578  47.106  -37.491  -0.069   28.548     TRUE 0.529 1.016  4213
# Sigma.B[3,6]   -0.663  64.497  -53.700  -0.067   48.966     TRUE 0.526 1.196  9733
# Sigma.B[4,6]   29.535 284.313  -54.959   1.499  270.533     TRUE 0.626 1.257   167
# Sigma.B[5,6]  -34.204 278.937 -337.925  -2.361   90.945     TRUE 0.629 1.230   109
# Sigma.B[6,6]   42.707 371.445    0.187   3.713  320.248    FALSE 1.000 1.286   155
# deviance       12.427   7.929    2.329  10.354   31.863    FALSE 1.000 1.388     9
# 
# **WARNING** Rhat values indicate convergence failure. 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 23.4 and DIC = 35.868 
# DIC is an estimate of expected predictive error (lower is better).


##run it to 80% again 

##by species here we go

quantile(out1$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.

#OFRA
quantile(out1$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

#PAST
quantile(out1$sims.list$BB[ ,2,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,6], probs = c(0.1, 0.9))




##make a logistic regression figure for all
View(df6$species)
table(df6$species, df6$worded_health)
df7<-na.omit(df6)


f1 <- ggplot(df6, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f1

f2 <- ggplot(df6, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f2

f3 <- ggplot(df6, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f3

f4 <- ggplot(df6, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f4

f5 <- ggplot(df6, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f5

f1/f2/f3/f4/f5


#################################
#Serial Fall 2021 SCTLD Samples
#################################

SCTLD_measures <- read_csv("data/Serial_measures_MIAMI.csv")

#step is VERY important or you get a bunch of NAs

SCTLD_measures$value <- as.numeric(SCTLD_measures$value)
SCTLD_measures$measurement <- as.factor(SCTLD_measures$measurement)


##get one value each
df1Tall <- SCTLD_measures %>% group_by(Individual, species, worded_health, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

#new df
df4 <- df1Tall %>% pivot_wider(id_cols = c(Individual, species, worded_health), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$ID <- paste(df4$species, df4$Individual, df4$worded_health, sep = "_")

View(df4)


# Group condition
##HD are included as DD because of previous research (they behave like DD)
## lets also remove any possible NAs before we start 
df4$state <- ifelse(df4$worded_health %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)
df5 <- df4[,c(1,2,3,4, 11, 5, 22, 23, 26, 36, 37)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so they're centered on zero and won't have huge outliers
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##complete.cases??


# Exopercent = x1
df6$x1 <- (df6$mean_prop_exo - mean(df6$mean_prop_exo))/sd(df6$mean_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$mean_gastro - mean(df6$mean_gastro))/sd(df6$mean_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$max_avg_symb - mean(df6$max_avg_symb))/sd(df6$max_avg_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #-0.08
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #0.30
cor(df6$mean_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.15
cor(df6$mean_gastro, df6$max_symb_vac, use='complete.obs') #0.08
cor(df6$mean_gastro, df6$mean_degr_symb, use='complete.obs') #0.29
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #0.02


##nothing is correlated!!!

# Enumerate species for indexing as random effect
df6$sp <- as.numeric(as.factor(df6$species))

# Number of species
J <- length(unique(df6$species))

# Number of parameters
##this is one more than params because of intercept
K <- 6

# Create identity matrix for Wishart dist'n
# (Number of parameters to estimate (K))
W <- diag(K)


# Load data
data2 <- list(y = df6$state, 
             group = df6$sp,
             n = dim(df6)[1], 
             J = J,
             x1 = df6$x1,
             x2 = df6$x2,
             x3 = df6$x3, 
             x4 = df6$x4,
             x5 = df6$x5,
             K = K, 
             W = W)


# Define the model in the BUGS language and write a text file
##
sink("model.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1 + exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i] + beta4[group[i]] * x4[i]	+ beta5[group[i]] * x5[i]# linear predictor    
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2] # Exo
    beta2[j] <- BB[j,3] # Gastro
    beta3[j] <- BB[j,4] # Zoox
    beta4[j] <- BB[j,5] # symb size
    beta5[j] <- BB[j,6] # degraded symb
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivariate normal
    
    BB.hat[j,1] <- mu.alpha
	  BB.hat[j,2] <- mu.beta1
	  BB.hat[j,3] <- mu.beta2
	  BB.hat[j,4] <- mu.beta3
	  BB.hat[j,5] <- mu.beta4
	  BB.hat[j,6] <- mu.beta5
	  # Use below for 2nd-level of model
    #BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    #BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    mu.beta4 ~ dnorm(0, 0.0001)
    mu.beta5 ~ dnorm(0, 0.0001)
    #gamma0b ~ dnorm(0, 0.0001)
    #gamma0b2 ~ dnorm(0, 0.0001)
    #gamma0b3 ~ dnorm(0, 0.0001)

    #gamma1b ~ dnorm(0, 0.0001)
    #gamma1b2 ~ dnorm(0, 0.0001)
    #gamma1b3 ~ dnorm(0, 0.0001)
 
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()


# Initial values
inits <- function(){
  list(mu.alpha = rnorm(1), 
       mu.beta1 = rnorm(1), 
       mu.beta2 = rnorm(1), 
       mu.beta3 = rnorm(1), 
       mu.beta4 = rnorm(1), 
       mu.beta5 = rnorm(1),
       BB = matrix(rnorm(J*K),nrow=J,ncol=K),
       Tau.B = rwish(K+1,diag(K))
       #gamma0b = rnorm(1), 
       #gamma1b=rnorm(1),
       #gamma0b2=rnorm(1),
       #gamma0b3=rnorm(1),
       #gamma1b2=rnorm(1),
       #gamma1b3=rnorm(1) 
  )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta1","mu.beta2","mu.beta3", "mu.beta4",
                "mu.beta5", "BB", "Sigma.B"
                #"gamma0b","gamma1b","gamma0b2",
                #"gamma0b3","gamma1b2","gamma1b3"
)


# MCMC settings
ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

data2$y <- as.numeric(data2$y)

out2 <- jags(data2, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out2)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 400 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.093 minutes at time 2026-01-13 11:26:21.711532.
# 
# mean      sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mu.alpha     -2.553   2.452 -6.679 -2.410  1.299     TRUE 0.919 1.203    15
# mu.beta1      0.011   1.902 -3.139  0.016  3.120     TRUE 0.507 1.005 52500
# mu.beta2     -0.237   1.958 -3.305 -0.231  2.880     TRUE 0.596 1.003 48998
# mu.beta3      0.389   2.039 -2.747  0.388  3.467     TRUE 0.653 1.005  7894
# mu.beta4      0.747   2.007 -2.330  0.735  3.852     TRUE 0.767 1.009 52500
# mu.beta5     21.057  10.208  2.298 19.242 41.316    FALSE 0.995 1.789     6
# BB[1,1]      -2.532   1.462 -5.483 -2.309  0.092     TRUE 0.970 1.757     6
# BB[1,2]       0.022   0.350 -0.667  0.021  0.707     TRUE 0.527 1.001  2787
# BB[1,3]      -0.231   0.423 -1.143 -0.209  0.533     TRUE 0.696 1.001  2367
# BB[1,4]       0.398   0.433 -0.392  0.377  1.306     TRUE 0.817 1.004   574
# BB[1,5]       0.738   0.415 -0.008  0.713  1.610     TRUE 0.972 1.001  1687
# BB[1,6]      21.069  10.026  2.374 19.163 41.016    FALSE 0.999 1.827     6
# Sigma.B[1,1]  4.192  52.740  0.136  0.721 18.991    FALSE 1.000 1.029 25877
# Sigma.B[2,1] -0.236  39.601 -6.130 -0.002  5.840     TRUE 0.502 1.156 52500
# Sigma.B[3,1]  0.296  47.547 -5.611  0.000  6.001     TRUE 0.500 1.100 17909
# Sigma.B[4,1]  0.785  97.845 -5.798  0.004  6.121     TRUE 0.504 1.122 42305
# Sigma.B[5,1]  0.038  34.174 -5.904 -0.002  6.113     TRUE 0.497 1.039 46997
# Sigma.B[6,1]  0.194  42.578 -6.064  0.000  5.778     TRUE 0.500 1.067 19267
# Sigma.B[1,2] -0.236  39.601 -6.130 -0.002  5.840     TRUE 0.502 1.156 52500
# Sigma.B[2,2]  4.099  60.750  0.134  0.717 19.012    FALSE 1.000 1.188 10166
# Sigma.B[3,2] -0.095  50.797 -5.979  0.001  5.917     TRUE 0.499 1.210 52500
# Sigma.B[4,2]  0.126  79.544 -5.887  0.000  5.951     TRUE 0.500 1.099 52500
# Sigma.B[5,2] -0.051  43.609 -5.828  0.006  5.966     TRUE 0.493 1.163 52500
# Sigma.B[6,2]  0.047  39.197 -5.886  0.000  5.941     TRUE 0.501 1.193 52500
# Sigma.B[1,3]  0.296  47.547 -5.611  0.000  6.001     TRUE 0.500 1.100 17909
# Sigma.B[2,3] -0.095  50.797 -5.979  0.001  5.917     TRUE 0.499 1.210 52500
# Sigma.B[3,3]  4.291  69.636  0.136  0.714 18.847    FALSE 1.000 1.108 16208
# Sigma.B[4,3]  0.701 109.170 -5.891  0.004  5.793     TRUE 0.505 1.110 52500
# Sigma.B[5,3]  0.217  44.634 -6.022 -0.003  5.541     TRUE 0.497 1.141 52500
# Sigma.B[6,3] -0.192  41.230 -6.196 -0.001  5.515     TRUE 0.501 1.166 52500
# Sigma.B[1,4]  0.785  97.845 -5.798  0.004  6.121     TRUE 0.504 1.122 42305
# Sigma.B[2,4]  0.126  79.544 -5.887  0.000  5.951     TRUE 0.500 1.099 52500
# Sigma.B[3,4]  0.701 109.170 -5.891  0.004  5.793     TRUE 0.505 1.110 52500
# Sigma.B[4,4]  6.339 333.951  0.138  0.732 19.865    FALSE 1.000 1.235 44403
# Sigma.B[5,4]  0.324  69.243 -5.587  0.004  6.287     TRUE 0.504 1.144 47276
# Sigma.B[6,4]  0.189  80.089 -6.086 -0.004  5.721     TRUE 0.496 1.101 22453
# Sigma.B[1,5]  0.038  34.174 -5.904 -0.002  6.113     TRUE 0.497 1.039 46997
# Sigma.B[2,5] -0.051  43.609 -5.828  0.006  5.966     TRUE 0.493 1.163 52500
# Sigma.B[3,5]  0.217  44.634 -6.022 -0.003  5.541     TRUE 0.497 1.141 52500
# Sigma.B[4,5]  0.324  69.243 -5.587  0.004  6.287     TRUE 0.504 1.144 47276
# Sigma.B[5,5]  4.556  68.074  0.136  0.722 18.940    FALSE 1.000 1.087 20988
# Sigma.B[6,5] -0.087  34.462 -5.819  0.001  5.916     TRUE 0.498 1.005 52500
# Sigma.B[1,6]  0.194  42.578 -6.064  0.000  5.778     TRUE 0.500 1.067 19267
# Sigma.B[2,6]  0.047  39.197 -5.886  0.000  5.941     TRUE 0.501 1.193 52500
# Sigma.B[3,6] -0.192  41.230 -6.196 -0.001  5.515     TRUE 0.501 1.166 52500
# Sigma.B[4,6]  0.189  80.089 -6.086 -0.004  5.721     TRUE 0.496 1.101 22453
# Sigma.B[5,6] -0.087  34.462 -5.819  0.001  5.916     TRUE 0.498 1.005 52500
# Sigma.B[6,6]  4.115  73.643  0.136  0.715 19.393    FALSE 1.000 1.138 36604
# deviance     65.674   3.382 61.128 65.026 74.063    FALSE 1.000 1.001  3079
# 
# **WARNING** Rhat values indicate convergence failure. 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 5.7 and DIC = 71.389 
# DIC is an estimate of expected predictive error (lower is better).
# ##run it to 80% again 

##by species here we go

quantile(out2$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.

#ONLY OFRA
quantile(out2$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out2$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out2$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out2$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out2$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))




##make a logistic regression figure for all
View(df6$species)
table(df6$species, df6$worded_health)
df7<-na.omit(df6)


f1 <- ggplot(df6, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("OFRA")), ncol = 9)

f1

f2 <- ggplot(df6, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("OFRA")), ncol = 9)
f2

f3 <- ggplot(df6, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("OFRA")), ncol = 9)

f3

f4 <- ggplot(df6, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("OFRA")), ncol = 9)

f4

f5 <- ggplot(df6, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("OFRA")), ncol = 9)
f5

f1/f2/f3/f4/f5



################################
#Fall 2022 Unknown (Correa/FGB)
################################


f22_measures <- read_csv("data/F22_Unk.csv")

#step is VERY important or you get a bunch of NAs

f22_measures$value <- as.numeric(f22_measures$value)
f22_measures$measurement <- as.factor(f22_measures$measurement)


##get one value each
df1Tall <- f22_measures %>% group_by(individual, species, condition_h, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

#new df
df4 <- df1Tall %>% pivot_wider(id_cols = c(individual, species, condition_h), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$individual <- paste(df4$species, df4$individual, df4$condition_h, sep = "_")

View(df4)


# Group condition
##HD are included as DD because of previous research (they behave like DD)
## lets also remove any possible NAs before we start 
df4$state <- ifelse(df4$condition_h %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)
df5 <- df4[,c(1,2,3,4, 11, 5, 22, 23, 26, 36)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so they're centered on zero and won't have huge outliers
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##complete.cases??


# Exopercent = x1
df6$x1 <- (df6$mean_prop_exo - mean(df6$mean_prop_exo))/sd(df6$mean_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$mean_gastro - mean(df6$mean_gastro))/sd(df6$mean_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$max_avg_symb - mean(df6$max_avg_symb))/sd(df6$max_avg_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #0.16
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #0.13
cor(df6$mean_prop_exo, df6$mean_degr_symb, use='complete.obs') #-0.03
cor(df6$mean_gastro, df6$max_symb_vac, use='complete.obs') #-0.25
cor(df6$mean_gastro, df6$mean_degr_symb, use='complete.obs') #0.17
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #-0.03


##nothing is correlated!!!

# Enumerate species for indexing as random effect
df6$sp <- as.numeric(as.factor(df6$species))

# Number of species
J <- length(unique(df6$species))

# Number of parameters
##this is one more than params because of intercept
K <- 6

# Create identity matrix for Wishart dist'n
# (Number of parameters to estimate (K))
W <- diag(K)


# Load data
data3 <- list(y = df6$state, 
             group = df6$sp,
             n = dim(df6)[1], 
             J = J,
             x1 = df6$x1,
             x2 = df6$x2,
             x3 = df6$x3, 
             x4 = df6$x4,
             x5 = df6$x5,
             K = K, 
             W = W)


# Define the model in the BUGS language and write a text file
##
sink("model.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1 + exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i] + beta4[group[i]] * x4[i]	+ beta5[group[i]] * x5[i]# linear predictor    
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2] # Exo
    beta2[j] <- BB[j,3] # Gastro
    beta3[j] <- BB[j,4] # Zoox
    beta4[j] <- BB[j,5] # symb size
    beta5[j] <- BB[j,6] # degraded symb
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivariate normal
    
    BB.hat[j,1] <- mu.alpha
	  BB.hat[j,2] <- mu.beta1
	  BB.hat[j,3] <- mu.beta2
	  BB.hat[j,4] <- mu.beta3
	  BB.hat[j,5] <- mu.beta4
	  BB.hat[j,6] <- mu.beta5
	  # Use below for 2nd-level of model
    #BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    #BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    mu.beta4 ~ dnorm(0, 0.0001)
    mu.beta5 ~ dnorm(0, 0.0001)
    #gamma0b ~ dnorm(0, 0.0001)
    #gamma0b2 ~ dnorm(0, 0.0001)
    #gamma0b3 ~ dnorm(0, 0.0001)

    #gamma1b ~ dnorm(0, 0.0001)
    #gamma1b2 ~ dnorm(0, 0.0001)
    #gamma1b3 ~ dnorm(0, 0.0001)
 
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()


# Initial values
inits <- function(){
  list(mu.alpha = rnorm(1), 
       mu.beta1 = rnorm(1), 
       mu.beta2 = rnorm(1), 
       mu.beta3 = rnorm(1), 
       mu.beta4 = rnorm(1), 
       mu.beta5 = rnorm(1),
       BB = matrix(rnorm(J*K),nrow=J,ncol=K),
       Tau.B = rwish(K+1,diag(K))
       #gamma0b = rnorm(1), 
       #gamma1b=rnorm(1),
       #gamma0b2=rnorm(1),
       #gamma0b3=rnorm(1),
       #gamma1b2=rnorm(1),
       #gamma1b3=rnorm(1) 
  )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta1","mu.beta2","mu.beta3", "mu.beta4",
                "mu.beta5", "BB", "Sigma.B"
                #"gamma0b","gamma1b","gamma0b2",
                #"gamma0b3","gamma1b2","gamma1b3"
)


# MCMC settings
ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

data3$y <- as.numeric(data3$y)

out3 <- jags(data3, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out3)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 300 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.254 minutes at time 2026-01-13 12:34:04.025738.
# 
# mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# mu.alpha      -0.894 0.539  -1.948  -0.888   0.133     TRUE 0.961 1.000  4448
# mu.beta1      -0.188 0.540  -1.251  -0.184   0.857     TRUE 0.653 1.000  8421
# mu.beta2      -0.051 0.522  -1.063  -0.053   0.976     TRUE 0.548 1.000  4838
# mu.beta3       0.079 0.544  -0.975   0.078   1.142     TRUE 0.567 1.000 24454
# mu.beta4       0.390 0.528  -0.627   0.387   1.425     TRUE 0.802 1.000 52500
# mu.beta5      -0.952 0.613  -2.201  -0.929   0.163     TRUE 0.958 1.002  1101
# BB[1,1]       -0.986 0.433  -1.889  -0.968  -0.183    FALSE 0.993 1.002  1322
# BB[2,1]       -0.990 0.412  -1.836  -0.983  -0.219    FALSE 0.994 1.002  1424
# BB[3,1]       -0.702 0.337  -1.381  -0.699  -0.048    FALSE 0.982 1.000 12929
# BB[1,2]       -0.354 0.445  -1.251  -0.346   0.492     TRUE 0.787 1.001  4803
# BB[2,2]       -0.169 0.437  -1.046  -0.157   0.676     TRUE 0.649 1.001  2725
# BB[3,2]       -0.050 0.302  -0.664  -0.044   0.529     TRUE 0.557 1.001  2680
# BB[1,3]        0.051 0.439  -0.786   0.038   0.949     TRUE 0.537 1.002  1179
# BB[2,3]       -0.141 0.372  -0.854  -0.145   0.605     TRUE 0.651 1.001  4724
# BB[3,3]       -0.071 0.331  -0.713  -0.071   0.590     TRUE 0.590 1.001  1746
# BB[1,4]       -0.130 0.424  -0.992  -0.123   0.687     TRUE 0.615 1.000  6814
# BB[2,4]        0.379 0.387  -0.337   0.363   1.209     TRUE 0.842 1.000 11134
# BB[3,4]       -0.013 0.335  -0.688  -0.005   0.631     TRUE 0.505 1.001  8658
# BB[1,5]        0.530 0.364  -0.150   0.526   1.264     TRUE 0.931 1.001  7395
# BB[2,5]        0.334 0.408  -0.459   0.326   1.154     TRUE 0.798 1.000 52500
# BB[3,5]        0.309 0.379  -0.428   0.304   1.064     TRUE 0.794 1.001  3288
# BB[1,6]       -1.096 0.639  -2.484  -1.039   0.029     TRUE 0.973 1.002  1091
# BB[2,6]       -0.925 0.443  -1.877  -0.894  -0.134    FALSE 0.990 1.005   465
# BB[3,6]       -0.834 0.409  -1.757  -0.801  -0.140    FALSE 0.993 1.003   784
# Sigma.B[1,1]   0.654 1.351   0.107   0.379   2.762    FALSE 1.000 1.017 30411
# Sigma.B[2,1]   0.039 0.748  -0.882   0.017   1.084     TRUE 0.543 1.001 52500
# Sigma.B[3,1]  -0.008 0.753  -0.974  -0.003   0.938     TRUE 0.507 1.007 52500
# Sigma.B[4,1]  -0.009 0.783  -1.017  -0.004   0.978     TRUE 0.509 1.002 52500
# Sigma.B[5,1]  -0.027 0.761  -1.035  -0.011   0.879     TRUE 0.530 1.013 17233
# Sigma.B[6,1]   0.007 0.812  -0.970   0.001   0.994     TRUE 0.502 1.005 20093
# Sigma.B[1,2]   0.039 0.748  -0.882   0.017   1.084     TRUE 0.543 1.001 52500
# Sigma.B[2,2]   0.645 1.165   0.106   0.376   2.798    FALSE 1.000 1.002 52500
# Sigma.B[3,2]  -0.019 0.716  -1.015  -0.008   0.893     TRUE 0.522 1.002 52500
# Sigma.B[4,2]   0.014 0.832  -0.965   0.003   1.018     TRUE 0.508 1.014 20291
# Sigma.B[5,2]  -0.012 0.721  -0.977  -0.006   0.901     TRUE 0.516 1.000 42566
# Sigma.B[6,2]   0.035 0.763  -0.878   0.009   1.090     TRUE 0.524 1.003 31358
# Sigma.B[1,3]  -0.008 0.753  -0.974  -0.003   0.938     TRUE 0.507 1.007 52500
# Sigma.B[2,3]  -0.019 0.716  -1.015  -0.008   0.893     TRUE 0.522 1.002 52500
# Sigma.B[3,3]   0.628 1.199   0.104   0.366   2.655    FALSE 1.000 1.002 52500
# Sigma.B[4,3]  -0.006 0.777  -0.975  -0.004   0.953     TRUE 0.509 1.010  6041
# Sigma.B[5,3]   0.022 0.735  -0.869   0.006   0.984     TRUE 0.517 1.002 21354
# Sigma.B[6,3]  -0.023 0.745  -1.042  -0.007   0.905     TRUE 0.519 1.008 10975
# Sigma.B[1,4]  -0.009 0.783  -1.017  -0.004   0.978     TRUE 0.509 1.002 52500
# Sigma.B[2,4]   0.014 0.832  -0.965   0.003   1.018     TRUE 0.508 1.014 20291
# Sigma.B[3,4]  -0.006 0.777  -0.975  -0.004   0.953     TRUE 0.509 1.010  6041
# Sigma.B[4,4]   0.688 1.493   0.110   0.398   2.899    FALSE 1.000 1.009 52500
# Sigma.B[5,4]  -0.031 0.812  -1.067  -0.012   0.906     TRUE 0.532 1.012 23944
# Sigma.B[6,4]   0.013 0.790  -0.985   0.003   1.064     TRUE 0.508 1.000 52500
# Sigma.B[1,5]  -0.027 0.761  -1.035  -0.011   0.879     TRUE 0.530 1.013 17233
# Sigma.B[2,5]  -0.012 0.721  -0.977  -0.006   0.901     TRUE 0.516 1.000 42566
# Sigma.B[3,5]   0.022 0.735  -0.869   0.006   0.984     TRUE 0.517 1.002 21354
# Sigma.B[4,5]  -0.031 0.812  -1.067  -0.012   0.906     TRUE 0.532 1.012 23944
# Sigma.B[5,5]   0.626 1.308   0.105   0.370   2.614    FALSE 1.000 1.061 12358
# Sigma.B[6,5]  -0.029 0.923  -1.037  -0.010   0.889     TRUE 0.525 1.034 52500
# Sigma.B[1,6]   0.007 0.812  -0.970   0.001   0.994     TRUE 0.502 1.005 20093
# Sigma.B[2,6]   0.035 0.763  -0.878   0.009   1.090     TRUE 0.524 1.003 31358
# Sigma.B[3,6]  -0.023 0.745  -1.042  -0.007   0.905     TRUE 0.519 1.008 10975
# Sigma.B[4,6]   0.013 0.790  -0.985   0.003   1.064     TRUE 0.508 1.000 52500
# Sigma.B[5,6]  -0.029 0.923  -1.037  -0.010   0.889     TRUE 0.525 1.034 52500
# Sigma.B[6,6]   0.687 1.576   0.108   0.391   2.982    FALSE 1.000 1.031 52500
# deviance     137.246 4.811 129.612 136.623 148.358    FALSE 1.000 1.001  4057
# 
# Successful convergence based on Rhat values (all < 1.1). 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 11.6 and DIC = 148.814 
# DIC is an estimate of expected predictive error (lower is better).

##run it to 80% again 

##by species here we go

quantile(out3$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.

#CNAT
quantile(out3$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

#OFAV
quantile(out3$sims.list$BB[ ,2,1], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,2,2], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,2,3], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,2,4], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,2,5], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,2,6], probs = c(0.1, 0.9))

##PSTR
quantile(out3$sims.list$BB[ ,3,1], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,3,2], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,3,3], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,3,4], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,3,5], probs = c(0.1, 0.9))

quantile(out3$sims.list$BB[ ,3,6], probs = c(0.1, 0.9))




##make a logistic regression figure for all
View(df6$species)
table(df6$species, df6$condition_h)
df7<-na.omit(df6)


f1 <- ggplot(df6, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "PSTR")), ncol = 9)

f1

f2 <- ggplot(df6, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "PSTR")), ncol = 9)
f2

f3 <- ggplot(df6, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "PSTR")), ncol = 9)

f3

f4 <- ggplot(df6, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "PSTR")), ncol = 9)

f4

f5 <- ggplot(df6, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "PSTR")), ncol = 9)
f5

f1/f2/f3/f4/f5


##########################################
#Spring 2023 FGB Unknown (Histo_FGB_March)
##########################################

s23_measures <- read_csv("data/s23_Unk.csv")

#step is VERY important or you get a bunch of NAs

s23_measures$value <- as.numeric(s23_measures$value)
s23_measures$measurement <- as.factor(s23_measures$measurement)


##get one value each
df1Tall <- s23_measures %>% group_by(individual, species, condition_h, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

#new df
df4 <- df1Tall %>% pivot_wider(id_cols = c(individual, species, condition_h), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$individual <- paste(df4$species, df4$individual, df4$condition_h, sep = "_")

View(df4)


# Group condition
##HD are included as DD because of previous research (they behave like DD)
## lets also remove any possible NAs before we start 
df4$state <- ifelse(df4$condition_h %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)
df5 <- df4[,c(1,2,3,4, 11, 5, 22, 23, 26, 36)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so they're centered on zero and won't have huge outliers
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##complete.cases??


# Exopercent = x1
df6$x1 <- (df6$mean_prop_exo - mean(df6$mean_prop_exo))/sd(df6$mean_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$mean_gastro - mean(df6$mean_gastro))/sd(df6$mean_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$max_avg_symb - mean(df6$max_avg_symb))/sd(df6$max_avg_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #0.32
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #-0.04
cor(df6$mean_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.23
cor(df6$mean_gastro, df6$max_symb_vac, use='complete.obs') #0.22
cor(df6$mean_gastro, df6$mean_degr_symb, use='complete.obs') #0.15
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #0.15


##nothing is correlated!!!

# Enumerate species for indexing as random effect
df6$sp <- as.numeric(as.factor(df6$species))

# Number of species
J <- length(unique(df6$species))

# Number of parameters
##this is one more than params because of intercept
K <- 6

# Create identity matrix for Wishart dist'n
# (Number of parameters to estimate (K))
W <- diag(K)


# Load data
data4 <- list(y = df6$state, 
              group = df6$sp,
              n = dim(df6)[1], 
              J = J,
              x1 = df6$x1,
              x2 = df6$x2,
              x3 = df6$x3, 
              x4 = df6$x4,
              x5 = df6$x5,
              K = K, 
              W = W)


# Define the model in the BUGS language and write a text file
##
sink("model.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1 + exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i] + beta4[group[i]] * x4[i]	+ beta5[group[i]] * x5[i]# linear predictor    
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2] # Exo
    beta2[j] <- BB[j,3] # Gastro
    beta3[j] <- BB[j,4] # Zoox
    beta4[j] <- BB[j,5] # symb size
    beta5[j] <- BB[j,6] # degraded symb
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivariate normal
    
    BB.hat[j,1] <- mu.alpha
	  BB.hat[j,2] <- mu.beta1
	  BB.hat[j,3] <- mu.beta2
	  BB.hat[j,4] <- mu.beta3
	  BB.hat[j,5] <- mu.beta4
	  BB.hat[j,6] <- mu.beta5
	  # Use below for 2nd-level of model
    #BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    #BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    mu.beta4 ~ dnorm(0, 0.0001)
    mu.beta5 ~ dnorm(0, 0.0001)
    #gamma0b ~ dnorm(0, 0.0001)
    #gamma0b2 ~ dnorm(0, 0.0001)
    #gamma0b3 ~ dnorm(0, 0.0001)

    #gamma1b ~ dnorm(0, 0.0001)
    #gamma1b2 ~ dnorm(0, 0.0001)
    #gamma1b3 ~ dnorm(0, 0.0001)
 
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()


# Initial values
inits <- function(){
  list(mu.alpha = rnorm(1), 
       mu.beta1 = rnorm(1), 
       mu.beta2 = rnorm(1), 
       mu.beta3 = rnorm(1), 
       mu.beta4 = rnorm(1), 
       mu.beta5 = rnorm(1),
       BB = matrix(rnorm(J*K),nrow=J,ncol=K),
       Tau.B = rwish(K+1,diag(K))
       #gamma0b = rnorm(1), 
       #gamma1b=rnorm(1),
       #gamma0b2=rnorm(1),
       #gamma0b3=rnorm(1),
       #gamma1b2=rnorm(1),
       #gamma1b3=rnorm(1) 
  )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta1","mu.beta2","mu.beta3", "mu.beta4",
                "mu.beta5", "BB", "Sigma.B"
                #"gamma0b","gamma1b","gamma0b2",
                #"gamma0b3","gamma1b2","gamma1b3"
)


# MCMC settings
ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

data4$y <- as.numeric(data4$y)

out4 <- jags(data4, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out4)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 200 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.172 minutes at time 2026-01-13 13:01:38.116723.
# 
# mean    sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mu.alpha     -3.215 1.732 -6.983 -2.984 -0.418    FALSE 0.990 1.506     8
# mu.beta1      0.287 0.886 -1.475  0.295  1.980     TRUE 0.650 1.000 52500
# mu.beta2     -0.096 0.783 -1.650 -0.092  1.391     TRUE 0.554 1.001  1888
# mu.beta3     -0.558 0.902 -2.321 -0.558  1.180     TRUE 0.745 1.001  3263
# mu.beta4     -0.036 0.838 -1.699 -0.028  1.591     TRUE 0.515 1.001  2312
# mu.beta5     11.002 7.494  0.587  9.815 27.714    FALSE 0.991 1.666     7
# BB[1,1]      -3.178 1.711 -6.904 -2.920 -0.454    FALSE 0.991 1.539     8
# BB[2,1]      -3.105 1.782 -7.034 -2.871 -0.197    FALSE 0.983 1.456     8
# BB[3,1]      -3.369 1.692 -7.069 -3.140 -0.730    FALSE 0.998 1.543     8
# BB[1,2]       0.160 0.573 -0.972  0.163  1.281     TRUE 0.613 1.002  1203
# BB[2,2]      -0.297 0.943 -2.473 -0.226  1.329     TRUE 0.607 1.003   959
# BB[3,2]       1.007 0.848 -0.499  0.936  2.898     TRUE 0.898 1.006   427
# BB[1,3]       0.424 0.554 -0.563  0.386  1.622     TRUE 0.780 1.002  1076
# BB[2,3]      -0.268 0.946 -2.349 -0.216  1.498     TRUE 0.602 1.009   295
# BB[3,3]      -0.459 0.625 -1.755 -0.433  0.714     TRUE 0.775 1.003  1675
# BB[1,4]      -0.503 0.887 -2.252 -0.504  1.205     TRUE 0.721 1.001  1732
# BB[2,4]      -1.006 0.796 -2.700 -0.955  0.428     TRUE 0.908 1.001  1744
# BB[3,4]      -0.172 0.894 -1.942 -0.171  1.668     TRUE 0.577 1.001  2581
# BB[1,5]       0.140 0.796 -1.430  0.143  1.728     TRUE 0.570 1.003   746
# BB[2,5]       0.101 0.823 -1.501  0.082  1.769     TRUE 0.539 1.002  1109
# BB[3,5]      -0.360 0.825 -2.108 -0.321  1.138     TRUE 0.660 1.001  1751
# BB[1,6]      11.049 7.522  0.500  9.802 27.915    FALSE 0.991 1.669     7
# BB[2,6]      10.940 7.495  0.624  9.718 27.643    FALSE 0.994 1.652     7
# BB[3,6]      11.023 7.486  0.446  9.839 27.797    FALSE 0.987 1.670     7
# Sigma.B[1,1]  0.916 2.525  0.120  0.481  4.205    FALSE 1.000 1.009 25561
# Sigma.B[2,1] -0.045 1.685 -2.023 -0.022  1.791     TRUE 0.535 1.009 13740
# Sigma.B[3,1]  0.004 1.455 -1.747  0.011  1.744     TRUE 0.519 1.001 10632
# Sigma.B[4,1] -0.007 1.446 -1.655 -0.007  1.690     TRUE 0.513 1.004  3287
# Sigma.B[5,1]  0.018 1.460 -1.490  0.012  1.581     TRUE 0.523 1.027  2833
# Sigma.B[6,1] -0.097 1.839 -2.189 -0.015  1.474     TRUE 0.526 1.024  2313
# Sigma.B[1,2] -0.045 1.685 -2.023 -0.022  1.791     TRUE 0.535 1.009 13740
# Sigma.B[2,2]  1.394 3.204  0.140  0.671  6.935    FALSE 1.000 1.016  1970
# Sigma.B[3,2] -0.118 1.819 -2.532 -0.050  1.895     TRUE 0.571 1.005  1807
# Sigma.B[4,2]  0.244 1.675 -1.393  0.085  2.727     TRUE 0.618 1.007 15029
# Sigma.B[5,2] -0.191 1.537 -2.510 -0.060  1.370     TRUE 0.591 1.009  3315
# Sigma.B[6,2]  0.103 2.057 -2.095  0.011  2.770     TRUE 0.517 1.008  6519
# Sigma.B[1,3]  0.004 1.455 -1.747  0.011  1.744     TRUE 0.519 1.001 10632
# Sigma.B[2,3] -0.118 1.819 -2.532 -0.050  1.895     TRUE 0.571 1.005  1807
# Sigma.B[3,3]  1.139 2.509  0.133  0.583  5.356    FALSE 1.000 1.003  8452
# Sigma.B[4,3] -0.011 1.821 -1.941 -0.005  2.027     TRUE 0.508 1.010  2196
# Sigma.B[5,3]  0.102 1.467 -1.483  0.036  2.026     TRUE 0.558 1.001 34472
# Sigma.B[6,3]  0.066 1.740 -1.884  0.012  2.379     TRUE 0.519 1.017  4890
# Sigma.B[1,4] -0.007 1.446 -1.655 -0.007  1.690     TRUE 0.513 1.004  3287
# Sigma.B[2,4]  0.244 1.675 -1.393  0.085  2.727     TRUE 0.618 1.007 15029
# Sigma.B[3,4] -0.011 1.821 -1.941 -0.005  2.027     TRUE 0.508 1.010  2196
# Sigma.B[4,4]  1.102 2.754  0.128  0.562  5.179    FALSE 1.000 1.018  3620
# Sigma.B[5,4] -0.218 1.519 -2.454 -0.059  1.130     TRUE 0.599 1.008  5350
# Sigma.B[6,4]  0.065 1.904 -1.853  0.005  2.322     TRUE 0.508 1.034  5995
# Sigma.B[1,5]  0.018 1.460 -1.490  0.012  1.581     TRUE 0.523 1.027  2833
# Sigma.B[2,5] -0.191 1.537 -2.510 -0.060  1.370     TRUE 0.591 1.009  3315
# Sigma.B[3,5]  0.102 1.467 -1.483  0.036  2.026     TRUE 0.558 1.001 34472
# Sigma.B[4,5] -0.218 1.519 -2.454 -0.059  1.130     TRUE 0.599 1.008  5350
# Sigma.B[5,5]  0.966 2.206  0.123  0.515  4.446    FALSE 1.000 1.002 52500
# Sigma.B[6,5] -0.047 1.672 -2.041 -0.007  1.696     TRUE 0.512 1.030 14670
# Sigma.B[1,6] -0.097 1.839 -2.189 -0.015  1.474     TRUE 0.526 1.024  2313
# Sigma.B[2,6]  0.103 2.057 -2.095  0.011  2.770     TRUE 0.517 1.008  6519
# Sigma.B[3,6]  0.066 1.740 -1.884  0.012  2.379     TRUE 0.519 1.017  4890
# Sigma.B[4,6]  0.065 1.904 -1.853  0.005  2.322     TRUE 0.508 1.034  5995
# Sigma.B[5,6] -0.047 1.672 -2.041 -0.007  1.696     TRUE 0.512 1.030 14670
# Sigma.B[6,6]  1.288 3.459  0.127  0.576  6.717    FALSE 1.000 1.010  4293
# deviance     58.903 4.584 51.488 58.339 69.301    FALSE 1.000 1.000 52500
# 
# **WARNING** Rhat values indicate convergence failure. 
# Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
# For each parameter, n.eff is a crude measure of effective sample size. 
# 
# overlap0 checks if 0 falls in the parameter's 95% credible interval.
# f is the proportion of the posterior with the same sign as the mean;
# i.e., our confidence that the parameter is positive or negative.
# 
# DIC info: (pD = var(deviance)/2) 
# pD = 10.5 and DIC = 69.41 
# DIC is an estimate of expected predictive error (lower is better).
 
##run it to 80% again 

##by species here we go

quantile(out4$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.

#OFRA
quantile(out4$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

#PAST
quantile(out4$sims.list$BB[ ,2,1], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,2,2], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,2,3], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,2,4], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,2,5], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,2,6], probs = c(0.1, 0.9))

##PSTR
quantile(out4$sims.list$BB[ ,3,1], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,3,2], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,3,3], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,3,4], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,3,5], probs = c(0.1, 0.9))

quantile(out4$sims.list$BB[ ,3,6], probs = c(0.1, 0.9))




##make a logistic regression figure for all
View(df6$species)
table(df6$species, df6$condition_h)
df7<-na.omit(df6)


f1 <- ggplot(df6, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f1

f2 <- ggplot(df6, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f2

f3 <- ggplot(df6, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f3

f4 <- ggplot(df6, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f4

f5 <- ggplot(df6, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f5

f1/f2/f3/f4/f5




