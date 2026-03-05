#####################
#WPH+FGB
#####################

library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)

#read in all data 

f22_measures <- read_csv("data/F22_Unk.csv")
s23_measures <- read_csv("data/s23_Unk.csv")
WPH_measures <- read_csv("data/WP_sonora.csv") 
SCTLD_measures <- read_csv("data/Serial_measures_MIAMI.csv")

##make sure all values are real

f22_measures$value <- as.numeric(f22_measures$value)
f22_measures$measurement <- as.factor(f22_measures$measurement)

s23_measures$value <- as.numeric(s23_measures$value)
s23_measures$measurement <- as.factor(s23_measures$measurement)

SCTLD_measures$value <- as.numeric(SCTLD_measures$value)
SCTLD_measures$measurement <- as.factor(SCTLD_measures$measurement)

WPH_measures$value <- as.numeric(WPH_measures$value)
WPH_measures$measurement <- as.factor(WPH_measures$measurement)

#Combine all datasets

combined_FGB <- bind_rows(f22_measures, s23_measures)

#clean SCTLD & WPH so columns match to combine + change condition_h to match in FGB

FGB_clean <- combined_FGB %>%
  mutate(
    health_status = condition_h
  ) %>%
  dplyr::select(
    individual,
    species,
    measurement,
    value,
    health_status
  ) %>% mutate(disease = "FGB")


     
WP_clean <- WPH_measures %>%
  mutate(
    health_status = worded_health
    
  ) %>%
  dplyr::select(
    individual = ID,
    species,
    measurement,
    value,
    health_status
  ) %>% mutate(disease = "WPH")


SCTLD_clean <- SCTLD_measures %>%
  mutate(
    health_status = worded_health
  ) %>%
  dplyr::select(
    individual = Individual,
    species,
    measurement,
    value,
    health_status
  ) %>% mutate(disease = "SCTLD")





#combine all data
combined <- bind_rows(
  FGB_clean,
  WP_clean,
  SCTLD_clean
)

#subset out model

WP_vs_FGB <- combined %>%
  filter(disease %in% c("WPH", "FGB"))



#########################
#lets do WP_vs_FGB first
#########################

##get one value each
df1Tall <- WP_vs_FGB %>% group_by(individual, species, health_status, measurement, disease) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

##making it wider and usable for the bayes model
df4 <- df1Tall %>% pivot_wider(id_cols = c(individual, species, health_status, disease), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 


df4$individual <- paste(df4$species, df4$individual, df4$health_status, df4$disease, sep = "_")

View(df4)

##lets make our disease states
df4$state <- ifelse(df4$health_status %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)

df5 <- df4[,c(1,2,3,4, 5, 8, 11, 12, 10, 22, 7, 17, 18, 36, 26, 23, 37, 27, 24)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so theyre centered on zero and won't have huge outlierts
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##zooxarea has been the name of vacuolization since the onset of this study and while it makes no sense, we keep it for nostalgia
# Exopercent = x1
df6$x1 <- (df6$max_prop_exo - mean(df6$max_prop_exo))/sd(df6$max_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$max_gastro - mean(df6$max_gastro))/sd(df6$max_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$min_min_symb - mean(df6$min_min_symb))/sd(df6$min_min_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


##what a great time to see if things are correlated
#values from -1 to 1, 0 means no linear relationship
cor(df6$max_prop_exo, df6$max_symb_vac, use='complete.obs') #0.14
cor(df6$max_prop_exo, df6$max_gastro, use='complete.obs') #0.02
cor(df6$max_prop_exo, df6$min_min_symb, use='complete.obs') #0.02
cor(df6$max_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.06
cor(df6$max_gastro, df6$max_symb_vac, use='complete.obs') #-0.19
cor(df6$max_gastro, df6$min_min_symb, use='complete.obs') #-0.17
cor(df6$max_gastro, df6$mean_degr_symb, use='complete.obs') #0.58
cor(df6$max_symb_vac, df6$min_min_symb, use='complete.obs') #0.57
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #-0.26
cor(df6$mean_degr_symb, df6$min_min_symb, use='complete.obs') #-0.12

##nothing is correlated!!!


##now we get into making things for the bayesian model itself
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

## model time!
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

## These can be changed as needed, up to what you need

ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 500 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.43 minutes at time 2026-01-26 20:55:10.117736.
# 
# mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# mu.alpha      -0.673 0.363  -1.389  -0.676   0.047     TRUE 0.967 1.002   823
# mu.beta1       0.026 0.350  -0.658   0.025   0.727     TRUE 0.530 1.001  3714
# mu.beta2       0.207 0.751  -1.182   0.167   1.796     TRUE 0.603 1.009   367
# mu.beta3      -0.270 0.489  -1.267  -0.262   0.673     TRUE 0.726 1.001  2190
# mu.beta4       0.215 0.405  -0.575   0.212   1.017     TRUE 0.716 1.001  1436
# mu.beta5      -0.538 0.609  -1.730  -0.540   0.699     TRUE 0.843 1.003  1130
# BB[1,1]       -0.699 0.416  -1.539  -0.693   0.113     TRUE 0.955 1.004   550
# BB[2,1]       -0.608 0.422  -1.435  -0.615   0.251     TRUE 0.924 1.002  1030
# BB[3,1]       -0.810 0.363  -1.536  -0.805  -0.121    FALSE 0.989 1.002  1379
# BB[4,1]       -0.596 0.505  -1.605  -0.594   0.385     TRUE 0.888 1.005   554
# BB[5,1]       -0.655 0.289  -1.238  -0.650  -0.106    FALSE 0.991 1.001  1994
# BB[1,2]       -0.151 0.368  -0.900  -0.141   0.542     TRUE 0.650 1.002  1027
# BB[2,2]       -0.179 0.355  -0.915  -0.172   0.491     TRUE 0.685 1.000  5308
# BB[3,2]        0.227 0.397  -0.545   0.218   1.021     TRUE 0.721 1.002  1676
# BB[4,2]        0.346 0.385  -0.352   0.322   1.180     TRUE 0.822 1.000 22476
# BB[5,2]       -0.117 0.273  -0.666  -0.113   0.408     TRUE 0.662 1.001  7561
# BB[1,3]        0.374 1.021  -1.492   0.311   2.644     TRUE 0.638 1.009   363
# BB[2,3]        0.804 0.875  -0.268   0.600   3.025     TRUE 0.891 1.029   147
# BB[3,3]       -0.406 0.543  -1.677  -0.337   0.445     TRUE 0.794 1.018   233
# BB[4,3]       -0.092 1.146  -2.685  -0.055   2.070     TRUE 0.522 1.002  2524
# BB[5,3]        0.367 0.895  -1.158   0.283   2.461     TRUE 0.646 1.017   198
# BB[1,4]       -0.282 0.552  -1.401  -0.273   0.781     TRUE 0.697 1.003   615
# BB[2,4]        0.405 0.513  -0.527   0.382   1.471     TRUE 0.783 1.002   958
# BB[3,4]       -0.482 0.511  -1.532  -0.473   0.488     TRUE 0.834 1.002  1866
# BB[4,4]       -1.162 0.578  -2.470  -1.097  -0.208    FALSE 0.994 1.001  3038
# BB[5,4]        0.166 0.317  -0.450   0.164   0.794     TRUE 0.701 1.001  1999
# BB[1,5]        0.458 0.420  -0.320   0.446   1.316     TRUE 0.867 1.001  2632
# BB[2,5]        0.405 0.489  -0.525   0.387   1.422     TRUE 0.799 1.002  1481
# BB[3,5]       -0.250 0.472  -1.242  -0.233   0.626     TRUE 0.702 1.003  1405
# BB[4,5]        0.481 0.474  -0.387   0.460   1.470     TRUE 0.858 1.004   489
# BB[5,5]       -0.010 0.339  -0.688  -0.004   0.653     TRUE 0.506 1.000 22376
# BB[1,6]       -0.891 0.778  -2.605  -0.833   0.479     TRUE 0.898 1.001  8640
# BB[2,6]       -1.212 0.652  -2.598  -1.169  -0.058    FALSE 0.982 1.003  1106
# BB[3,6]        0.183 0.533  -0.620   0.105   1.460     TRUE 0.602 1.022   223
# BB[4,6]        0.104 0.952  -1.441   0.005   2.321     TRUE 0.503 1.004   543
# BB[5,6]       -0.868 0.517  -1.990  -0.826   0.048     TRUE 0.967 1.003   665
# Sigma.B[1,1]   0.397 0.421   0.095   0.287   1.367    FALSE 1.000 1.004 10515
# Sigma.B[2,1]  -0.004 0.278  -0.512  -0.002   0.481     TRUE 0.508 1.004 52500
# Sigma.B[3,1]   0.002 0.495  -0.851   0.006   0.826     TRUE 0.515 1.011 11578
# Sigma.B[4,1]   0.034 0.447  -0.688   0.013   0.878     TRUE 0.531 1.002 32819
# Sigma.B[5,1]   0.014 0.323  -0.531   0.008   0.594     TRUE 0.527 1.003  4298
# Sigma.B[6,1]  -0.063 0.510  -1.031  -0.029   0.705     TRUE 0.566 1.003 11522
# Sigma.B[1,2]  -0.004 0.278  -0.512  -0.002   0.481     TRUE 0.508 1.004 52500
# Sigma.B[2,2]   0.424 0.462   0.099   0.304   1.470    FALSE 1.000 1.004 52500
# Sigma.B[3,2]  -0.119 0.535  -1.152  -0.054   0.556     TRUE 0.628 1.015  2000
# Sigma.B[4,2]  -0.168 0.459  -1.162  -0.093   0.391     TRUE 0.713 1.003 52500
# Sigma.B[5,2]  -0.009 0.320  -0.602  -0.005   0.550     TRUE 0.514 1.000 16875
# Sigma.B[6,2]   0.173 0.537  -0.466   0.090   1.303     TRUE 0.694 1.002 23142
# Sigma.B[1,3]   0.002 0.495  -0.851   0.006   0.826     TRUE 0.515 1.011 11578
# Sigma.B[2,3]  -0.119 0.535  -1.152  -0.054   0.556     TRUE 0.628 1.015  2000
# Sigma.B[3,3]   1.087 2.239   0.132   0.567   5.214    FALSE 1.000 1.066   376
# Sigma.B[4,3]   0.278 0.986  -0.821   0.124   2.230     TRUE 0.679 1.024  1189
# Sigma.B[5,3]   0.097 0.614  -0.797   0.051   1.239     TRUE 0.605 1.018  1886
# Sigma.B[6,3]  -0.482 1.402  -3.306  -0.185   0.561     TRUE 0.758 1.041   346
# Sigma.B[1,4]   0.034 0.447  -0.688   0.013   0.878     TRUE 0.531 1.002 32819
# Sigma.B[2,4]  -0.168 0.459  -1.162  -0.093   0.391     TRUE 0.713 1.003 52500
# Sigma.B[3,4]   0.278 0.986  -0.821   0.124   2.230     TRUE 0.679 1.024  1189
# Sigma.B[4,4]   0.854 0.999   0.148   0.569   3.366    FALSE 1.000 1.002  5319
# Sigma.B[5,4]  -0.094 0.523  -1.153  -0.042   0.672     TRUE 0.588 1.001  3044
# Sigma.B[6,4]  -0.391 0.890  -2.411  -0.206   0.518     TRUE 0.776 1.002  1591
# Sigma.B[1,5]   0.014 0.323  -0.531   0.008   0.594     TRUE 0.527 1.003  4298
# Sigma.B[2,5]  -0.009 0.320  -0.602  -0.005   0.550     TRUE 0.514 1.000 16875
# Sigma.B[3,5]   0.097 0.614  -0.797   0.051   1.239     TRUE 0.605 1.018  1886
# Sigma.B[4,5]  -0.094 0.523  -1.153  -0.042   0.672     TRUE 0.588 1.001  3044
# Sigma.B[5,5]   0.524 0.582   0.110   0.367   1.886    FALSE 1.000 1.001 40386
# Sigma.B[6,5]  -0.098 0.578  -1.247  -0.055   0.790     TRUE 0.607 1.001  8584
# Sigma.B[1,6]  -0.063 0.510  -1.031  -0.029   0.705     TRUE 0.566 1.003 11522
# Sigma.B[2,6]   0.173 0.537  -0.466   0.090   1.303     TRUE 0.694 1.002 23142
# Sigma.B[3,6]  -0.482 1.402  -3.306  -0.185   0.561     TRUE 0.758 1.041   346
# Sigma.B[4,6]  -0.391 0.890  -2.411  -0.206   0.518     TRUE 0.776 1.002  1591
# Sigma.B[5,6]  -0.098 0.578  -1.247  -0.055   0.790     TRUE 0.607 1.001  8584
# Sigma.B[6,6]   1.085 1.672   0.147   0.628   4.808    FALSE 1.000 1.003  1255
# deviance     226.639 6.391 215.818 226.087 240.625    FALSE 1.000 1.001  1873
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
# pD = 20.4 and DIC = 247.041 
# DIC is an estimate of expected predictive error (lower is better).




#back again at 80%
##by species here we go
#CNAT
quantile(out1$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.
quantile(out1$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

#OFAV
quantile(out1$sims.list$BB[ ,2,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,6], probs = c(0.1, 0.9))

#OFRA
quantile(out1$sims.list$BB[ ,3,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,6], probs = c(0.1, 0.9))

#PAST
quantile(out1$sims.list$BB[ ,4,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,6], probs = c(0.1, 0.9))

##PSTR
quantile(out1$sims.list$BB[ ,5,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,6], probs = c(0.1, 0.9))


View(df6$species)
table(df6$species, df6$health_status)
df7<-na.omit(df6)


######################
#Viewing WPH vs FGB
######################

##make a figure for everything
# 1️⃣ Define your species vector
species_names <- c("CNAT", "OFAV", "OFRA", "PAST",
                   "PSTR")

colors <- setNames(hcl.colors(length(species_names), "Spectral"), species_names)

q1 <- ggplot(df7, aes(x = mean_prop_exo, y = state, color =species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = colors) +
  scale_x_continuous(limits=c(0.0,0.7)) + theme(legend.position = "none") +
  xlab("Proportion Exocytosis") + ylab("Disease State")

q1

q2 <- ggplot(df7, aes(x = max_symb_vac, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  xlab("Vacuolization") + ylab("Disease State")

q2

q3 <- ggplot(df7, aes(x = mean_gastro, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = colors) +
  scale_x_continuous(limits=c(0, 50)) + theme(legend.position = "none") +
  xlab("Gastrodermal Separation (microns)") + theme(axis.title.y = element_blank())


q4 <- ggplot(df7, aes(x = max_avg_symb, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  xlab("Maximum Symbiont Size (microns)") + theme(axis.title.y = element_blank())
q4

(q2 | q4) / (q1 | q3)

library(patchwork)

combined_plot <- (q2 | q4) / (q1 | q3) +
  plot_layout(guides = "collect") +      # collect shared legend
  plot_annotation() & 
  theme(legend.position = "right",       # put legend to right
        legend.justification = "center") # vertically centered
combined_plot

df7_orb$disease <- ifelse(df7_orb$disease == "FGB", "FGB", "WPH")
df7_vac$disease <- ifelse(df7_vac$disease == "FGB","FGB", "WPH")
df7_sub$disease <- ifelse(df7_sub$disease == "FGB", "FGB", "WPH")

df7

r1 <- ggplot(df7_orb, aes(x = mean_prop_exo, y = state, color = species, shape = disease)) +
  geom_point(alpha = 0.7, size = 3, position = position_jitter(width=0.02, height=0)) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8, aes(linetype=disease)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  scale_shape_manual(values = c("FGB" = 16, "WP" = 17)) +
  xlab("Proportion Exocytosis") + ylab("Disease State") +
  theme(legend.position = "right")


r1

r2 <- ggplot(df7_vac, aes(x = max_symb_vac, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  theme(legend.position = "none") +
  xlab("Symbiont-To-Vacuole Ratio") + ylab("Disease State")

r2

r3 <- ggplot(df7_sub, aes(x = mean_gastro, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  scale_x_continuous(limits=c(0, 50)) + 
  xlab("Gastrodermal Separation (microns)") + theme(axis.title.y = element_blank())


r4 <- ggplot(df7_sub, aes(x = max_avg_symb, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  theme(legend.position = "none") +
  xlab("Maximum Symbiont Size (microns)") + theme(axis.title.y = element_blank())

(r2 | r1 | r4) +
  plot_layout(guides = "collect") +      # collect shared legend
  plot_annotation() & 
  theme(legend.position = "right",       # put legend to right
        legend.justification = "center") 
(r2 | r4) / (r1 | r3)
combined_plot2 <- (r2 | r4) / (r1 | r3) +
  plot_layout(guides = "collect") +      # collect shared legend
  plot_annotation() & 
  theme(legend.position = "right",       # put legend to right
        legend.justification = "center") # vertically centered
combined_plot2
