#####################
#SCTLD+FGB 
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
SCTLD_vs_FGB <- combined %>%
  filter(disease %in% c("SCTLD", "FGB"))

#########################
#now do SCTLD_vs_FGB 
#########################

##get one value each
df1Tall <- SCTLD_vs_FGB %>% group_by(individual, species, health_status, measurement, disease) %>%
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
cor(df6$max_prop_exo, df6$max_symb_vac, use='complete.obs') #-0.04
cor(df6$max_prop_exo, df6$max_gastro, use='complete.obs') #0.02
cor(df6$max_prop_exo, df6$min_min_symb, use='complete.obs') #0.05
cor(df6$max_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.07
cor(df6$max_gastro, df6$max_symb_vac, use='complete.obs') #0.002
cor(df6$max_gastro, df6$min_min_symb, use='complete.obs') #-0.11
cor(df6$max_gastro, df6$mean_degr_symb, use='complete.obs') #-0.02
cor(df6$max_symb_vac, df6$min_min_symb, use='complete.obs') #0.38
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #0.02
cor(df6$mean_degr_symb, df6$min_min_symb, use='complete.obs') #0.16

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
# adaptation = 400 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.34 minutes at time 2026-01-27 09:46:11.260969.
# 
# mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# mu.alpha      -1.334 0.599  -2.574  -1.268  -0.312    FALSE 0.994 1.776     6
# mu.beta1      -0.137 0.343  -0.818  -0.134   0.534     TRUE 0.668 1.000 33302
# mu.beta2      10.196 7.017   0.483   7.909  25.746    FALSE 0.984 2.580     4
# mu.beta3      -0.107 0.428  -0.981  -0.102   0.717     TRUE 0.609 1.001  3037
# mu.beta4       0.452 0.363  -0.249   0.449   1.171     TRUE 0.905 1.001  3405
# mu.beta5      -0.185 0.478  -1.067  -0.205   0.820     TRUE 0.691 1.001  2887
# BB[1,1]       -1.444 0.603  -2.718  -1.377  -0.414    FALSE 0.998 1.788     6
# BB[2,1]       -1.458 0.613  -2.742  -1.400  -0.406    FALSE 0.998 1.692     6
# BB[3,1]       -0.993 0.560  -2.168  -0.922  -0.047    FALSE 0.981 1.929     5
# BB[4,1]       -1.313 0.745  -2.840  -1.276   0.061     TRUE 0.970 1.448     8
# BB[5,1]       -1.454 0.549  -2.620  -1.370  -0.558    FALSE 1.000 1.978     5
# BB[1,2]       -0.225 0.358  -0.950  -0.216   0.456     TRUE 0.739 1.000 12766
# BB[2,2]       -0.173 0.348  -0.900  -0.160   0.473     TRUE 0.685 1.000  7380
# BB[3,2]       -0.010 0.277  -0.552  -0.012   0.527     TRUE 0.516 1.002  1245
# BB[4,2]       -0.114 0.560  -1.237  -0.114   0.979     TRUE 0.591 1.001 14632
# BB[5,2]       -0.162 0.273  -0.713  -0.159   0.356     TRUE 0.719 1.001  2820
# BB[1,3]       10.216 7.037   0.418   7.997  25.806    FALSE 0.983 2.574     4
# BB[2,3]       10.227 6.976   0.713   8.072  25.491    FALSE 0.997 2.588     4
# BB[3,3]       10.200 7.036   0.369   7.908  25.653    FALSE 0.980 2.555     4
# BB[4,3]       10.135 7.136   0.085   7.988  25.874    FALSE 0.976 2.537     4
# BB[5,3]       10.192 6.990   0.576   7.967  25.648    FALSE 0.985 2.568     4
# BB[1,4]       -0.163 0.495  -1.147  -0.157   0.795     TRUE 0.627 1.001 52500
# BB[2,4]        0.413 0.461  -0.423   0.385   1.396     TRUE 0.826 1.002  1558
# BB[3,4]       -0.113 0.338  -0.765  -0.118   0.562     TRUE 0.639 1.001  4210
# BB[4,4]       -0.827 0.538  -2.023  -0.780   0.077     TRUE 0.961 1.004   554
# BB[5,4]        0.153 0.271  -0.388   0.154   0.696     TRUE 0.715 1.001  3722
# BB[1,5]        0.475 0.360  -0.207   0.469   1.208     TRUE 0.911 1.000  6150
# BB[2,5]        0.538 0.423  -0.254   0.522   1.423     TRUE 0.907 1.001  8155
# BB[3,5]        0.694 0.325   0.075   0.683   1.358    FALSE 0.987 1.001  3805
# BB[4,5]        0.479 0.508  -0.472   0.459   1.561     TRUE 0.838 1.003  1318
# BB[5,5]        0.076 0.291  -0.508   0.081   0.643     TRUE 0.610 1.000 18695
# BB[1,6]       -0.564 0.532  -1.742  -0.519   0.378     TRUE 0.874 1.004   694
# BB[2,6]       -0.688 0.345  -1.411  -0.669  -0.076    FALSE 0.987 1.001  3727
# BB[3,6]        0.423 0.503  -0.388   0.360   1.566     TRUE 0.809 1.003  1329
# BB[4,6]        0.396 0.737  -0.724   0.293   2.237     TRUE 0.702 1.009   602
# BB[5,6]       -0.486 0.292  -1.135  -0.458   0.020     TRUE 0.969 1.000  9955
# Sigma.B[1,1]   0.426 0.454   0.099   0.306   1.463    FALSE 1.000 1.001 16372
# Sigma.B[2,1]   0.043 0.277  -0.388   0.024   0.586     TRUE 0.587 1.001 17316
# Sigma.B[3,1]  -0.031 0.435  -0.857  -0.009   0.660     TRUE 0.526 1.003 12113
# Sigma.B[4,1]  -0.020 0.388  -0.721  -0.014   0.639     TRUE 0.537 1.000 14309
# Sigma.B[5,1]   0.014 0.292  -0.497   0.011   0.540     TRUE 0.538 1.000 52500
# Sigma.B[6,1]   0.085 0.435  -0.568   0.054   0.893     TRUE 0.639 1.002 52500
# Sigma.B[1,2]   0.043 0.277  -0.388   0.024   0.586     TRUE 0.587 1.001 17316
# Sigma.B[2,2]   0.383 0.399   0.091   0.276   1.317    FALSE 1.000 1.002 52500
# Sigma.B[3,2]   0.003 0.432  -0.693   0.002   0.712     TRUE 0.506 1.023  5020
# Sigma.B[4,2]  -0.013 0.369  -0.669  -0.007   0.612     TRUE 0.520 1.002 19722
# Sigma.B[5,2]   0.011 0.284  -0.463   0.007   0.524     TRUE 0.524 1.002 52500
# Sigma.B[6,2]   0.035 0.403  -0.622   0.021   0.771     TRUE 0.559 1.004  5879
# Sigma.B[1,3]  -0.031 0.435  -0.857  -0.009   0.660     TRUE 0.526 1.003 12113
# Sigma.B[2,3]   0.003 0.432  -0.693   0.002   0.712     TRUE 0.506 1.023  5020
# Sigma.B[3,3]   0.785 1.524   0.116   0.444   3.543    FALSE 1.000 1.107   979
# Sigma.B[4,3]   0.016 0.740  -1.027   0.006   1.097     TRUE 0.512 1.066  1001
# Sigma.B[5,3]   0.006 0.458  -0.729   0.002   0.792     TRUE 0.507 1.006  2939
# Sigma.B[6,3]  -0.046 0.969  -1.344  -0.006   1.091     TRUE 0.511 1.117  1011
# Sigma.B[1,4]  -0.020 0.388  -0.721  -0.014   0.639     TRUE 0.537 1.000 14309
# Sigma.B[2,4]  -0.013 0.369  -0.669  -0.007   0.612     TRUE 0.520 1.002 19722
# Sigma.B[3,4]   0.016 0.740  -1.027   0.006   1.097     TRUE 0.512 1.066  1001
# Sigma.B[4,4]   0.635 0.789   0.119   0.429   2.372    FALSE 1.000 1.003 44367
# Sigma.B[5,4]  -0.063 0.401  -0.841  -0.029   0.518     TRUE 0.580 1.000 21614
# Sigma.B[6,4]  -0.224 0.634  -1.506  -0.114   0.430     TRUE 0.731 1.024  3962
# Sigma.B[1,5]   0.014 0.292  -0.497   0.011   0.540     TRUE 0.538 1.000 52500
# Sigma.B[2,5]   0.011 0.284  -0.463   0.007   0.524     TRUE 0.524 1.002 52500
# Sigma.B[3,5]   0.006 0.458  -0.729   0.002   0.792     TRUE 0.507 1.006  2939
# Sigma.B[4,5]  -0.063 0.401  -0.841  -0.029   0.518     TRUE 0.580 1.000 21614
# Sigma.B[5,5]   0.436 0.482   0.102   0.312   1.505    FALSE 1.000 1.001 31943
# Sigma.B[6,5]   0.054 0.446  -0.619   0.031   0.861     TRUE 0.578 1.002 11929
# Sigma.B[1,6]   0.085 0.435  -0.568   0.054   0.893     TRUE 0.639 1.002 52500
# Sigma.B[2,6]   0.035 0.403  -0.622   0.021   0.771     TRUE 0.559 1.004  5879
# Sigma.B[3,6]  -0.046 0.969  -1.344  -0.006   1.091     TRUE 0.511 1.117  1011
# Sigma.B[4,6]  -0.224 0.634  -1.506  -0.114   0.430     TRUE 0.731 1.024  3962
# Sigma.B[5,6]   0.054 0.446  -0.619   0.031   0.861     TRUE 0.578 1.002 11929
# Sigma.B[6,6]   0.783 1.117   0.131   0.490   3.229    FALSE 1.000 1.020  1425
# deviance     261.493 5.931 251.607 260.869 274.623    FALSE 1.000 1.002  1043
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
# pD = 17.6 and DIC = 279.049 
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
#Viewing STCLD vs FGB
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

df7_sub <- df7 %>%
  filter(species %in% c("OFRA", "PAST"))

df7_orb <- df7 %>%
  filter(species %in% c("OFRA", "PAST"))

df7_vac <- df7 %>%
  filter(species %in% c("OFRA", "PAST"))

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

r2 <- ggplot(df7_vac, aes(x = max_symb_vac, y = state, color = species, shape = disease)) +
  geom_point(alpha = 0.7, size = 3, position = position_jitter(width=0.02, height=0)) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8, aes(linetype=disease)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  scale_shape_manual(values = c("FGB" = 16, "WP" = 17)) +
  xlab("Symbiont-To-Vacuole Ratio") + ylab("Disease State") +
  theme(legend.position = "right")

r2

r3 <- ggplot(df7_sub, aes(x = mean_gastro, y = state, color = species, shape = disease)) +
  geom_point(alpha = 0.7, size = 3, position = position_jitter(width=0.02, height=0)) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8, aes(linetype=disease)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  scale_shape_manual(values = c("FGB" = 16, "WP" = 17)) +
  xlab("Gastrodermal Separation (microns)") + ylab("Disease State") +
  theme(legend.position = "right")

r4 <- ggplot(df7_sub, aes(x = max_avg_symb, y = state, color = species, shape = disease)) +
  geom_point(alpha = 0.7, size = 3, position = position_jitter(width=0.02, height=0)) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8, aes(linetype=disease)) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580")) +
  scale_shape_manual(values = c("FGB" = 16, "WP" = 17)) +
  xlab("Maximum Symbiont Size (microns)") + ylab("Disease State") +
  theme(legend.position = "right")

(r2 | r1 | r4) +
  plot_layout(guides = "collect") +      # collect shared legend
  plot_annotation() & 
  theme(legend.position = "right",       # put legend to right
        legend.justification = "center") 
(r2 | r4) / (r1 | r3)
pcombined_plot2 <- (r2 | r4) / (r1 | r3) +
  plot_layout(guides = "collect") +      # collect shared legend
  plot_annotation() & 
  theme(legend.position = "right",       # put legend to right
        legend.justification = "center") # vertically centered
combined_plot2
