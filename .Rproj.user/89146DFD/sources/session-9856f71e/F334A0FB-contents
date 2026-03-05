library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)

#################
#WPH
#################

WPH_measures <- read_csv("data/WP_sonora.csv") 


## we do all of this first to make sure all values are real

WPH_measures$value <- as.numeric(WPH_measures$value)
WPH_measures$measurement <- as.factor(WPH_measures$measurement)

##get one value each
df1Tall <- WPH_measures %>% group_by(ID, species, worded_health, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

##making it wider and usable for the bayes model
df4 <- df1Tall %>% pivot_wider(id_cols = c(ID, species, worded_health), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 


df4$ID <- paste(df4$species, df4$ID, df4$worded_health, sep = "_")

View(df4)

##lets make our disease states
df4$state <- ifelse(df4$worded_health %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)

df5 <- df4[,c(1,2,3,4, 11, 10, 22, 7, 17, 36, 26, 23)]
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
cor(df6$max_prop_exo, df6$max_symb_vac, use='complete.obs') #0.108
cor(df6$max_prop_exo, df6$max_gastro, use='complete.obs') #0.08
cor(df6$max_prop_exo, df6$min_min_symb, use='complete.obs') #0.011
cor(df6$max_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.154
cor(df6$max_gastro, df6$max_symb_vac, use='complete.obs') #-0.36
cor(df6$max_gastro, df6$min_min_symb, use='complete.obs') #-0.06
cor(df6$max_gastro, df6$mean_degr_symb, use='complete.obs') #0.79
cor(df6$max_symb_vac, df6$min_min_symb, use='complete.obs') #0.425
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #-0.46
cor(df6$mean_degr_symb, df6$min_min_symb, use='complete.obs') #-0.102

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
## if it is false, the mean doesn't cross 0, so it's what we would normally call significant, but here is "strongly correlated"

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 200 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.102 minutes at time 2026-01-27 12:35:56.500449.
# 
# mean      sd     2.5%    50%   97.5% overlap0     f  Rhat n.eff
# mu.alpha      -4.538   3.222  -11.869 -3.824  -0.015    FALSE 0.975 1.025    98
# mu.beta1       2.196   2.149   -0.825  1.763   7.277     TRUE 0.915 1.004  1078
# mu.beta2      -4.904   5.204  -16.773 -4.264   3.717     TRUE 0.878 1.274    12
# mu.beta3      -2.466   4.586  -12.597 -1.964   5.005     TRUE 0.829 1.003  2499
# mu.beta4      -1.491   2.932   -8.126 -1.070   2.838     TRUE 0.759 1.004  3243
# mu.beta5       4.708   5.790   -5.582  4.226  17.530     TRUE 0.821 1.187    18
# BB[1,1]       -4.621   2.423  -10.447 -4.085  -1.204    FALSE 1.000 1.051    52
# BB[2,1]       -4.438   3.819  -13.506 -3.328   1.123     TRUE 0.939 1.067    44
# BB[1,2]        1.726   1.603   -0.767  1.433   5.587     TRUE 0.896 1.023   911
# BB[2,2]        2.659   2.397   -0.062  1.867   9.290     TRUE 0.970 1.019   179
# BB[1,3]       -5.581   4.202  -16.394 -4.954  -0.031    FALSE 0.976 1.327    12
# BB[2,3]       -4.223   6.106  -17.880 -3.520   6.945     TRUE 0.777 1.347    10
# BB[1,4]        0.846   2.231   -2.838  0.546   5.685     TRUE 0.602 1.007  2786
# BB[2,4]       -5.821   4.279  -16.542 -4.344  -0.793    FALSE 0.997 1.014   256
# BB[1,5]       -3.473   2.937  -10.545 -2.699   0.437     TRUE 0.940 1.024  6632
# BB[2,5]        0.520   1.175   -1.559  0.424   3.022     TRUE 0.659 1.016   172
# BB[1,6]        5.376   4.968   -1.348  4.627  18.064     TRUE 0.896 1.307    13
# BB[2,6]        4.021   6.490   -9.872  4.068  17.195     TRUE 0.760 1.134    23
# Sigma.B[1,1]   6.683  43.537    0.151  1.225  41.196    FALSE 1.000 1.032 10131
# Sigma.B[2,1]  -1.329  21.479  -15.771 -0.061   6.165     TRUE 0.556 1.119   824
# Sigma.B[3,1]  -4.913  55.122  -42.771 -0.192   7.426     TRUE 0.615 1.115 52500
# Sigma.B[4,1]   1.319  98.773  -39.758 -0.063  49.985     TRUE 0.467 1.090   530
# Sigma.B[5,1]   0.754  42.807  -22.059  0.065  28.390     TRUE 0.540 1.024   470
# Sigma.B[6,1]  -3.918  56.438  -37.543 -0.125  10.598     TRUE 0.590 1.013  3236
# Sigma.B[1,2]  -1.329  21.479  -15.771 -0.061   6.165     TRUE 0.556 1.119   824
# Sigma.B[2,2]   3.461  17.776    0.143  0.840  21.249    FALSE 1.000 1.077   827
# Sigma.B[3,2]   2.328  35.600   -9.689  0.059  26.635     TRUE 0.548 1.166   997
# Sigma.B[4,2]  -5.922  60.808  -51.090 -0.306   9.419     TRUE 0.629 1.126  1072
# Sigma.B[5,2]   2.409  24.776   -7.656  0.193  24.505     TRUE 0.614 1.040  1143
# Sigma.B[6,2]  -1.246  29.257  -19.269 -0.027   9.671     TRUE 0.523 1.085  9589
# Sigma.B[1,3]  -4.913  55.122  -42.771 -0.192   7.426     TRUE 0.615 1.115 52500
# Sigma.B[2,3]   2.328  35.600   -9.689  0.059  26.635     TRUE 0.548 1.166   997
# Sigma.B[3,3]  15.805 109.413    0.159  1.836 106.459    FALSE 1.000 1.075  2958
# Sigma.B[4,3]  -9.428 160.177 -103.058 -0.183  42.154     TRUE 0.574 1.135   427
# Sigma.B[5,3]   5.086  70.866  -24.118  0.153  57.700     TRUE 0.576 1.026   302
# Sigma.B[6,3]  -0.574  77.268  -36.422  0.006  33.436     TRUE 0.496 1.066   649
# Sigma.B[1,4]   1.319  98.773  -39.758 -0.063  49.985     TRUE 0.467 1.090   530
# Sigma.B[2,4]  -5.922  60.808  -51.090 -0.306   9.419     TRUE 0.629 1.126  1072
# Sigma.B[3,4]  -9.428 160.177 -103.058 -0.183  42.154     TRUE 0.574 1.135   427
# Sigma.B[4,4]  39.012 303.626    0.213  6.320 239.071    FALSE 1.000 1.070 31018
# Sigma.B[5,4] -20.593 140.936 -131.166 -3.260   0.790     TRUE 0.878 1.017 16030
# Sigma.B[6,4]   9.086 136.151  -29.737  0.142  94.331     TRUE 0.559 1.018  8844
# Sigma.B[1,5]   0.754  42.807  -22.059  0.065  28.390     TRUE 0.540 1.024   470
# Sigma.B[2,5]   2.409  24.776   -7.656  0.193  24.505     TRUE 0.614 1.040  1143
# Sigma.B[3,5]   5.086  70.866  -24.118  0.153  57.700     TRUE 0.576 1.026   302
# Sigma.B[4,5] -20.593 140.936 -131.166 -3.260   0.790     TRUE 0.878 1.017 16030
# Sigma.B[5,5]  14.191  85.902    0.188  2.825  86.540    FALSE 1.000 1.040 12325
# Sigma.B[6,5]  -4.294  62.387  -49.900 -0.107  18.658     TRUE 0.558 1.013  5484
# Sigma.B[1,6]  -3.918  56.438  -37.543 -0.125  10.598     TRUE 0.590 1.013  3236
# Sigma.B[2,6]  -1.246  29.257  -19.269 -0.027   9.671     TRUE 0.523 1.085  9589
# Sigma.B[3,6]  -0.574  77.268  -36.422  0.006  33.436     TRUE 0.496 1.066   649
# Sigma.B[4,6]   9.086 136.151  -29.737  0.142  94.331     TRUE 0.559 1.018  8844
# Sigma.B[5,6]  -4.294  62.387  -49.900 -0.107  18.658     TRUE 0.558 1.013  5484
# Sigma.B[6,6]  13.279 107.231    0.158  1.499  87.486    FALSE 1.000 1.011 33394
# deviance      26.193   5.805   17.103 25.600  38.949    FALSE 1.000 1.007   307
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
# pD = 16.7 and DIC = 42.932 
# DIC is an estimate of expected predictive error (lower is better).

View(df6$species)
table(df6$species, df6$worded_health)
df7<-na.omit(df6)


f1 <- ggplot(df7, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f1

f2 <- ggplot(df7, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f2

f3 <- ggplot(df7, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f3

f4 <- ggplot(df7, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)

f4

f5 <- ggplot(df7, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("OFRA", "PAST", "PSTR")), ncol = 9)
f5

f1/f2/f3/f4/f5

##by species here we go
#OFRA
quantile(out1$sims.list$BB[ ,1,1], probs = c(0.1, 0.9)) #In Bayesian terms, this gives you an 80% credible interval for the parameter.

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



######################
#Viewing WPH
######################

##make a figure for everything
# 1️⃣ Define your species vector
species_names <- c("PAST",
                   "OFRA")

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

r1 <- ggplot(df7_orb, aes(x = mean_prop_exo, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554")) +
  scale_x_continuous(limits=c(0.0,0.7)) + theme(legend.position = "none") +
  xlab("Proportion Exocytosis") + ylab("Disease State")

r1

r2 <- ggplot(df7_vac, aes(x = max_symb_vac, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554")) +
  theme(legend.position = "none") +
  xlab("Symbiont-To-Vacuole Ratio") + ylab("Disease State")

r2

r3 <- ggplot(df7_sub, aes(x = mean_gastro, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554")) +
  scale_x_continuous(limits=c(0, 50)) + 
  xlab("Gastrodermal Separation (microns)") + theme(axis.title.y = element_blank())


r4 <- ggplot(df7_sub, aes(x = max_avg_symb, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554")) +
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



