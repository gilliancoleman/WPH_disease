library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)

###############################
#Serial Fall 2021 SCTLD Samples
###############################

SCTLD_measures <- read_csv("data/Serial_measures_MIAMI.csv")



##ok this step is VERY important or you get a bunch of NAs, definitely something to remember

SCTLD_measures$value <- as.numeric(SCTLD_measures$value)
SCTLD_measures$measurement <- as.factor(SCTLD_measures$measurement)

##get one value each
df1Tall <- SCTLD_measures %>% group_by(Individual, species, worded_health, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

##I am number 4
df4 <- df1Tall %>% pivot_wider(id_cols = c(Individual, species, worded_health), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$ID <- paste(df4$species, df4$Individual, df4$worded_health, sep = "_")

View(df4)



# Group condition
##HD are included as DD because of previous research (they behave like DD idfk)
## lets also remove any possible NAs before we start so thats not the start of tears
df4$state <- ifelse(df4$worded_health %in% c("HD", "DD"), 0, 1)
View(df4)
df5 <- df4[,c(1,2,3,4, 11, 5, 22, 23, 26, 36, 37)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)
# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so theyre centered on zero and won't have huge outlierts
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


##what a great time to see if things are correlated
#values from -1 to 1, 0 means no linear relationship
cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #-0.08
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #0.302
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

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1)

# JAGS output for model 'model.txt', generated by jagsUI.
# Estimates based on 3 chains of 40000 iterations,
# adaptation = 400 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.092 minutes at time 2026-01-27 12:45:18.457652.
# 
# mean     sd   2.5%    50%  97.5% overlap0     f  Rhat n.eff
# mu.alpha     -2.966  3.330 -9.732 -2.195  1.786     TRUE 0.861 2.017     5
# mu.beta1      0.003  1.996 -3.140  0.019  3.121     TRUE 0.509 1.001 47112
# mu.beta2     -0.219  1.922 -3.311 -0.218  2.940     TRUE 0.589 1.001 14983
# mu.beta3      0.403  1.879 -2.662  0.402  3.506     TRUE 0.658 1.008 52500
# mu.beta4      0.725  2.011 -2.396  0.716  3.897     TRUE 0.766 1.001 16162
# mu.beta5     24.179 19.481  0.850 18.051 66.257    FALSE 0.989 3.014     4
# BB[1,1]      -2.974  2.795 -9.125 -2.065  0.466     TRUE 0.907 2.985     4
# BB[1,2]       0.029  0.346 -0.657  0.029  0.697     TRUE 0.534 1.003   768
# BB[1,3]      -0.215  0.422 -1.095 -0.200  0.571     TRUE 0.692 1.002  1039
# BB[1,4]       0.397  0.416 -0.353  0.371  1.278     TRUE 0.834 1.001  8874
# BB[1,5]       0.724  0.403 -0.030  0.708  1.556     TRUE 0.970 1.000  9287
# BB[1,6]      24.184 19.398  1.176 18.096 65.967    FALSE 0.999 3.064     4
# Sigma.B[1,1]  3.871 44.749  0.136  0.721 20.429    FALSE 1.000 1.056 52500
# Sigma.B[2,1] -0.230 37.463 -6.335  0.000  6.041     TRUE 0.500 1.087 18861
# Sigma.B[3,1] -0.091 25.467 -6.048  0.000  6.307     TRUE 0.500 1.038 31137
# Sigma.B[4,1] -0.131 31.192 -6.152  0.001  5.888     TRUE 0.499 1.074 52500
# Sigma.B[5,1] -0.162 35.121 -6.099 -0.005  5.889     TRUE 0.507 1.053 34410
# Sigma.B[6,1]  0.212 31.822 -6.104  0.001  6.324     TRUE 0.501 1.049 52500
# Sigma.B[1,2] -0.230 37.463 -6.335  0.000  6.041     TRUE 0.500 1.087 18861
# Sigma.B[2,2]  4.255 49.803  0.135  0.714 20.600    FALSE 1.000 1.015 52410
# Sigma.B[3,2]  0.295 28.379 -6.042  0.003  6.457     TRUE 0.504 1.070 24178
# Sigma.B[4,2] -0.006 32.863 -6.238  0.001  6.047     TRUE 0.499 1.064 52500
# Sigma.B[5,2]  0.007 40.528 -6.078 -0.001  6.127     TRUE 0.499 1.015 52500
# Sigma.B[6,2] -0.157 30.792 -6.401 -0.001  6.113     TRUE 0.501 1.002 35805
# Sigma.B[1,3] -0.091 25.467 -6.048  0.000  6.307     TRUE 0.500 1.038 31137
# Sigma.B[2,3]  0.295 28.379 -6.042  0.003  6.457     TRUE 0.504 1.070 24178
# Sigma.B[3,3]  3.965 40.423  0.135  0.728 19.608    FALSE 1.000 1.057 27658
# Sigma.B[4,3] -0.141 27.239 -6.168  0.002  6.060     TRUE 0.498 1.189 52500
# Sigma.B[5,3] -0.011 33.647 -6.129 -0.002  6.126     TRUE 0.503 1.051 52500
# Sigma.B[6,3] -0.061 33.445 -6.079  0.002  6.244     TRUE 0.498 1.062 52500
# Sigma.B[1,4] -0.131 31.192 -6.152  0.001  5.888     TRUE 0.499 1.074 52500
# Sigma.B[2,4] -0.006 32.863 -6.238  0.001  6.047     TRUE 0.499 1.064 52500
# Sigma.B[3,4] -0.141 27.239 -6.168  0.002  6.060     TRUE 0.498 1.189 52500
# Sigma.B[4,4]  3.917 46.554  0.134  0.718 20.118    FALSE 1.000 1.122  8059
# Sigma.B[5,4]  0.155 38.663 -6.126 -0.001  6.170     TRUE 0.499 1.069 52500
# Sigma.B[6,4]  0.045 30.315 -6.205  0.000  6.303     TRUE 0.500 1.065  7277
# Sigma.B[1,5] -0.162 35.121 -6.099 -0.005  5.889     TRUE 0.507 1.053 34410
# Sigma.B[2,5]  0.007 40.528 -6.078 -0.001  6.127     TRUE 0.499 1.015 52500
# Sigma.B[3,5] -0.011 33.647 -6.129 -0.002  6.126     TRUE 0.503 1.051 52500
# Sigma.B[4,5]  0.155 38.663 -6.126 -0.001  6.170     TRUE 0.499 1.069 52500
# Sigma.B[5,5]  4.484 74.019  0.136  0.725 20.126    FALSE 1.000 1.054 52500
# Sigma.B[6,5]  0.051 76.227 -6.425 -0.003  6.269     TRUE 0.497 1.133 26402
# Sigma.B[1,6]  0.212 31.822 -6.104  0.001  6.324     TRUE 0.501 1.049 52500
# Sigma.B[2,6] -0.157 30.792 -6.401 -0.001  6.113     TRUE 0.501 1.002 35805
# Sigma.B[3,6] -0.061 33.445 -6.079  0.002  6.244     TRUE 0.498 1.062 52500
# Sigma.B[4,6]  0.045 30.315 -6.205  0.000  6.303     TRUE 0.500 1.065  7277
# Sigma.B[5,6]  0.051 76.227 -6.425 -0.003  6.269     TRUE 0.497 1.133 26402
# Sigma.B[6,6]  4.499 94.124  0.134  0.722 20.615    FALSE 1.000 1.144 52500
# deviance     65.567  3.358 61.084 64.915 73.867    FALSE 1.000 1.001  4049
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
# pD = 5.6 and DIC = 71.203 
# DIC is an estimate of expected predictive error (lower is better).

## i think we run it to 80% again because I clearly hate myself

##by species here we go
##OFRA
quantile(out1$sims.list$BB[ ,1,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

View(df6$species)
table(df6$species, df6$worded_health)
df7<-na.omit(df6)


######################
#Viewing SCTLD
######################
species_names <- c("OFRA")

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
  filter(species %in% c("OFRA"))

df7_orb <- df7 %>%
  filter(species %in% c("OFRA"))

df7_vac <- df7 %>%
  filter(species %in% c("OFRA"))

df7

r1 <- ggplot(df7_orb, aes(x = mean_prop_exo, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e"
                                )) +
  scale_x_continuous(limits=c(0.0,0.7)) + theme(legend.position = "none") +
  xlab("Proportion Exocytosis") + ylab("Disease State")

r1

r2 <- ggplot(df7_vac, aes(x = max_symb_vac, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e"
                                )) +
  theme(legend.position = "none") +
  xlab("Symbiont-To-Vacuole Ratio") + ylab("Disease State")

r2

r3 <- ggplot(df7_sub, aes(x = mean_gastro, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("OFRA" = "#bf8f5e"
                                )) +
  scale_x_continuous(limits=c(0, 50)) + 
  xlab("Gastrodermal Separation (microns)") + theme(axis.title.y = element_blank())


r4 <- ggplot(df7_sub, aes(x = max_avg_symb, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e"
                                )) +
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


