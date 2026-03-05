library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)


##############
#Unknown
##############

f22_measures <- read_csv("data/F22_Unk.csv")
s23_measures <- read_csv("data/s23_Unk.csv")


## we do all of this first to make sure all values are real

f22_measures$value <- as.numeric(f22_measures$value)
f22_measures$measurement <- as.factor(f22_measures$measurement)

s23_measures$value <- as.numeric(s23_measures$value)
s23_measures$measurement <- as.factor(s23_measures$measurement)

#combine data
all_measures <- bind_rows(f22_measures, s23_measures)

##get one value each
df1Tall <- all_measures %>% group_by(individual, species, condition_h, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

#new df
df4 <- df1Tall %>% pivot_wider(id_cols = c(individual, species, condition_h), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$ID <- paste(df4$species, df4$individual, df4$condition_h, sep = "_")

View(df4)


##lets make our disease states
df4$state <- ifelse(df4$condition_h %in% c("HD", "DD"), 0, 1) #group HD/DD and assign a 1 for both
View(df4)

df5 <- df4[,c(1,2,3,4, 11, 10, 22, 7, 17, 36, 37, 26, 23)]
na.omit(df5)
View(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)

# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so theyre centered on zero and won't have huge outliers
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
cor(df6$max_prop_exo, df6$max_gastro, use='complete.obs') #0.007
cor(df6$max_prop_exo, df6$min_min_symb, use='complete.obs') #-0.007
cor(df6$max_prop_exo, df6$mean_degr_symb, use='complete.obs') #0.04
cor(df6$max_gastro, df6$max_symb_vac, use='complete.obs') #0.02
cor(df6$max_gastro, df6$min_min_symb, use='complete.obs') #-0.13
cor(df6$max_gastro, df6$mean_degr_symb, use='complete.obs') #-0.02
cor(df6$max_symb_vac, df6$min_min_symb, use='complete.obs') #0.47
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #0.06
cor(df6$mean_degr_symb, df6$min_min_symb, use='complete.obs') #0.15

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
# adaptation = 500 iterations (sufficient),
# burn-in = 5000 iterations and thin rate = 2,
# yielding 52500 total samples from the joint posterior. 
# MCMC ran for 0.295 minutes at time 2026-01-27 12:53:17.118342.
# 
# mean    sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# mu.alpha      -1.864 0.781  -3.373  -1.841  -0.540    FALSE 0.998 2.020     5
# mu.beta1      -0.104 0.373  -0.836  -0.107   0.641     TRUE 0.623 1.000 36151
# mu.beta2      13.761 8.852   1.303  13.589  30.489    FALSE 0.997 2.568     4
# mu.beta3      -0.129 0.460  -1.080  -0.121   0.748     TRUE 0.616 1.002  1295
# mu.beta4       0.336 0.395  -0.421   0.329   1.128     TRUE 0.819 1.001  2036
# mu.beta5      -0.318 0.464  -1.197  -0.332   0.635     TRUE 0.787 1.001  2195
# BB[1,1]       -1.882 0.777  -3.373  -1.849  -0.558    FALSE 0.999 2.016     5
# BB[2,1]       -1.845 0.786  -3.351  -1.812  -0.477    FALSE 0.999 1.983     5
# BB[3,1]       -1.873 0.796  1-3.42  -1.845  -0.497    FALSE 0.999 1.901     5
# BB[4,1]       -1.906 0.902  -3.666  -1.876  -0.275    FALSE 0.990 1.671     6
# BB[5,1]       -1.814 0.751  -3.234  -1.794  -0.597    FALSE 1.000 2.243     5
# BB[1,2]       -0.211 0.357  -0.928  -0.204   0.470     TRUE 0.725 1.000 52500
# BB[2,2]       -0.162 0.339  -0.864  -0.148   0.469     TRUE 0.677 1.000 17195
# BB[3,2]        0.085 0.480  -0.832   0.082   1.068     TRUE 0.571 1.003  1284
# BB[4,2]       -0.093 0.569  -1.225  -0.088   1.007     TRUE 0.566 1.001  4202
# BB[5,2]       -0.146 0.265  -0.680  -0.141   0.354     TRUE 0.707 1.000  7729
# BB[1,3]       13.760 8.860   1.253  13.588  30.372    FALSE 0.996 2.564     4
# BB[2,3]       13.755 8.829   1.330  13.544  30.453    FALSE 0.999 2.558     4
# BB[3,3]       13.746 8.856   1.298  13.512  30.716    FALSE 0.994 2.558     4
# BB[4,3]       13.782 8.917   1.050  13.539  30.560    FALSE 0.992 2.557     4
# BB[5,3]       13.768 8.854   1.391  13.649  30.571    FALSE 0.997 2.555     4
# BB[1,4]       -0.198 0.508  -1.224  -0.187   0.783     TRUE 0.650 1.000  6318
# BB[2,4]        0.389 0.457  -0.453   0.368   1.342     TRUE 0.805 1.001  3404
# BB[3,4]       -0.187 0.621  -1.434  -0.182   1.041     TRUE 0.622 1.007   360
# BB[4,4]       -0.810 0.552  -2.072  -0.756   0.125     TRUE 0.953 1.001  3511
# BB[5,4]        0.163 0.284  -0.404   0.162   0.724     TRUE 0.720 1.001  2552
# BB[1,5]        0.446 0.387  -0.294   0.439   1.235     TRUE 0.879 1.003   819
# BB[2,5]        0.491 0.435  -0.331   0.478   1.371     TRUE 0.876 1.001  5188
# BB[3,5]        0.258 0.527  -0.769   0.252   1.297     TRUE 0.697 1.002  2751
# BB[4,5]        0.457 0.528  -0.524   0.435   1.578     TRUE 0.820 1.001  3178
# BB[5,5]        0.036 0.304  -0.569   0.035   0.634     TRUE 0.547 1.001  3277
# BB[1,6]       -0.628 0.537  -1.797  -0.597   0.340     TRUE 0.895 1.005   450
# BB[2,6]       -0.751 0.363  -1.522  -0.732  -0.096    FALSE 0.989 1.000 16371
# BB[3,6]        0.089 0.601  -0.902   0.027   1.453     TRUE 0.521 1.005   892
# BB[4,6]        0.233 0.640  -0.828   0.165   1.714     TRUE 0.621 1.001  4060
# BB[5,6]       -0.531 0.319  -1.249  -0.497   0.002     TRUE 0.974 1.002 21967
# Sigma.B[1,1]   0.399 0.456   0.094   0.284   1.384    FALSE 1.000 1.001 28962
# Sigma.B[2,1]   0.010 0.302  -0.467   0.004   0.509     TRUE 0.512 1.004 39448
# Sigma.B[3,1]  -0.037 0.432  -0.806  -0.013   0.586     TRUE 0.540 1.001 31552
# Sigma.B[4,1]   0.038 0.420  -0.593   0.017   0.759     TRUE 0.548 1.013 16528
# Sigma.B[5,1]  -0.030 0.308  -0.602  -0.013   0.447     TRUE 0.546 1.002 14472
# Sigma.B[6,1]  -0.025 0.410  -0.723  -0.012   0.599     TRUE 0.533 1.008 52500
# Sigma.B[1,2]   0.010 0.302  -0.467   0.004   0.509     TRUE 0.512 1.004 39448
# Sigma.B[2,2]   0.412 0.470   0.096   0.292   1.433    FALSE 1.000 1.007  6623
# Sigma.B[3,2]  -0.008 0.420  -0.720  -0.002   0.691     TRUE 0.507 1.001  8120
# Sigma.B[4,2]  -0.024 0.403  -0.738  -0.013   0.630     TRUE 0.534 1.004 52500
# Sigma.B[5,2]   0.003 0.315  -0.528   0.002   0.528     TRUE 0.507 1.003 51275
# Sigma.B[6,2]   0.038 0.406  -0.597   0.019   0.760     TRUE 0.553 1.002  6973
# Sigma.B[1,3]  -0.037 0.432  -0.806  -0.013   0.586     TRUE 0.540 1.001 31552
# Sigma.B[2,3]  -0.008 0.420  -0.720  -0.002   0.691     TRUE 0.507 1.001  8120
# Sigma.B[3,3]   0.729 1.230   0.113   0.429   3.166    FALSE 1.000 1.011  3351
# Sigma.B[4,3]  -0.017 0.631  -1.122  -0.006   0.995     TRUE 0.513 1.004  2653
# Sigma.B[5,3]  -0.002 0.468  -0.780   0.001   0.750     TRUE 0.498 1.008 14325
# Sigma.B[6,3]   0.009 0.655  -1.013   0.003   1.073     TRUE 0.506 1.002  2032
# Sigma.B[1,4]   0.038 0.420  -0.593   0.017   0.759     TRUE 0.548 1.013 16528
# Sigma.B[2,4]  -0.024 0.403  -0.738  -0.013   0.630     TRUE 0.534 1.004 52500
# Sigma.B[3,4]  -0.017 0.631  -1.122  -0.006   0.995     TRUE 0.513 1.004  2653
# Sigma.B[4,4]   0.670 0.832   0.123   0.453   2.523    FALSE 1.000 1.006 52500
# Sigma.B[5,4]  -0.067 0.437  -0.890  -0.032   0.573     TRUE 0.580 1.005 36522
# Sigma.B[6,4]  -0.207 0.565  -1.356  -0.110   0.421     TRUE 0.728 1.011  4775
# Sigma.B[1,5]  -0.030 0.308  -0.602  -0.013   0.447     TRUE 0.546 1.002 14472
# Sigma.B[2,5]   0.003 0.315  -0.528   0.002   0.528     TRUE 0.507 1.003 51275
# Sigma.B[3,5]  -0.002 0.468  -0.780   0.001   0.750     TRUE 0.498 1.008 14325
# Sigma.B[4,5]  -0.067 0.437  -0.890  -0.032   0.573     TRUE 0.580 1.005 36522
# Sigma.B[5,5]   0.454 0.515   0.102   0.322   1.603    FALSE 1.000 1.004 52500
# Sigma.B[6,5]  -0.014 0.411  -0.756  -0.006   0.670     TRUE 0.515 1.000 22741
# Sigma.B[1,6]  -0.025 0.410  -0.723  -0.012   0.599     TRUE 0.533 1.008 52500
# Sigma.B[2,6]   0.038 0.406  -0.597   0.019   0.760     TRUE 0.553 1.002  6973
# Sigma.B[3,6]   0.009 0.655  -1.013   0.003   1.073     TRUE 0.506 1.002  2032
# Sigma.B[4,6]  -0.207 0.565  -1.356  -0.110   0.421     TRUE 0.728 1.011  4775
# Sigma.B[5,6]  -0.014 0.411  -0.756  -0.006   0.670     TRUE 0.515 1.000 22741
# Sigma.B[6,6]   0.671 0.929   0.116   0.433   2.637    FALSE 1.000 1.020 52500
# deviance     196.226 5.409 187.210 195.657 208.336    FALSE 1.000 1.001  6290
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
# pD = 14.6 and DIC = 210.85 
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
table(df6$species, df6$condition_h)
df7<-na.omit(df6)


#######################
#Viewing FGB Unknowns 
#######################
species_names <- c("PAST",
                   "OFRA", "CNAT","OFAV","PSTR")

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
  filter(species %in% c("CNAT", "OFAV", "OFRA", "PAST", "PSTR"))

df7_orb <- df7 %>%
  filter(species %in% c("CNAT", "OFAV", "OFRA", "PAST", "PSTR"))

df7_vac <- df7 %>%
  filter(species %in% c("CNAT", "OFAV", "OFRA", "PAST", "PSTR"))

df7

r1 <- ggplot(df7_orb, aes(x = mean_prop_exo, y = state, color = species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 15) + 
  scale_color_manual(values = c("OFRA" = "#bf8f5e",
                                "PAST" = "#bcb554", "CNAT" = "#FFC000", "PSTR" = "#F08000","OFAV"= "#FFD580" )) +
  scale_x_continuous(limits=c(0.0,0.7)) + theme(legend.position = "none") +
  xlab("Proportion Exocytosis") + ylab("Disease State")

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
