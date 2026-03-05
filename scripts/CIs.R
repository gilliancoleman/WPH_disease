##################
#Comparing CIs 
##################

library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)

#load em up
WPH_FGB80 <- read.csv("./data/CIs/WPvsFGB80.csv")
WPH_FGB95 <- read.csv("./data/CIs/WPvsFGB95.csv")
SCTLD_FGB95 <- read.csv("./data/CIs/SCTLDvsFGB95.csv")
SCTLD_FGB80 <- read.csv("./data/CIs/SCTLDvsFGB80.csv")


#function to see if CI crosses 0
crosses_zero <- function(lower, upper) {
  lower <= 0 & upper >= 0
}

#do 95% CI first

wp_vs_fgb95 <- WPH_FGB95 %>%
  mutate(crosses0 = crosses_zero(X2.50., X97.5.))

sctld_vs_fgb95 <- SCTLD_FGB95 %>%
  mutate(crosses0 = crosses_zero(X2.50., X97.5.))

#merge the two comparisons
compare_fgb95 <- wp_vs_fgb95 %>%
  dplyr::select(species, parameter, crosses0) %>%
  rename(crosses0_wp = crosses0) %>%
  left_join(
    sctld_vs_fgb95 %>%
      dplyr::select(species, parameter, crosses0) %>%
      rename(crosses0_sctld = crosses0),
    by = c("species", "parameter")
  )

#which is it closer to?
compare_fgb95 <- compare_fgb95 %>%
  mutate(closer_to = case_when(
    crosses0_wp & !crosses0_sctld ~ "WP",
    !crosses0_wp & crosses0_sctld ~ "SCTLD",
    crosses0_wp & crosses0_sctld  ~ "Both",
    TRUE                           ~ "Distinct from both"
  ))

#drummmmm rollll
compare_fgb95 %>% arrange(species, parameter)


#Compute distance between posterior means
wp_vs_fgb95 <- wp_vs_fgb95 %>%
  mutate(distance = abs(mean))

sctld_vs_fgb95 <- sctld_vs_fgb95 %>%
  mutate(distance = abs(mean))

#reorder my parameters

compare_fgb95 <- compare_fgb95 %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("intercept", "prop_exo", "gastro", "symb_vac", "avg_symb", "degr_symb")
    )
  )


#visualize it 

ggplot(compare_fgb95,
       aes(x = parameter, y = species, fill = closer_to)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "WP" = "#FFB6C1",
      "SCTLD" = "#FF00FF",
      "Both" = "#A95C68",
      "Distinct from both" = "#C9A9A6"
    )
  ) +
  labs(
    x = "Parameter",
    y = "Species",
    fill = "FGB is closer to"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


#visualize posterior means
distance_95 <- wp_vs_fgb95 %>%
  dplyr::select(species, parameter, distance) %>%
  mutate(comparison = "WP") %>%
  bind_rows(
    sctld_vs_fgb95 %>%
      dplyr::select(species, parameter, distance) %>%
      mutate(comparison = "SCTLD")
  )

ggplot(distance_95,
       aes(x = parameter, y = distance, fill = comparison)) +
  geom_col(position = "dodge") +
  facet_wrap(~ species) +
  theme_classic(base_size = 14) +
  labs(
    x = "Parameter",
    y = "Effect size difference from FGB",
    fill = "Compared to"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################################
#go ahead and do 80% while we're at it
########################################

wp_vs_fgb80 <- WPH_FGB80 %>%
  mutate(crosses0 = crosses_zero(X10.00., X90.0.))

sctld_vs_fgb80 <- SCTLD_FGB80 %>%
  mutate(crosses0 = crosses_zero(X10.00., X90.0.))

#merge the two comparisons
compare_fgb80 <- wp_vs_fgb80 %>%
  dplyr::select(species, parameter, crosses0) %>%
  rename(crosses0_wp = crosses0) %>%
  left_join(
    sctld_vs_fgb80 %>%
      dplyr::select(species, parameter, crosses0) %>%
      rename(crosses0_sctld = crosses0),
    by = c("species", "parameter")
  )

#which is it closer to?
compare_fgb80 <- compare_fgb80 %>%
  mutate(closer_to = case_when(
    crosses0_wp & !crosses0_sctld ~ "WP",
    !crosses0_wp & crosses0_sctld ~ "SCTLD",
    crosses0_wp & crosses0_sctld  ~ "Both",
    TRUE                           ~ "Distinct from both"
  ))

#drummmmm rollll
compare_fgb80 %>% arrange(species, parameter)

#rename mean column so next section works
wp_vs_fgb80 <- wp_vs_fgb80 %>% rename(post_mean = mean)
sctld_vs_fgb80 <- sctld_vs_fgb80 %>% rename(post_mean = mean)

#Compute distance between posterior means
wp_vs_fgb80 <- wp_vs_fgb80 %>%
  mutate(distance = abs(post_mean))

sctld_vs_fgb80 <- sctld_vs_fgb80 %>%
  mutate(distance = abs(post_mean))


#reorder my parameters

compare_fgb80 <- compare_fgb80 %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("intercept", "prop_exo", "gastro", "symb_vac", "avg_symb", "degr_symb")
    )
  )


#visualize it 

ggplot(compare_fgb80,
       aes(x = parameter, y = species, fill = closer_to)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "WP" = "#FFB6C1",
      "SCTLD" = "#FF00FF",
      "Both" = "#A95C68",
      "Distinct from both" = "#C9A9A6"
    )
  ) +
  labs(
    x = "Parameter",
    y = "Species",
    fill = "FGB is closer to"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )




#visualize posterior means
distance_80 <- wp_vs_fgb80 %>%
  dplyr::select(species, parameter, distance) %>%
  mutate(comparison = "WP") %>%
  bind_rows(
    sctld_vs_fgb80 %>%
      dplyr::select(species, parameter, distance) %>%
      mutate(comparison = "SCTLD")
  )

ggplot(distance_80,
       aes(x = parameter, y = distance, fill = comparison)) +
  geom_col(position = "dodge") +
  facet_wrap(~ species) +
  theme_classic(base_size = 14) +
  labs(
    x = "Parameter",
    y = "Effect size difference from FGB",
    fill = "Compared to"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#osterior means are on the log-odds (logit) scale

#that means:

# 0 = no effect on log-odds
# 
# ±1 = moderate effect
# 
# ±5 = enormous effect
# 
# ±10 = essentially deterministic (probability ≈ 0 or 1)

#we can transform the log-odds value
plogis(post_mean)
