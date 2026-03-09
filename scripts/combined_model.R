#######################################
#combined datasets
#Question: Can we use histology to classify unknown FGB samples as
#  Healthy
#  WPH-like disease
#  SCTLD-like disease
#######################################


#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)


WPH <- read_csv("./data/WPH_HistoScore.csv")
SCTLD <- read_csv("./data/Serial_HistoScore.csv")
FGB <- read_csv("./data/FGB_histo_score.csv")

#get FGB to match others 

#rename fungus/sponge
FGB <- FGB %>%
  rename(fungus_sponge = `fungus/sponge`)

#rename health state & gastro/symb columns 
FGB <- FGB %>%
  rename(health_state = `health state`)

FGB <- FGB %>%
  rename( gastro_sep_C = `gastro sep_C`)

FGB <- FGB %>%
  rename( degraded_symb_C = `degraded symb_C`)

FGB <- FGB %>%
  rename( gastro_sep_I = `gastro sep_I`)

FGB <- FGB %>%
  rename( degraded_symb_I = `degraded symb_I`)

ordinal_vars <- c(
  "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C",
  "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I",
  "gastro_sep_I", "degraded_symb_I"
)

binary_vars <- c(
  "fungus_sponge",
  "blown_out_gastro",
  "amoebocytes",
  "loss_of_eosin",
  "loss_of_structure"
)


FGB <- FGB %>%
  mutate(across(
    all_of(ordinal_vars),
    ~ case_when(
      . == "L" ~ 1,
      . == "M" ~ 2,
      . == "H" ~ 3,
      TRUE ~ as.numeric(.)
    )
  )) %>%
  mutate(across(
    all_of(binary_vars),
    ~ case_when(
      . == "Y" ~ 1,
      . == "N" ~ 0,
      TRUE ~ as.numeric(.)
    )
  ))


#rename health state & gastro/symb columns 
WPH <- WPH %>%
  rename(health_state = `Sample_Health_State`)

#rename health state & gastro/symb columns 
SCTLD <- SCTLD %>%
  rename(health_state = `worded_health`)

#create disease label columns
FGB$dataset <- "FGB"

FGB$class <- ifelse(FGB$health_state == "HH",
                    "Healthy",
                    "FGB_Diseased")

WPH$dataset <- "WPH"

WPH$class <- ifelse(WPH$health_state == "HH",
                    "Healthy",
                    "WPH")

SCTLD$dataset <- "SCTLD"

SCTLD$class <- ifelse(SCTLD$health_state == "HH",
                    "Healthy",
                    "SCTLD")

#combine all datasets
combined <- rbind(FGB, WPH, SCTLD)

#####################################################################
#First we need to combine WPH & SCTLD to build training dataset
#####################################################################

