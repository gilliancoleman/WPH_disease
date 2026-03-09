######################
#Trees Combined 
#03/08/2026
#####################

#You want to train a model that can classify a sample as:
# Healthy
# WPH-type disease
# SCTLD-type disease

#And then apply that trained model to FGB samples to see where they fall.

#Does an FGB diseased sample look more like WPH or SCTLD histologically?

#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)
library(caret)
library(data.table)
library(Ckmeans.1d.dp)

############################################################
#First we need to combine :WPH dataset + SCTLD dataset
############################################################

WPH <- read_csv("./data/WPH_HistoScore.csv")
SCTLD <- read_csv("./data/Serial_HistoScore.csv")
FGB <- read_csv("./data/FGB_histo_score.csv")

str(FGB)
str(WPH)
str(SCTLD)

#standardize the health column
WPH$health_state <- WPH$Sample_Health_State
SCTLD$health_state <- SCTLD$worded_health
FGB$health_state <- FGB$'health state'

#rename columns to match
FGB <- FGB %>%
  rename(fungus_sponge = `fungus/sponge`)


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

#redo LMH to be 123

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

#rename SCTLD columns to match WPH: 
SCTLD <- SCTLD %>%
  rename( gastro_sep_C = `gastro_sep_c`)

SCTLD <- SCTLD %>%
  rename( degraded_symb_C = `degraded_symb_c`)

SCTLD <- SCTLD %>%
  rename( gastro_sep_I = `gastro_sep_i`)

SCTLD <- SCTLD %>%
  rename( degraded_symb_I = `degraded_symb_i`)

SCTLD <- SCTLD %>%
  rename( necrosis_C = `necrosis_c`)

SCTLD <- SCTLD %>%
  rename( vacuolization_C = `vacuolization_c`)

SCTLD <- SCTLD %>%
  rename( exocytosis_C = `exocytosis_c`)

SCTLD <- SCTLD %>%
  rename( necrosis_I = `necrosis_i`)

SCTLD <- SCTLD %>%
  rename( vacuolization_I = `vacuolization_i`)

SCTLD <- SCTLD %>%
  rename( exocytosis_I = `exocytosis_i`)
####################################################################################################
#First we need to combine WPH & SCTLD to build training dataset *FGB* will not be included in training 
####################################################################################################


#theres some different columns but not really ones that matter just extra so
common_cols <- intersect(colnames(WPH), colnames(SCTLD))
common_cols

#subset both based on common columns
WPH_clean <- WPH[, common_cols]
SCTLD_clean <- SCTLD[, common_cols]

#now combine
training_data <- rbind(WPH_clean, SCTLD_clean)

training_data$class <- as.factor(training_data$class)

#check it out
dim(training_data)
table(training_data$class)

#assign predictor variables; there's more than not so we'll just say which not to include 
predictor_vars <- training_data |>
  select(-class, -health_state, -Individual, -health_state, -Site, -Initials, -Sample_Health_State, -date_completed, -Sample_Date, -Pics_taken, -score, -measurements, -Notes, -DONE, -dataset) |>
  colnames()

#lets check/ Should be 15
predictor_vars #checks out 


#now for the model:

set.seed(123)

rf_model <- randomForest(
  x = training_data[, predictor_vars],
  y = training_data$class,
  ntree = 1000,
  importance = TRUE
)

print(rf_model)
# Call:
#   randomForest(x = training_data[, predictor_vars], y = training_data$class,      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 51.22%
# Confusion matrix:
#   Healthy SCTLD WPH class.error
# Healthy      35     4   5   0.2045455
# SCTLD        13     5   2   0.7500000
# WPH          14     4   0   1.0000000

##############################################################
#lets break it down now because I know I'll forget if I don't 
# OOB error = 51.22%
# Accuracy ≈ 48.8%

#This is a 3 class model so random guessing would be ~33.33%
# so the model is predicting but the signal is weak 

#the model does recognize healthy histo really well! : Out of the Health ones 35 were predicted correctly so it has 80% accuracy (100%-20% error)

#most SCTLD samples are predicted as healthy (only 5/20 correct)

#model never gets WPH right 

#could either mean that WPH overlaps with Healthy or SCTLD a lot or that sample size is too small (it is kind of small and that could also explain why SCTLD accuracy is lower too)

#Earlier results hinted at this: WPH binary model → moderate signal & SCTLD binary model → weak signal and now we're trying to add a class with only ~80 total samples with a huge class bias towards healthy

##############################################################

#look at which variables are important 
varImpPlot(rf_model)

##############################################################
#on the brightside the first 5 predictors are the same on these plots

#If the same predictors appear in both, it means those variables consistently carry signal. Those are likely our key histological indicators.
#############################################################

#before we just pull those 5 predictors out, since we have such a large class imbalance lets run a. model with WPH & SCTLD combined and see if it can pick between healthy and diseased --> and then see if the two diseases can be separated 

#create 2 step class columns

# Keep only common columns & make combined & not training dataset

#add a species column to SCTLD 

SCTLD$species <- "OFRA"

# Numeric mapping for histology scores
score_cols <- c("necrosis_C","vacuolization_C","exocytosis_C","gastro_sep_C",
                "degraded_symb_C","necrosis_I","vacuolization_I","exocytosis_I",
                "gastro_sep_I","degraded_symb_I")

# Binary mapping for presence/absence columns
bin_cols <- c("fungus_sponge", "blown_out_gastro", "amoebocytes", 
              "loss_of_eosin", "loss_of_structure")


common_cols <- c("health_state", score_cols, bin_cols, "species")
combined <- bind_rows(
  WPH[, intersect(common_cols, names(WPH))],
  SCTLD[, intersect(common_cols, names(SCTLD))],
  FGB[, intersect(common_cols, names(FGB))]
)

# Remove rows with NA
combined <- na.omit(combined)

# Step 1: Binary classification (Healthy vs Diseased)
combined <- combined %>%
  mutate(class_bin = ifelse(health_state %in% c("HH"), "Healthy", "Diseased") %>% factor())

# Step 2: Diseased subclass (only for diseased samples)
diseased_subset <- combined %>%
  filter(class_bin == "Diseased") %>%
  mutate(class_disease = species %>% factor())  # WPH vs SCTLD vs FGB diseased

#define predictor variables
predictor_vars <- c(score_cols, bin_cols)


#now run bunary rf model (healthy vs diseased)

set.seed(123)
rf_bin <- randomForest(
  x = combined[, predictor_vars],
  y = combined$class_bin,
  ntree = 1000,
  importance = TRUE
)

print(rf_bin)
varImpPlot(rf_bin)

# Call:
#   randomForest(x = combined[, predictor_vars], y = combined$class_bin,      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 31.68%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       92      15   0.1401869
# Healthy        36      18   0.6666667

######################
#interpret
#model misclassifies ~32% (OOB = 32%)
#about 86% of disease is predicted correctly 
#about 33% correct for healthy

#The model is biased toward predicting Diseased because there are more diseased samples (92+15=107 vs 36+18=54).

######################

#let's see if it's the same variables driving these differences as last time 
varImpPlot(rf_bin) #spoiler...it's not 

#the next option: You could improve accuracy for Healthy by: Using class weights in Random Forest or XGBoost or Possibly undersampling Diseased or oversampling Healthy during training.

#let's try weighting classes 

#XGB model accounting for the class imbalance observed in RF.

#make pred vars
predictor_vars <- c(
  "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C",
  "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I",
  "gastro_sep_I", "degraded_symb_I",
  "fungus_sponge", "blown_out_gastro", "amoebocytes",
  "loss_of_eosin", "loss_of_structure"
)

#Make sure combined is a data.table
combined <- as.data.table(combined)

#build x matrix and label vector
X_bin <- as.matrix(combined[, ..predictor_vars])
y_bin <- as.integer(combined$class_bin) - 1  # 0 = Healthy, 1 = Diseased

#create d matrix
dtrain_bin <- xgb.DMatrix(data = X, label = y)

#Handle class imbalance
# Calculate scale_pos_weight = #Healthy / #Diseased
num_healthy <- sum(y == 0)
num_diseased <- sum(y == 1)
scale_pos_weight <- num_healthy / num_diseased

#Define parameters
params_bin <- list(
  objective = "binary:logistic",
  eval_metric = "error",
  max_depth = 3,
  eta = 0.1,
  scale_pos_weight = scale_pos_weight
)


# Run cross-validation to choose nrounds
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  maximize = FALSE
)

# Get the best number of rounds safely
best_round <- cv_results$best_iteration
if (is.null(best_round) || length(best_round) == 0) {
  best_round <- 100  # fallback if xgb.cv didn't find a best iteration
}

# Train final model using best_round
xgb_bin <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_round,
  verbose = 0
)

# Feature importance
imp_bin <- xgb.importance(feature_names = predictor_vars, model = xgb_bin)
print(imp_bin)

# Feature       Gain      Cover  Frequency
# <char>      <num>      <num>      <num>
#   1: loss_of_structure 0.11394933 0.05932323 0.06914894
# 2:   vacuolization_I 0.11301053 0.11905036 0.08865248
# 3:      gastro_sep_C 0.10429409 0.13082603 0.11347518
# 4:     loss_of_eosin 0.09451132 0.07697404 0.11879433
# 5:       amoebocytes 0.08898916 0.10047277 0.06205674
# 6:   degraded_symb_C 0.08312436 0.06761663 0.04787234
# 7:   degraded_symb_I 0.07334368 0.03873052 0.07446809
# 8:        necrosis_C 0.06655113 0.09073689 0.07092199
# 9:   vacuolization_C 0.06196558 0.04377976 0.07092199
# 10:      exocytosis_I 0.05465647 0.06513089 0.11347518
# 11:  blown_out_gastro 0.05107706 0.03181488 0.03723404
# 12:     fungus_sponge 0.03000836 0.06923258 0.03723404
# 13:      gastro_sep_I 0.02824084 0.05617114 0.04432624
# 14:        necrosis_I 0.02037157 0.01724840 0.03014184
# 15:      exocytosis_C 0.01590652 0.03289189 0.02127660

################################
#Interpretation:
#loss_of_structure, vacuolization_I, gastro_sep_C, and loss_of_eosin are the strongest predictors for distinguishing Healthy vs Diseased in combined dataset.

#This roughly matches what the Random Forest OOB suggested: the top histology predictors are consistent.
################################
xgb.plot.importance(imp_bin)


#now do the rf multipstep (disease type only)
set.seed(123)
rf_disease <- randomForest(
  x = diseased_subset[, predictor_vars],
  y = diseased_subset$class_disease,
  ntree = 1000,
  importance = TRUE
)

print(rf_disease)
#Call:
# randomForest(x = diseased_subset[, predictor_vars], y = diseased_subset$class_disease,      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 47.66%
# Confusion matrix:
#   CNAT OFAV OFRA PAST PSTR class.error
# CNAT    7    3    0    3    9   0.6818182
# OFAV    5   12    2    0    3   0.4545455
# OFRA    0    7    2    0    2   0.8181818
# PAST    3    0    0    6    1   0.4000000
# PSTR    5    6    0    2   29   0.3095238

#############################
#Interpretation:
#OOB error is 48% which is kind of high and means the model is only correct close to half the time & is having a hard time separating the disease types

#OFRA has the highest error (~82%) — the model struggles to correctly classify OFRA diseased samples.

#PSTR is the easiest to classify (~31% error). 

#############################

varImpPlot(rf_disease)

#now XGB multiclass

X_d <- as.matrix(diseased_subset[, predictor_vars])
y_d <- as.integer(diseased_subset$class_disease) - 1

dtrain_d <- xgb.DMatrix(data = X_d, label = y_d)

params_d <- list(
  objective = "multi:softprob",
  num_class = length(levels(diseased_subset$class_disease)),
  eval_metric = "mlogloss",
  max_depth = 3,
  eta = 0.1
)

cv_results_d <- xgb.cv(
  params = params_d,
  data = dtrain_d,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  maximize = FALSE
)

best_round_d <- cv_results_d$best_iteration
if (is.null(best_round_d) || length(best_round_d) == 0) best_round_d <- 100

xgb_d <- xgb.train(
  params = params_d,
  data = dtrain_d,
  nrounds = best_round,
  verbose = 1
)


imp_d <- xgb.importance(feature_names = predictor_vars, model = xgb_d)
print(imp_d)

#            Feature       Gain       Cover  Frequency
# <char>      <num>       <num>      <num>
#   1:   vacuolization_C 0.21901938 0.196944674 0.14250513
# 2:        necrosis_C 0.14313239 0.091479434 0.07761807
# 3:      gastro_sep_C 0.12425191 0.069555312 0.11129363
# 4:        necrosis_I 0.10026344 0.109431477 0.13223819
# 5:   degraded_symb_C 0.07667869 0.084602020 0.07022587
# 6:      exocytosis_C 0.06616761 0.047786332 0.07843943
# 7:      gastro_sep_I 0.05445591 0.080949834 0.06570842
# 8:   vacuolization_I 0.04660307 0.070082428 0.07433265
# 9:       amoebocytes 0.04515793 0.062843166 0.04928131
# 10: loss_of_structure 0.04104962 0.054955833 0.05420945
# 11:   degraded_symb_I 0.02865992 0.050745197 0.06242300
# 12:      exocytosis_I 0.01950295 0.028358556 0.04024641
# 13:     loss_of_eosin 0.01879559 0.043814344 0.02505133
# 14:  blown_out_gastro 0.01626161 0.008451392 0.01642710

#############################
#Interpretation:
#vacuolization_C and necrosis_C are the strongest predictors for separating the diseased types.
#Most of the top 5 predictors are consistent with the binary Healthy/Diseased XGB. This is expected because the same histology features are informative for both binary and multiclass distinctions.
#############################
xgb.plot.importance(imp_d)

##################################
#BIG PIC
#The XGBoost and Random Forest models are highlighting similar histology features, which reinforces that these are key disease markers.

#Binary model (Healthy vs Diseased):
  #OOB error ~32% (from previous step)
  #Most important predictors: loss_of_structure, vacuolization_I, gastro_sep_C, loss_of_eosin

#Diseased-only model (WPH vs SCTLD):
  #OOB error ~48%
  #Top predictors: vacuolization_C, necrosis_C, gastro_sep_C, necrosis_I
  #Classes with fewer samples (OFRA) are hardest to predict

#Feature consistency:
  #Both RF and XGBoost agree on the top histology predictors.
  #Variables with low gain/frequency could be candidates to drop in future model simplifications.

#Class imbalance:
  # handled it in XGB with scale_pos_weight for binary. For multiclass, might consider sampling techniques or weight adjustments if we want to improve accuracy for rarer disease classes.

##################################

###################################################################################
#I think maybe now we should try keeping the top 5 predictors and running a model
#this will potentially increase accuracy if low-importance features were causing overfitting
###################################################################################

#we're going to do this binary again and look at healthy vs diseased first

#for binary XGB (healthy vs diseased)


#Define top 5 predictors

top5_vars <- c("loss_of_structure", "vacuolization_I", "gastro_sep_C", 
               "loss_of_eosin", "amoebocytes")


#Binary model: Healthy vs Diseased

# Prepare data
X_bin <- as.matrix(combined[, ..top5_vars])
y_bin <- as.integer(combined$class_bin) - 1  # 0 = Healthy, 1 = Diseased
dtrain_bin <- xgb.DMatrix(data = X_bin, label = y_bin)

# Handle class imbalance
num_healthy <- sum(y_bin == 0)
num_diseased <- sum(y_bin == 1)
scale_pos_weight <- num_healthy / num_diseased

# XGBoost parameters
params_bin <- list(
  objective = "binary:logistic",
  eval_metric = "error",
  max_depth = 3,
  eta = 0.1,
  scale_pos_weight = scale_pos_weight
)

# Cross-validation for nrounds
cv_bin <- xgb.cv(
  params = params_bin,
  data = dtrain_bin,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  maximize = FALSE
)

best_round_bin <- cv_bin$best_iteration
if (is.null(best_round_bin) || length(best_round_bin) == 0) best_round_bin <- 100

# Train final binary XGBoost model
xgb_bin <- xgb.train(
  params = params_bin,
  data = dtrain_bin,
  nrounds = best_round_bin,
  verbose = 1
)

# Feature importance
imp_bin <- xgb.importance(feature_names = top5_vars, model = xgb_bin)

print(imp_bin)
#             Feature      Gain     Cover Frequency
#<char>     <num>     <num>     <num>
#  1: loss_of_structure 0.2856236 0.1590569 0.2030568
#2:     loss_of_eosin 0.2315248 0.1321076 0.2314410
#3:      gastro_sep_C 0.1998020 0.2075365 0.1943231
#4:       amoebocytes 0.1615393 0.2125325 0.1462882
#5:   vacuolization_I 0.1215103 0.2887666 0.2248908

xgb.plot.importance(imp_bin)

# Random Forest for binary
rf_bin <- randomForest(
  x = combined[, ..top5_vars],
  y = combined$class_bin,
  ntree = 1000,
  importance = TRUE
)

print(rf_bin)
# Call:
#   randomForest(x = combined[, ..top5_vars], y = combined$class_bin,      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 34.16%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       90      17   0.1588785
# Healthy        38      16   0.7037037

varImpPlot(rf_bin)



#Diseased-only model: classify WPH vs SCTLD

# Subset only diseased samples
diseased_subset <- combined[class_bin == "Diseased"]

# XGBoost multiclass
X_d <- as.matrix(diseased_subset[, ..top5_vars])
y_d <- as.integer(diseased_subset$class_disease) - 1
dtrain_d <- xgb.DMatrix(data = X_d, label = y_d)

`params_d <- list(
  objective = "multi:softprob",
  num_class = length(levels(diseased_subset$class_disease)),
  eval_metric = "mlogloss",
  max_depth = 3,
  eta = 0.1
)

# CV for best nrounds
cv_d <- xgb.cv(
  params = params_d,
  data = dtrain_d,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  maximize = FALSE
)

best_round_d <- cv_d$best_iteration
if (is.null(best_round_d) || length(best_round_d) == 0) best_round_d <- 100

# Train final XGBoost diseased model
xgb_d <- xgb.train(
  params = params_d,
  data = dtrain_d,
  nrounds = best_round_d,
  verbose = 1
)

# Feature importance for diseased
imp_d <- xgb.importance(feature_names = top5_vars, model = xgb_d)
print(imp_d)
xgb.plot.importance(imp_d)

# Random Forest for diseased-only
rf_disease <- randomForest(
  x = diseased_subset[, ..top5_vars],
  y = diseased_subset$class_disease,
  ntree = 1000,
  importance = TRUE
)

print(rf_disease)
varImpPlot(rf_disease)