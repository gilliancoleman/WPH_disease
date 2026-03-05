#################################
#Let's build some decision trees 
#FGB
#################################

#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)

#which disease signs matter?
FGB <- read_csv("./data/FGB_histo_score.csv")

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
 
# Convert categorical variables to factors
FGB$health_state <- as.factor(FGB$health_state)  
FGB$Individual <- as.factor(FGB$Individual)  

#check out some tables
table(FGB$species, FGB$health_state)
table(FGB$health_state, FGB$sample_date)


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

#run rf and xgb just to see where everything's at before combining all diseased and healthy

# Select relevant columns for analysis
FGB_subset <- FGB[, c("health_state",
                        "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                        "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                        "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                        "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure")]

# Remove rows with missing values
FGB_subset_clean <- na.omit(FGB_subset) #wooo didn't lose anybody


#Train Random Forest model

set.seed(123)
rf_model <- randomForest(health_state ~ ., data =FGB_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Call:
#   randomForest(formula = health_state ~ ., data = FGB_subset_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 54.04%
# Confusion matrix:
#   DD DH DL HD HH class.error
# DD 21  0  0 11 20   0.5961538
# DH  0  0  0  0  1   1.0000000
# DL  0  0  0  0  1   1.0000000
# HD 14  0  0 19 20   0.6415094
# HH 11  0  0  9 34   0.3703704

#probably going to need to combine disease states because DH/DL only have n=1

# Plot feature importance
varImpPlot(rf_model)


#now let's do the XGB model
predictor_vars <- c(
  "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C",
  "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I",
  "gastro_sep_I", "degraded_symb_I",
  "fungus_sponge", "blown_out_gastro", "amoebocytes",
  "loss_of_eosin", "loss_of_structure"
)

#build x & y
X <- as.matrix(FGB_subset_clean[, predictor_vars])

y <- as.integer(FGB_subset_clean$health_state) - 1

#train xgb boost

dtrain <- xgb.DMatrix(data = X, label = y)

params <- list(
  objective = "multi:softprob",
  num_class = length(levels(FGB_subset_clean$health_state)),
  max_depth = 3,
  eta = 0.1,
  eval_metric = "mlogloss"
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#see what's important 
imp <- xgb.importance(model = xgb_model)
imp

# Feature       Gain      Cover  Frequency
# <char>      <num>      <num>      <num>
#   1:     loss_of_eosin 0.10538391 0.05421004 0.08349783
# 2:      exocytosis_C 0.09708431 0.07377848 0.07639163
# 3:   degraded_symb_C 0.09618468 0.06047375 0.08941966
# 4:      gastro_sep_C 0.08757712 0.11601915 0.10501382
# 5:   vacuolization_I 0.08383315 0.08067672 0.07718121
# 6:   vacuolization_C 0.07430262 0.06542039 0.07165417
# 7:       amoebocytes 0.06993038 0.05350555 0.03750493
# 8:        necrosis_C 0.06190699 0.09450964 0.09139360
# 9:      exocytosis_I 0.06123850 0.06763285 0.08369522
# 10:   degraded_symb_I 0.05872982 0.05488539 0.05073036
# 11: loss_of_structure 0.05355623 0.05357498 0.05013818
# 12:  blown_out_gastro 0.05301939 0.02833257 0.03414923
# 13:      gastro_sep_I 0.04809191 0.10124167 0.06928543
# 14:        necrosis_I 0.03270358 0.05838962 0.06040268
# 15:     fungus_sponge 0.01645742 0.03734921 0.01954205

#compare training accuracy

pred <- predict(xgb_model, dtrain)

pred_matrix <- matrix(pred,
                      ncol = length(levels(FGB_subset_clean$health_state)),
                      byrow = TRUE)

pred_class <- max.col(pred_matrix) - 1

mean(pred_class == y)
# 0.1863354 #super low


#now let's model HH vs All Disease States

#create a binary health variable 
FGB$health_binary <- ifelse(
  FGB$health_state %in% c("HD", "DD", "DL", "DH"),
  "Diseased",
  "Healthy"
)

FGB$health_binary <- as.factor(FGB$health_binary)

#check it
table(FGB$health_binary) #twice as much diseased samples now instead of health but we'll see 

#re subset for this one 

FGB_subset_bin <- FGB[, c(
  "health_binary",
  "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C",
  "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I",
  "gastro_sep_I", "degraded_symb_I",
  "fungus_sponge", "blown_out_gastro", "amoebocytes",
  "loss_of_eosin", "loss_of_structure"
)]

FGB_subset_bin_clean <- na.omit(FGB_subset_bin) #kept all the folks

#random forest model

set.seed(123)

rf_binary <- randomForest(
  health_binary ~ .,
  data = FGB_subset_bin_clean,
  importance = TRUE,
  ntree = 500
)

print(rf_binary)

# Call:
#   randomForest(formula = health_binary ~ ., data = FGB_subset_bin_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 32.92%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       90      17   0.1588785
# Healthy        36      18   0.6666667

varImpPlot(rf_binary)


#build xgb model



#build matrix 
predictor_vars <- colnames(FGB_subset_bin_clean)
predictor_vars <- predictor_vars[predictor_vars != "health_binary"]

X <- as.matrix(FGB_subset_bin_clean[, predictor_vars])

y <- ifelse(FGB_subset_bin_clean$health_binary == "Diseased", 1, 0)

#train model

dtrain <- xgb.DMatrix(data = X, label = y)

params_bin <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 3,
  eta = 0.1
)

xgb_binary <- xgb.train(
  params = params_bin,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)


#check accuracy
pred_prob <- predict(xgb_binary, dtrain)

pred_class <- ifelse(pred_prob > 0.5, 1, 0)

mean(pred_class == y) #0.8944099  #accuracy jumped wayyyyyyy up 


#importance
imp_bin <- xgb.importance(model = xgb_binary)
imp_bin

# Feature       Gain      Cover   Frequency
# <char>      <num>      <num>       <num>
#   1: loss_of_structure 0.12919342 0.06291550 0.074289876
# 2:     loss_of_eosin 0.10890407 0.07217573 0.108521486
# 3:   vacuolization_I 0.09054739 0.13410833 0.104151493
# 4:   degraded_symb_C 0.08727499 0.07173579 0.060451566
# 5:      gastro_sep_C 0.08614500 0.13341393 0.101966497
# 6:       amoebocytes 0.07560784 0.08527326 0.057538237
# 7:  blown_out_gastro 0.07248395 0.02875659 0.038601602
# 8:      exocytosis_I 0.07208925 0.05184441 0.083029862
# 9:   degraded_symb_I 0.05828456 0.03500056 0.053896577
# 10:   vacuolization_C 0.05437933 0.06379012 0.075018208
# 11:        necrosis_I 0.04504637 0.05866473 0.071376548
# 12:        necrosis_C 0.04380125 0.08401543 0.070648216
# 13:      gastro_sep_I 0.03660265 0.05452375 0.045156591
# 14:      exocytosis_C 0.02153128 0.04545861 0.045884924
# 15:     fungus_sponge 0.01810864 0.01832325 0.009468318

###################################################################################
#from what I've read....the rf stuff is good at downweighing less important variables, so maybe leave everything and it looks like there's a couple ways we can test the stability of the model without dropping variables....it also wasn't recommended to drop variables if n is small unless it's just for a publication pic
###################################################################################


#so let's increase ntree to test stability of rf model

rf_model <- randomForest(
  health_binary ~ .,
  data = FGB_subset_bin_clean,
  importance = TRUE,
  ntree = 1000
)

print(rf_model)

# Call:
#   randomForest(formula = health_binary ~ ., data = FGB_subset_bin_clean,      importance = TRUE, ntree = 1000) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 31.06%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       95      12   0.1121495
# Healthy        38      16   0.7037037




#now let's cross validate xgb model using xgb.cv(): Splits data into 5 folds, Trains on 4, and then test on the other 
set.seed(123)

cv <- xgb.cv(
  params = list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 3,
    eta = 0.1
  ),
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 0,
  early_stopping_rounds = 20
)

#see what the best round was 
best_nrounds <- which.min(cv$evaluation_log$test_logloss_mean)
best_nrounds


#now train the model off of the best rounds 
xgb_model_best <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 3,
    eta = 0.1
  ),
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

#check accuracy

cv_acc <- xgb.cv(
  params = list(
    objective = "binary:logistic",
    eval_metric = "error",
    max_depth = 3,
    eta = 0.1
  ),
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 0,
  early_stopping_rounds = 20
)

head(cv_acc$evaluation_log)

#accuracy is 1-0.336 = 0.664 so 66% which is ~ equal to RF OOB (~69%) --> reassuring 

# RF OOB ≈ 69%
# 
# XGB 5-fold CV ≈ 66%
# 
# Two completely different algorithms converging on the same performance.
# 
# That strongly suggests:
#   
# Your real predictive ability is ~65–70%

