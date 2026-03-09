#################################
#Let's build some decision trees 
#SCTLD
#################################

#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)

#which disease signs matter?
SCTLD <- read_csv("./data/Serial_HistoScore.csv")

# Convert categorical variables to factors
SCTLD$worded_health <- as.factor(SCTLD$worded_health)  
SCTLD$Individual <- as.factor(SCTLD$Individual)  

#add species column
SCTLD$species <- "OFRA"

table(SCTLD$species, SCTLD$worded_health)
table(SCTLD$worded_health, SCTLD$Sample_Date)

#rename fungus/sponge
SCTLD <- SCTLD %>%
  rename(fungus_sponge = `fungus/sponge`)


# Select relevant columns for analysis #not using Sample_Health State
SCTLD_subset <- SCTLD[, c("worded_health", 
                        "necrosis_c", "vacuolization_c", "exocytosis_c", "gastro_sep_c", 
                        "degraded_symb_c", "necrosis_i", "vacuolization_i", "exocytosis_i", 
                        "gastro_sep_i", "degraded_symb_i", "fungus_sponge",
                        "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure")]


# Remove rows with missing values
SCTLD_subset_clean <- na.omit(SCTLD_subset) #wooo didn't lose anybody


#Train Random Forest model

set.seed(123)
rf_model <- randomForest(worded_health ~ ., data =SCTLD_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Call:
#   randomForest(formula = worded_health ~ ., data = SCTLD_subset_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 46.94%
# Confusion matrix:
#   DD HD HH class.error
# DD  1  2  7    0.900000
# HD  1  0  9    1.000000
# HH  1  3 25    0.137931

####The model cannot separate HD at all. --> so maybe we go ahead and combine HD & DD and only have two classes before we keep going but let's see what the XGB model does first

# Plot feature importance
varImpPlot(rf_model)

#now let's do the XGB model
predictor_vars <- c(
  "necrosis_c", "vacuolization_c", "exocytosis_c", "gastro_sep_c",
  "degraded_symb_c", "necrosis_i", "vacuolization_i", "exocytosis_i",
  "gastro_sep_i", "degraded_symb_i",
  "fungus_sponge", "blown_out_gastro", "amoebocytes",
  "loss_of_eosin", "loss_of_structure"
)

#build x & y
X <- as.matrix(SCTLD_subset_clean[, predictor_vars])

y <- as.integer(SCTLD_subset_clean$worded_health) - 1

#train xgb boost

dtrain <- xgb.DMatrix(data = X, label = y)

params <- list(
  objective = "multi:softprob",
  num_class = length(levels(SCTLD_subset_clean$worded_health)),
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
#   1:   degraded_symb_c 0.15937906 0.16280591 0.12375720
# 2:      gastro_sep_i 0.14003121 0.08689177 0.10125589
# 3:   vacuolization_i 0.13913622 0.12871695 0.12401884
# 4:      exocytosis_c 0.11771828 0.10436534 0.10125589
# 5:        necrosis_c 0.10232639 0.11621006 0.08686552
# 6:      gastro_sep_c 0.09737083 0.04789189 0.07456829
# 7:      exocytosis_i 0.07026675 0.10457345 0.12244898
# 8:     loss_of_eosin 0.05967752 0.09086875 0.08451073
# 9:   degraded_symb_i 0.05465882 0.08481523 0.08765044
# 10:   vacuolization_c 0.03331289 0.02538644 0.04029304
# 11: loss_of_structure 0.01487420 0.02375927 0.03349032
# 12:        necrosis_i 0.01124783 0.02371494 0.01988488


#compare training accuracy

pred <- predict(xgb_model, dtrain)

pred_matrix <- matrix(pred,
                      ncol = length(levels(SCTLD_subset_clean$worded_health)),
                      byrow = TRUE)

pred_class <- max.col(pred_matrix) - 1

mean(pred_class == y)
# 0.3673469

#I think we need to combine HD & DD because neither model predicts the differences between those very well and each only has n = 10 but HH has n = 29. So maybe if we combine those we'll get a better fit? 

#now let's model HH vs DD+HD 

#create a binary health variable 
SCTLD$health_binary <- ifelse(
  SCTLD$worded_health %in% c("HD", "DD"),
  "Diseased",
  "Healthy"
)

SCTLD$health_binary <- as.factor(SCTLD$health_binary)

#check it
table(SCTLD$health_binary) #looks right 

#re subset for this one 

SCTLD_subset_bin <- SCTLD[, c(
  "health_binary",
  "necrosis_c", "vacuolization_c", "exocytosis_c", "gastro_sep_c",
  "degraded_symb_c", "necrosis_i", "vacuolization_i", "exocytosis_i",
  "gastro_sep_i", "degraded_symb_i",
  "fungus_sponge", "blown_out_gastro", "amoebocytes",
  "loss_of_eosin", "loss_of_structure"
)]

SCTLD_subset_bin_clean <- na.omit(SCTLD_subset_bin)

#random forest model

set.seed(123)

rf_binary <- randomForest(
  health_binary ~ .,
  data = SCTLD_subset_bin_clean,
  importance = TRUE,
  ntree = 500
)

print(rf_binary)

# Call:
#   randomForest(formula = health_binary ~ ., data = SCTLD_subset_bin_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 42.86%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased        9      11   0.5500000
# Healthy        10      19   0.3448276

varImpPlot(rf_binary)


#build xgb model

#build matrix 
predictor_vars <- colnames(SCTLD_subset_bin_clean)
predictor_vars <- predictor_vars[predictor_vars != "health_binary"]

X <- as.matrix(SCTLD_subset_bin_clean[, predictor_vars])

y <- ifelse(SCTLD_subset_bin_clean$health_binary == "Diseased", 1, 0)

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

mean(pred_class == y) # 0.9795918  #accuracy jumped wayyyyyyy up 


#importance
imp_bin <- xgb.importance(model = xgb_binary)
imp_bin

# Feature       Gain      Cover  Frequency
# <char>      <num>      <num>      <num>
#   1:    gastro_sep_i 0.29818773 0.17473098 0.19357977
# 2: degraded_symb_c 0.24223024 0.19455249 0.16926070
# 3:      necrosis_c 0.11116192 0.16454130 0.11867704
# 4: vacuolization_i 0.07897509 0.07733368 0.09922179
# 5:    exocytosis_c 0.05898106 0.04478143 0.07003891
# 6:   loss_of_eosin 0.05271669 0.10304827 0.08171206
# 7:    exocytosis_i 0.04768967 0.05015682 0.04863813
# 8: vacuolization_c 0.04387358 0.03937138 0.06517510
# 9: degraded_symb_i 0.03408971 0.05139604 0.06809339
# 10:      necrosis_i 0.02012226 0.08225405 0.06517510
# 11:    gastro_sep_c 0.01197204 0.01783356 0.02042802



#so let's increase ntree to test stability of rf model

rf_model <- randomForest(
  health_binary ~ .,
  data = SCTLD_subset_bin_clean,
  importance = TRUE,
  ntree = 1000
)

print(rf_model)

# Call:
#   randomForest(formula = health_binary ~ ., data = SCTLD_subset_bin_clean,      importance = TRUE, ntree = 1000) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 42.86%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased        8      12   0.6000000
# Healthy         9      20   0.3103448




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


# iter train_error_mean train_error_std test_error_mean test_error_std
# <int>            <num>           <num>           <num>          <num>
#   1:     1        0.3771795      0.05745249       0.4244444      0.2442423
# 2:     2        0.3617949      0.07855324       0.4644444      0.2126261
# 3:     3        0.3166667      0.07949752       0.5044444      0.2162532
# 4:     4        0.3166667      0.07949752       0.5044444      0.2162532
# 5:     5        0.3012821      0.06247945       0.5044444      0.2162532
# 6:     6        0.2911538      0.06577373       0.4844444      0.2323524

# train_error_mean ≈ 0.38
# test_error_mean  ≈ 0.42
#train error is decreasing and test error is increasing -> overfitting probably due to a small sample size 
table(SCTLD_subset_bin_clean$health_binary)
#CV accuracy ≈ 50–58%--> also consistent with RF OOB 57%; Test error unstable (high std ~0.21–0.24)

##############
#We could say: 
#Distinguish WPH moderately well; NOT clearly distinguish SCTLD from healthy

################
#so can we combine with total dataset if it doesn't clearly distinguish SCTLD from healthy......yes? maybe if the question is "Is there a shared histological signature of disease"