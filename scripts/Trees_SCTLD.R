#################################
#Let's build some decision trees 
#SCTLD
#################################

#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)
library(caret)
library(pROC)

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

################################################################################
#Now we're going to add in DRTO samples to up sample size and see if that helps
#03/10/2026
################################################################################

DRTO <- read_csv("./data/DRTO_histoscore.csv")
 

#Rename columns to match existing SCTLD dataset
DRTO <- DRTO %>%
  rename(
    Sample_Health_State = `sample health state`,
    gastro_sep_c = `gastro sep_C`,
    degraded_symb_c = `degraded symb_C`,
    gastro_sep_i = `gastro sep_I`,
    degraded_symb_i = `degraded symb_I`,
    exocytosis_i = `exocytosis_I`,
    vacuolization_i = `vacuolization_I`,
    necrosis_i =`necrosis_I`,
    necrosis_c = `necrosis_C`,
    vacuolization_c = `vacuolization_C`,
    exocytosis_c = `exocytosis_C`
    
  )

#check
colnames(SCTLD)
colnames(DRTO)

#keep only common columns 
common_cols <- intersect(names(SCTLD), names(DRTO))

DRTO <- DRTO[, common_cols]
SCTLD <- SCTLD[, common_cols]

#combine datasets
SCTLD <- bind_rows(SCTLD, DRTO)

#check class balance 
table(SCTLD$worded_health)

#make species column
SCTLD$species <- substr(SCTLD$Individual, 1, 4) #pulls out the first 4 letters since they all use 4 letter species code 

#make all health cases in worded_health match 
SCTLD <- SCTLD %>%
  mutate(
    worded_health = case_when(
      worded_health %in% c("HH", "naive") ~ "Healthy",
      worded_health %in% c("DD", "HD", "H_exposed") ~ "Diseased",
      TRUE ~ worded_health
    )
  )

#convert to factor
SCTLD$worded_health <- as.factor(SCTLD$worded_health)
table(SCTLD$worded_health) #114 diseased & 65 healthy


#######################################################
#now that everything is combined let's rerun the model
#######################################################
SCTLD <- SCTLD %>%
  rename(fungus_sponge = `fungus/sponge`)

predictor_vars <- c(
  "necrosis_c","vacuolization_c","exocytosis_c","gastro_sep_c",
  "degraded_symb_c","necrosis_i","vacuolization_i","exocytosis_i",
  "gastro_sep_i","degraded_symb_i",
  "fungus_sponge","blown_out_gastro","amoebocytes",
  "loss_of_eosin","loss_of_structure"
)

#remove NAs
SCTLD_clean <- na.omit(SCTLD[, c("worded_health", predictor_vars)]) #lost 3 

#Random forest model

set.seed(123)

rf_SCTLD <- randomForest(
  x = SCTLD_clean[, predictor_vars],
  y = as.factor(SCTLD_clean$worded_health),
  ntree = 1000,
  importance = TRUE
)

print(rf_SCTLD)
# Call:
#   randomForest(x = SCTLD_clean[, predictor_vars], y = as.factor(SCTLD_clean$worded_health),      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 39.77%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       91      20   0.1801802
# Healthy        50      15   0.7692308


varImpPlot(rf_SCTLD)

################################################################################
#Interpret
#this model has ~60% accuracy which is better than before! 
#Diseased Sample: 91 correct, 20 incorrect -> 91 / (91+20) = 82%
#Healthy Samples: 15 correct, 50 incorrect -> 15 / (15+50) = 23% 
#not so good on the healthy samples 
#Diseased = 111, Healthy = 65 -> total sample count so model is favoring disease 
#model is really good at detecting disease which could be good for combined model
################################################################################

##############################
#Now let's try balancing the classes to force the model to use the same number of healthy & diseased samples
##############################

set.seed(123)

rf_SCTLD_bal <- randomForest(
  x = SCTLD_clean[, predictor_vars],
  y = SCTLD_clean$worded_health,
  ntree = 1000,
  importance = TRUE,
  sampsize = c(65, 65)   # equal samples from each class
)

# Call:
#   randomForest(x = SCTLD_clean[, predictor_vars], y = SCTLD_clean$worded_health,      ntree = 1000, sampsize = c(65, 65), importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 41.48%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       77      34   0.3063063
# Healthy        39      26   0.6000000

print(rf_SCTLD_bal)
varImpPlot(rf_SCTLD_bal)


#########################################
#Interpretation
#More balanced and ~60% chance of being right
#maybe use this one and duck out some of the parameters that may not be contributing as much 
#########################################


#############################################################################################################################################
#Now we're going to make a loop that will test top 2,3,4,5...etc to top 15 (all) predictors and see which is the better fit 
#Still going to balance classes since there is some size difference between healthy and diseased 
#############################################################################################################################################

#run full model
response_var <- "worded_health"   # response variable

predictors <- setdiff(names(SCTLD_clean), response_var)

rf_full <- randomForest(
  as.formula(paste(response_var, "~ .")),
  data = SCTLD_clean,
  importance = TRUE
)

imp <- importance(rf_full)

imp_df <- data.frame(
  variable = rownames(imp),
  importance = imp[,1]
)

imp_df <- imp_df[order(-imp_df$importance), ] #ranked predictor list

#Set cross-validation 10-fold CV repeated 3 times

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Make sure outcome factor order works with ROC
SCTLD_clean$worded_health <- factor(
  SCTLD_clean$worded_health,
  levels = c("Healthy","Diseased")
)

#predictor subset loop

results <- data.frame()

for(i in 2:min(15, nrow(imp_df))){
  
  print(paste("Running model with", i, "predictors"))
  
  top_vars <- imp_df$variable[1:i]
  
  formula <- as.formula(
    paste(response_var, "~", paste(top_vars, collapse="+"))
  )
  
  rf_model <- train(
    formula,
    data = SCTLD_clean,
    method = "rf",
    metric = "ROC",
    trControl = ctrl
  )
  
  results <- rbind(
    results,
    data.frame(
      predictors = i,
      ROC = rf_model$results$ROC,
      Sensitivity = rf_model$results$Sens,
      Specificity = rf_model$results$Spec
    )
  )
}

print(results)



#####################################################
#Hopefully the last model
#now we'll 
#####################################################

SCTLD_clean$worded_health <- factor(
  SCTLD_clean$worded_health,
  levels = c("Healthy", "Diseased")
)

#define predictor variables

predictor_vars <- c(
  "necrosis_c","vacuolization_c","exocytosis_c","gastro_sep_c",
  "degraded_symb_c","necrosis_i","vacuolization_i","exocytosis_i",
  "gastro_sep_i","degraded_symb_i",
  "fungus_sponge","blown_out_gastro","amoebocytes",
  "loss_of_eosin","loss_of_structure"
)

#run recursive feature elimination --> find the best combination
#use Random Forest inside RFE


set.seed(123)

control <- rfeControl(
  functions = rfFuncs,
  method = "repeatedcv",
  number = 10,
  repeats = 3
)

rfe_model <- rfe(
  x = SCTLD_clean[, predictor_vars],
  y = SCTLD_clean$worded_health,
  sizes = 2:length(predictor_vars),
  rfeControl = control
)

print(rfe_model)

#Recursive feature selection

# Outer resampling method: Cross-Validated (10 fold, repeated 3 times) 
# 
# Resampling performance over subset size:
#   
#   Variables Accuracy    Kappa AccuracySD KappaSD Selected
# 2   0.6065 -0.03208    0.05601 0.08418         
# 3   0.6141 -0.02008    0.04822 0.08029         
# 4   0.5923  0.02620    0.10906 0.23384         
# 5   0.6090  0.05996    0.10563 0.23333         
# 6   0.6261  0.09104    0.09834 0.22024         
# 7   0.6347  0.12156    0.10016 0.22155        * 7 variables Accuracy = 0.6347
#   8   0.6094  0.06418    0.09928 0.20702         
# 9   0.6028  0.07545    0.09952 0.21707         
# 10   0.6141  0.09784    0.11005 0.23780         
# 11   0.6166  0.10412    0.09100 0.18737         
# 12   0.6305  0.12782    0.10140 0.21300         
# 13   0.6279  0.12065    0.09310 0.19427         
# 14   0.6107  0.07476    0.09432 0.20502         
# 15   0.6186  0.08272    0.09688 0.21240         
# 
# The top 5 variables (out of 7):
#   gastro_sep_i, loss_of_eosin, amoebocytes, exocytosis_i, loss_of_structure

rfe_model$optVariables
#"gastro_sep_i"      "loss_of_eosin"     "amoebocytes"       "exocytosis_i"      "loss_of_structure" "degraded_symb_c"   "necrosis_i"   
#model is best with this group of predictors 

################################
#now let's see if xgb improves 
################################

#pull out best variables
best_vars <- rfe_model$optVariables
best_vars

#build xgb 

# Create matrix with only selected variables
X <- as.matrix(SCTLD_clean[, best_vars])

# Convert outcome to 0/1
y <- ifelse(SCTLD_clean$worded_health == "Diseased", 1, 0)

dtrain <- xgb.DMatrix(data = X, label = y)

#train xgb with cv
set.seed(123)


cv <- xgb.cv(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8
  ),
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 20,
  verbose = 0
)

best_nrounds <- which.max(cv$evaluation_log$test_auc_mean)
best_nrounds

#train final xgb model
xgb_final <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8
  ),
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

#get cross validation AUC 
max(cv$evaluation_log$test_auc_mean)

#which predictors mattered 
importance_matrix <- xgb.importance(model = xgb_final)
print(importance_matrix)


######################################################################
#Interpret
#It is not outperforming RFE + Random Forest (~0.63 accuracy earlier)
#XGBoost does NOT improve performance over RF approach.
#RFE selected a stable 7-variable model -> RF performed around ~0.63 accuracy -> XGB AUC is ~0.60 -> Increasing model complexity did NOT improve performance
#suggest that the signal is real but moderate -> likely limited by sample size 
#there is overlap between RFE and XGB between predictors -> That overlap is VERY strong evidence of stability.
######################################################################


#####################################################################
#The simpler model is out RFE 7 variable model so we'll use that one and bootstrap
#How often each of the 7 predictors is important across 500 resamples.
#####################################################################

#This version records mean permutation importance (more stable than “top variable only”)

set.seed(123)

n_boot <- 500
importance_matrix <- matrix(0, nrow = n_boot, ncol = length(best_vars))
colnames(importance_matrix) <- best_vars

for(i in 1:n_boot){
  
  boot_index <- sample(1:nrow(SCTLD_clean), replace = TRUE)
  boot_data <- SCTLD_clean[boot_index, ]
  
  rf_boot <- randomForest(
    x = boot_data[, best_vars],
    y = boot_data$worded_health,
    ntree = 500,
    importance = TRUE
  )
  
  importance_matrix[i, ] <- importance(rf_boot)[,1]  # MeanDecreaseAccuracy
}

# Average importance across bootstraps
boot_mean <- colMeans(importance_matrix)

# Sort results
sort(boot_mean, decreasing = TRUE)

# gastro_sep_i   degraded_symb_c     loss_of_eosin       amoebocytes      exocytosis_i loss_of_structure        necrosis_i 
# 17.970157         17.855237         16.736840         14.949747         13.273175         11.444624          7.590522 

##################################
#Interpret
#top 4 variables consistently carry the signal
#we'll keep the 7 variable model to be more conservative 
##################################




####################################################################
#Now we clean it up to make one final clean model 
####################################################################

#refit final RF with the 7 variables 

set.seed(123)

final_rf <- randomForest(
  x = SCTLD_clean[, best_vars],
  y = SCTLD_clean$worded_health,
  ntree = 1000,
  importance = TRUE
)

final_rf 
# Call:
#   randomForest(x = SCTLD_clean[, best_vars], y = SCTLD_clean$worded_health,      ntree = 1000, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 32.39%
# Confusion matrix:
#   Healthy Diseased class.error
# Healthy       20       45   0.6923077
# Diseased      12       99   0.1081081


#get final variable importance 
importance(final_rf)
# Healthy   Diseased MeanDecreaseAccuracy MeanDecreaseGini
# gastro_sep_i      14.669336 -1.9164420             8.495903         6.298485
# loss_of_eosin     11.935699 -0.1091841             7.995638         4.239905
# amoebocytes        8.672226  5.4371679             9.328986         4.976324
# exocytosis_i       4.575044  4.7859070             6.691658         5.817685
# loss_of_structure  3.455470  8.9792406             9.600801         3.896807
# degraded_symb_c   11.867237  5.4117081            11.645803         6.690140
# necrosis_i        -6.122842  9.8812446             4.863348         4.862462
varImpPlot(final_rf)


#report cross validation performance

control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

set.seed(123)

final_model <- train(
  worded_health ~ .,
  data = SCTLD_clean[, c("worded_health", best_vars)],
  method = "rf",
  metric = "ROC",
  trControl = control,
  ntree = 1000
)

final_model

#Random Forest 

# 176 samples
# 7 predictor
# 2 classes: 'Healthy', 'Diseased' 
# 
# No pre-processing
# Resampling: Cross-Validated (10 fold, repeated 3 times) 
# Summary of sample sizes: 157, 159, 158, 159, 158, 159, ... 
# Resampling results across tuning parameters:
#   
#   mtry  ROC        Sens       Spec     
# 2     0.6071849  0.3087302  0.8472222
# 4     0.5901094  0.3341270  0.7962121
# 7     0.5867544  0.3690476  0.7871212
# 
# ROC was used to select the optimal model using the largest value.
# The final value used for the model was mtry = 2.

