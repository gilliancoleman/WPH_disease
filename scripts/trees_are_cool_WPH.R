#################################
#Let's build some decision trees 
#################################

#load libraries
library(tidyverse)
library(randomForest)
library(xgboost)

############
#WPH First
############

#which disease signs matter?
WPH <- read_csv("./data/WPH_HistoScore.csv")

#Convert categorical variables to factors
WPH$Sample_Health_State <- as.factor(WPH$Sample_Health_State)  
WPH$Individual <- as.factor(WPH$Individual)  

#do some tabling
table(WPH$species, WPH$Site)
table(WPH$species, WPH$Sample_Health_State)
table(WPH$Sample_Health_State, DRTO$sample_date) #do i need this one?

#rename fungus/sponge
WPH <- WPH %>%
  rename(fungus_sponge = `fungus/sponge`)


# Select relevant columns for analysis
WPH_subset <- WPH[, c("Sample_Health_State",
                        "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                        "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                        "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                        "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure")]

#Remove rows with missing values
WPH_subset_clean <- na.omit(WPH_subset) #wooooo none 

#Train Random Forest model
set.seed(123)
rf_model <- randomForest(Sample_Health_State ~ ., data = WPH_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Call:
#   randomForest(formula = Sample_Health_State ~ ., data = WPH_subset_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 54.55%
# Confusion matrix:
#   DD HD HH class.error
# DD  0  2  5   1.0000000
# HD  2  4  5   0.6363636
# HH  0  4 11   0.2666667

# Plot feature importance
varImpPlot(rf_model)

#going to use xgboost instead of gbm because gmb is out of date



# Ensure factor with explicit levels
WPH_subset_clean$Sample_Health_State <- factor(
  WPH_subset_clean$Sample_Health_State
)

X <- model.matrix(Sample_Health_State ~ . - 1,
                  data = WPH_subset_clean)

y <- as.integer(WPH_subset_clean$Sample_Health_State) - 1


#Uses params explicitly → no guessing
dtrain <- xgb.DMatrix(
  data = X,
  label = y
)

params <- list(
  objective = "multi:softprob",
  num_class = length(levels(WPH_subset_clean$Sample_Health_State)),
  max_depth = 3,
  learning_rate = 0.1,
  eval_metric = "mlogloss"
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#Get model summary
xgb_model #just a header

imp <- xgb.importance(model = xgb_model)
print(imp) #tells us which predictors are driving splits
#1.) feature: predictor variable used in the model
#2.) Gain How much this variable improves model performance when it is used *very important* Higher = more informative for distinguishing health states
#3.) Cover: How much of the data is affected by splits on this variable
#4.) Frequency: How often the variable is used in trees

#lets looks at this one
# Feature        Gain       Cover   Frequency
# <char>       <num>       <num>       <num>
#   1:   vacuolization_C 0.186161243 0.101332862 0.111877395
# 2:   degraded_symb_I 0.160217129 0.092702682 0.075478927
# 3:      exocytosis_C 0.120776391 0.094410954 0.084291188
# 4:   degraded_symb_C 0.120629128 0.120208438 0.111877395
# 5:     loss_of_eosin 0.098696048 0.130083426 0.122605364
# 6:        necrosis_C 0.081623479 0.098646618 0.101149425
# 7:      gastro_sep_I 0.073659197 0.112491979 0.130268199
# 8:   vacuolization_I 0.050985140 0.072616575 0.080076628
# 9:      gastro_sep_C 0.040734991 0.085736952 0.080459770
# 10:      exocytosis_I 0.037650735 0.053142061 0.045210728
# 11:        necrosis_I 0.020993157 0.030443578 0.047509579
# 12: loss_of_structure 0.004449397 0.006694434 0.007279693
# 13:       amoebocytes 0.003423967 0.001489441 0.001915709

###########################
#Interpretation
#Vacuolization_C strongly separates health states Gain ~0.19
#degraded_symb_I also strongly associated with shifts in health state Gain ~0.16
#exocytosis_C & degraded_symb_C both cell degradation processes Gain ~0.12

#loss_of_structure & amoebocytes minimal contribution

#also notice that most of the _C variables contribute more than _I variables so maybe condition consistency is more predictive than mere intensity?
xgb.plot.importance(imp)
################################

#lets make a new subset without loss_of_structure & amoebocytes

WPH_subset_2 <- WPH[, c("Sample_Health_State",
                      "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                      "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                      "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                      "blown_out_gastro", "loss_of_eosin")]

#Remove rows with missing values
WPH_subset_clean_2 <- na.omit(WPH_subset_2) #lost 2 


# Train Random Forest model
set.seed(123)
rf_model_2 <- randomForest(Sample_Health_State ~ ., data =WPH_subset_clean_2, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model_2)

#Call:
# randomForest(formula = Sample_Health_State ~ ., data = WPH_subset_clean_2,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 60.61%
# Confusion matrix:
#   DD HD HH class.error
# DD  0  4  3   1.0000000
# HD  2  4  5   0.6363636
# HH  2  4  9   0.4000000


# Plot feature importance
varImpPlot(rf_model_2)


# Ensure factor with explicit levels
WPH_subset_clean_2$Sample_Health_State <- factor(
  WPH_subset_clean_2$Sample_Health_State
)

X2 <- model.matrix(Sample_Health_State ~ . - 1,
                  data = WPH_subset_clean_2)

y2 <- as.integer(WPH_subset_clean_2$Sample_Health_State) - 1


#Uses params explicitly → no guessing
dtrain <- xgb.DMatrix(
  data = X2,
  label = y2
)

params2 <- list(
  objective = "multi:softprob",
  num_class = length(levels(WPH_subset_clean$Sample_Health_State)),
  max_depth = 3,
  learning_rate = 0.1,
  eval_metric = "mlogloss"
)

xgb_model2 <- xgb.train(
  params2 = params2,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#Get model summary
xgb_model2 #just a header

imp2 <- xgb.importance(model = xgb_model2)
print(imp2) #tells us which predictors are driving splits

#           Feature        Gain       Cover   Frequency
# <char>       <num>       <num>       <num>
#   1: degraded_symb_I 0.221967999 0.075175358 0.049707602
# 2:      necrosis_C 0.220668360 0.138456847 0.239766082
# 3: vacuolization_C 0.187416021 0.197163769 0.137426901
# 4: degraded_symb_C 0.118046254 0.122445868 0.086257310
# 5:   loss_of_eosin 0.090991300 0.060231778 0.046783626
# 6: vacuolization_I 0.038359349 0.041323574 0.040935673
# 7:    gastro_sep_I 0.037049383 0.091338823 0.073099415
# 8:    gastro_sep_C 0.036705462 0.080664837 0.109649123
# 9:    exocytosis_C 0.027002861 0.075937786 0.135964912
# 10:    exocytosis_I 0.018426032 0.091033852 0.046783626
# 11:      necrosis_I 0.003246596 0.023330284 0.027777778
# 12:   fungus_sponge 0.000120383 0.002897225 0.005847953

###########################
#Interpretation
#necrosis_I & fungus_sponge have minimal contribution
###########################

#remove those and do a 3rd 

# Select relevant columns for analysis
WPH_subset_3 <- WPH[, c("Sample_Health_State",
                      "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                      "degraded_symb_C", "vacuolization_I", "exocytosis_I", 
                      "gastro_sep_I", "degraded_symb_I", 
                      "blown_out_gastro","loss_of_eosin")]

#Remove rows with missing values
WPH_subset_clean_3 <- na.omit(WPH_subset_3) #wooooo none 

#Train Random Forest model
set.seed(123)
rf_model_3 <- randomForest(Sample_Health_State ~ ., data = WPH_subset_clean_3, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model_3)

# Call:
#   randomForest(formula = Sample_Health_State ~ ., data = WPH_subset_clean_3,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 60.61%
# Confusion matrix:
#   DD HD HH class.error
# DD  0  3  4   1.0000000
# HD  2  3  6   0.7272727
# HH  2  3 10   0.3333333


# Plot feature importance
varImpPlot(rf_model_3)


# Ensure factor with explicit levels
WPH_subset_clean_3$Sample_Health_State <- factor(
  WPH_subset_clean_3$Sample_Health_State
)

X3 <- model.matrix(Sample_Health_State ~ . - 1,
                  data = WPH_subset_clean_3)

y3 <- as.integer(WPH_subset_clean_3$Sample_Health_State) - 1


#Uses params explicitly → no guessing
dtrain <- xgb.DMatrix(
  data = X3,
  label = y3
)

params3 <- list(
  objective = "multi:softprob",
  num_class = length(levels(WPH_subset_clean_3$Sample_Health_State)),
  max_depth = 3,
  learning_rate = 0.1,
  eval_metric = "mlogloss"
)

xgb_model3 <- xgb.train(
  params = params3,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#Get model summary
xgb_model3 #just a header

imp3 <- xgb.importance(model = xgb_model3)
print(imp3)

# Feature       Gain      Cover  Frequency
# <char>      <num>      <num>      <num>
#   1: degraded_symb_I 0.22245780 0.08507440 0.04619970
# 2:      necrosis_C 0.22139527 0.13422516 0.23994039
# 3: vacuolization_C 0.18673435 0.20622276 0.15201192
# 4: degraded_symb_C 0.12124834 0.09935368 0.10134128
# 5:   loss_of_eosin 0.08990643 0.07380129 0.05961252
# 6:    gastro_sep_C 0.03890584 0.10882309 0.11177347
# 7: vacuolization_I 0.03819806 0.07004359 0.06110283
# 8:    gastro_sep_I 0.03665844 0.08417255 0.09388972
# 9:    exocytosis_C 0.02611837 0.10626785 0.11028316
# 10:    exocytosis_I 0.01837710 0.03201563 0.02384501

################################
#Interpretation
#exocytosis_I minimal contribution
#################################

#do for PAST and OFRA seperately 

#OFRA first

treeOFRAsub <- subset(WPH, species=="OFRA")
View(treeOFRAsub)

# Select relevant columns for analysis
OFRA_subset <- treeOFRAsub[, c("Sample_Health_State",
                               "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                               "degraded_symb_C", "vacuolization_I", "exocytosis_I", 
                               "gastro_sep_I", "degraded_symb_I", 
                               "blown_out_gastro")]

# Remove rows with missing values
OFRA_subset_clean <- na.omit(OFRA_subset)

# Train Random Forest model
set.seed(123)
rf_model_4 <- randomForest(Sample_Health_State ~ ., data =OFRA_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model_4)

#Call:
# randomForest(formula = Sample_Health_State ~ ., data = OFRA_subset_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 80%
# Confusion matrix:
#   DD HD HH class.error
# DD  0  1  4         1.0
# HD  1  0  4         1.0
# HH  2  4  4         0.6

# Plot feature importance
varImpPlot(rf_model_4)

# Ensure factor with explicit levels
OFRA_subset_clean$Sample_Health_State <- factor(
  OFRA_subset_clean$Sample_Health_State
)

# train eXtreme Gradient Boosting model
X4 <- model.matrix(Sample_Health_State ~ . - 1,
                   data = OFRA_subset_clean)

y4 <- as.integer(OFRA_subset_clean$Sample_Health_State) - 1


#Uses params explicitly → no guessing
dtrain <- xgb.DMatrix(
  data = X4,
  label = y4
)

params4 <- list(
  objective = "multi:softprob",
  num_class = length(levels(OFRA_subset_clean$Sample_Health_State)),
  max_depth = 3,
  learning_rate = 0.1,
  eval_metric = "mlogloss"
)

xgb_model4 <- xgb.train(
  params = params4,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#Get model summary
xgb_model4 #just a header

imp4 <- xgb.importance(model = xgb_model4)
print(imp4)

# Feature        Gain      Cover  Frequency
# <char>       <num>      <num>      <num>
#   1:      necrosis_C 0.236025970 0.20658282 0.24689165
# 2:    gastro_sep_C 0.188751418 0.10275622 0.12078153
# 3: vacuolization_C 0.178479858 0.18972438 0.24156306
# 4:    exocytosis_C 0.114208162 0.15226117 0.06927176
# 5:    gastro_sep_I 0.113287065 0.07546160 0.09413854
# 6: degraded_symb_I 0.072484843 0.07251806 0.03552398
# 7:    exocytosis_I 0.046800748 0.02087236 0.01953819
# 8: degraded_symb_C 0.045995407 0.13727589 0.11900533
# 9: vacuolization_I 0.003966528 0.04254750 0.05328597

###########################################
#Interpretation
#vacuolization_I has minimal contribution
###########################################

#go ahead and do it for PAST just to see if there's any major differences

treePASTsub <- subset(WPH, species=="PAST")
View(treePASTsub)

# Select relevant columns for analysis
PAST_subset <- treePASTsub[, c("Sample_Health_State",
                               "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                               "degraded_symb_C", "vacuolization_I", "exocytosis_I", 
                               "gastro_sep_I", "degraded_symb_I", 
                               "blown_out_gastro")]

# Remove rows with missing values
PAST_subset_clean <- na.omit(PAST_subset)

# Train Random Forest model
set.seed(123)
rf_model_5 <- randomForest(Sample_Health_State ~ ., data =PAST_subset_clean, importance = TRUE, ntree = 500)



# Print model summary
print(rf_model_5)

#Call:
# randomForest(formula = Sample_Health_State ~ ., data = OFRA_subset_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 80%
# Confusion matrix:
#   DD HD HH class.error
# DD  0  1  4         1.0
# HD  1  0  4         1.0
# HH  2  4  4         0.6

# Plot feature importance
varImpPlot(rf_model_5)

# Ensure factor with explicit levels
PAST_subset_clean$Sample_Health_State <- factor(
  PAST_subset_clean$Sample_Health_State
)

# train eXtreme Gradient Boosting model
X5 <- model.matrix(Sample_Health_State ~ . - 1,
                   data = PAST_subset_clean)

y5 <- as.integer(PAST_subset_clean$Sample_Health_State) - 1


#Uses params explicitly → no guessing
dtrain <- xgb.DMatrix(
  data = X5,
  label = y5
)

params5 <- list(
  objective = "multi:softprob",
  num_class = length(levels(PAST_subset_clean$Sample_Health_State)),
  max_depth = 3,
  learning_rate = 0.1,
  eval_metric = "mlogloss"
)

xgb_model5 <- xgb.train(
  params = params5,
  data = dtrain,
  nrounds = 300,
  verbose = 0
)

#Get model summary
xgb_model5 #just a header

imp5 <- xgb.importance(model = xgb_model5)
print(imp5)

# Feature         Gain      Cover  Frequency
# <char>        <num>      <num>      <num>
#   1: degraded_symb_I 0.4468022117 0.40311005 0.23371648
# 2: degraded_symb_C 0.3352424389 0.15251196 0.21839080
# 3:    exocytosis_I 0.1140321845 0.05382775 0.06896552
# 4:    exocytosis_C 0.0700790674 0.11782297 0.17241379
# 5: vacuolization_C 0.0234351536 0.08433014 0.09578544
# 6: vacuolization_I 0.0077626710 0.09928230 0.12643678
# 7:      necrosis_C 0.0018346649 0.07476077 0.07279693
# 8:    gastro_sep_I 0.0008116081 0.01435407 0.01149425

###########################################
#Interpretation
#gastro_sep_I has minimal contribution
###########################################

#looks like there are some differences between species but is it enough to matter? - Maybe but I think sample sizes might be too small



##############################################
#let's try starting from scratch and running the original model to test accuracy and then see if the binary one makes the most sense 
##############################################

#create a binary health variable 
WPH$health_binary <- ifelse(
  WPH$Sample_Health_State %in% c("HD", "DD"),
  "Diseased",
  "Healthy"
)

WPH$health_binary <- as.factor(WPH$health_binary)

#check it
table(WPH$health_binary) #looks right & more balanced now too

# Select relevant columns for analysis
WPH_subset_bin <- WPH[, c("health_binary",
                      "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                      "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                      "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                      "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure")]

#Remove rows with missing values
WPH_subset_bin_clean <- na.omit(WPH_subset_bin) #wooooo none 

#Train Random Forest model
set.seed(123)

rf_binary <- randomForest(
  health_binary ~ .,
  data = WPH_subset_bin_clean,
  importance = TRUE,
  ntree = 500
)

print(rf_binary)

# Call:
#   randomForest(formula = health_binary ~ ., data = WPH_subset_bin_clean,      importance = TRUE, ntree = 500) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 42.42%
# Confusion matrix:
#   Diseased Healthy class.error
# Diseased       12       6   0.3333333
# Healthy         8       7   0.5333333


# Plot feature importance
varImpPlot(rf_model)

#build xgb model

#build matrix 
predictor_vars <- colnames(WPH_subset_bin_clean)
predictor_vars <- predictor_vars[predictor_vars != "health_binary"]

X <- as.matrix(WPH_subset_bin_clean[, predictor_vars])

y <- ifelse(WPH_subset_bin_clean$health_binary == "Diseased", 1, 0)

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

mean(pred_class == y) # 0.9393939  #accuracy jumped wayyyyyyy up 


#importance
imp_bin <- xgb.importance(model = xgb_binary)
imp_bin

# Feature       Gain      Cover  Frequency
# <char>      <num>      <num>      <num>
#   1: degraded_symb_C 0.21557734 0.24231944 0.24485126
# 2:      necrosis_C 0.18234801 0.11671510 0.11899314
# 3:    exocytosis_C 0.14634751 0.08395266 0.09382151
# 4:    gastro_sep_I 0.11029883 0.11338053 0.11670481
# 5: vacuolization_I 0.09925477 0.11589783 0.11098398
# 6: vacuolization_C 0.07451302 0.11169531 0.08123570
# 7:   loss_of_eosin 0.04714246 0.03225444 0.04347826
# 8:    exocytosis_I 0.04281158 0.07724291 0.05377574
# 9: degraded_symb_I 0.03930996 0.05259681 0.06178490
# 10:    gastro_sep_C 0.02571260 0.03105129 0.04462243
# 11:      necrosis_I 0.01668393 0.02289369 0.02974828