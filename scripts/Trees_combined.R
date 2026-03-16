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


############################################################
#First we need to combine :WPH dataset + SCTLD dataset
############################################################

WPH <- read_csv("./data/WPH_HistoScore.csv")
SCTLD <- read_csv("./data/Serial_HistoScore.csv")
DRTO <- read_csv("./data/DRTO_histoscore.csv")
FGB <- read_csv("./data/FGB_histo_score.csv")

str(FGB)
str(WPH)
str(SCTLD)

#standardize the health column
WPH$health_state <- WPH$Sample_Health_State
SCTLD$health_state <- SCTLD$worded_health
FGB$health_state <- FGB$'health state'

#combine DRTO & SCTLD
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

FGB <- FGB %>%
  rename( health_state = `health state`)

WPH <- WPH %>%
  rename( health_state = `Sample_Health_State`)

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

#redo LMH to be 123 for FGB

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

SCTLD$class <- ifelse(SCTLD$worded_health == "Healthy", "Healthy", "SCTLD")

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
# Healthy   SCTLD     WPH 
# 80     114      18 

#assign predictor variables; there's more than not so we'll just say which not to include 
predictor_vars <- training_data |>
  select(-class, -Individual, -Initials, -date_completed,  -Pics_taken,  -measurements, -dataset, -DONE,-Notes, -score) |>
  colnames()

#lets check/ Should be 15
predictor_vars #checks out 

#remove Nas
training_data <- training_data %>%
  drop_na(all_of(predictor_vars)) #only lost 3 

#recheck everything 
dim(training_data)
table(training_data$class) #lost the 3 from SCTLD
# Healthy   SCTLD     WPH 
# 80     111      18 

##########################################################################
#Now let's build the final model based on the SCTLD and WPH only models
##########################################################################

#we'll use the top 7 predictors from the SCTLD dataset 
#we may change this after we get more WPH samples and top 7 may be different or more predictors may be necessary

top7_vars <- c(
  "loss_of_structure",
  "vacuolization_I",
  "gastro_sep_C",
  "loss_of_eosin",
  "amoebocytes",
  "degraded_symb_C",
  "degraded_symb_I"
)

#train the model

rf_model <- randomForest(
  x = training_data[, top7_vars],
  y = training_data$class,
  ntree = 1000,
  importance = TRUE
)

#Now apply model to FGB
FGB_clean <- FGB %>%
  drop_na(all_of(top7_vars)) #only lost 1 

#predict
pred <- predict(rf_model, FGB_clean[, top7_vars])

#attach predictions
FGB_clean$predicted_class <- pred

#ahhh lets see 
table(FGB_clean$predicted_class)

# Healthy   SCTLD     WPH 
# 64      96       1 

#model is favoring SCTLD because the sample size is so small 
#so we'll run a model that balances classes and see how much better we get 

#retrain model
rf_model <- randomForest(
  x = training_data[, top7_vars],
  y = training_data$class,
  ntree = 1000,
  importance = TRUE,
  classwt = c(
    Healthy = 1,
    SCTLD = 1,
    WPH = 4 #start with 4 and see how that goes 
  )
)

#check predictions
pred <- predict(rf_model, training_data[, top7_vars])

table(pred)
#pred
# Healthy   SCTLD     WPH 
# 50      75      84 

#now lets check model performance
confusionMatrix(pred, training_data$class)

#Confusion Matrix and Statistics

# Reference
# Prediction Healthy SCTLD WPH
# Healthy      37    13   0
# SCTLD        11    64   0
# WPH          32    34  18
# 
# Overall Statistics
# 
# Accuracy : 0.5694          
# 95% CI : (0.4993, 0.6375)
# No Information Rate : 0.5311          
# P-Value [Acc > NIR] : 0.1492          
# 
# Kappa : 0.3697          
# 
# Mcnemar's Test P-Value : 2.823e-14       
# 
# Statistics by Class:
# 
#                      Class: Healthy Class: SCTLD Class: WPH
# Sensitivity                  0.4625       0.5766    1.00000
# Specificity                  0.8992       0.8878    0.65445
# Pos Pred Value               0.7400       0.8533    0.21429
# Neg Pred Value               0.7296       0.6493    1.00000
# Prevalence                   0.3828       0.5311    0.08612
# Detection Rate               0.1770       0.3062    0.08612
# Detection Prevalence         0.2392       0.3589    0.40191
# Balanced Accuracy            0.6809       0.7322    0.82723


###########################################################
#Interpretation
#model is actively predicting WPH. So the weighting did its job.
#Problem: The model predicts WPH for many Healthy and SCTLD samples.
#Sensitivity/Specificity/Precision show -> The model finds all WPH samples -> But it over-predicts WPH -> happens because the class weight forces the model to prefer WPH predictions.
#WPH histology overlaps strongly with both Healthy and SCTLD features.
#Balanced Accuracy: means the model does detect WPH signals, but they are not unique.
###########################################################


###############################################################################
#So we'll keep the weight but reduce them to balance sensitivity and precision
###############################################################################

#retrain model
rf_model <- randomForest(
  x = training_data[, top7_vars],
  y = training_data$class,
  ntree = 1000,
  importance = TRUE,
  classwt = c(Healthy = 1, SCTLD = 1, WPH = 2)
)

#check predictions
pred <- predict(rf_model, training_data[, top7_vars])

table(pred)

# pred
# Healthy   SCTLD     WPH 
# 63      77      69 

#now lets check model performance
confusionMatrix(pred, training_data$class)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction Healthy SCTLD WPH
# Healthy      45    18   0
# SCTLD        11    66   0
# WPH          24    27  18
# 
# Overall Statistics
# 
# Accuracy : 0.6172          
# 95% CI : (0.5476, 0.6834)
# No Information Rate : 0.5311          
# P-Value [Acc > NIR] : 0.007391        
# 
# Kappa : 0.4205          
# 
# Mcnemar's Test P-Value : 2.135e-11       
# 
# Statistics by Class:
# 
#                      Class: Healthy Class: SCTLD Class: WPH
# Sensitivity                  0.5625       0.5946    1.00000
# Specificity                  0.8605       0.8878    0.73298
# Pos Pred Value               0.7143       0.8571    0.26087
# Neg Pred Value               0.7603       0.6591    1.00000
# Prevalence                   0.3828       0.5311    0.08612
# Detection Rate               0.2153       0.3158    0.08612
# Detection Prevalence         0.3014       0.3684    0.33014
# Balanced Accuracy            0.7115       0.7412    0.86649

###########################################################
#Interpretation
#This is better! -> we improved healthy & SCTLD sensitivity but still predict all WPH samples 
# #True class	Predicted WPH
# Healthy	24
# SCTLD	27
# WPH	18
#so it misclassified 24 healthy and 27 SCTLD but in the big pic that's not too many 
#still means WPH histology overlaps strongly with both healthy and SCTLD features
###########################################################

#lets make a plot now 

#run a PCA


X <- training_data[, top7_vars]

pca <- prcomp(X, scale. = TRUE)
pca_df <- data.frame(pca$x)

pca_df$class <- training_data$class

#plot PCA
ggplot(pca_df, aes(PC1, PC2, color = class)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Histological Feature Space",
    x = "PC1",
    y = "PC2"
  )

#lots of overlap!

#now lets add FGB samples to PCA and see where they fall

#first we have to run on training data only 


X_train <- training_data[, top7_vars]

pca <- prcomp(X_train, scale. = TRUE)

pca_train <- as.data.frame(pca$x)
pca_train$class <- training_data$class
pca_train$type <- "Training"

#project FGB samples on PCA plot
X_fgb <- FGB_clean[, top7_vars]

pca_fgb <- predict(pca, newdata = X_fgb)

pca_fgb <- as.data.frame(pca_fgb)
pca_fgb$class <- "FGB"
pca_fgb$type <- "FGB"

#combine 
pca_all <- bind_rows(pca_train, pca_fgb)

#plot
ggplot() +
  geom_point(
    data = pca_train,
    aes(PC1, PC2, color = class),
    shape = 16,
    alpha = 0.4,
    size = 2
  ) +
  geom_point(
    data = pca_fgb,
    aes(PC1, PC2),
    color = "black",
    shape = 17,
    size = 4
  ) +
  stat_ellipse(
    data = pca_train,
    aes(PC1, PC2, color = class),
    level = 0.95,
    linewidth = 1
  ) +
  theme_minimal() +
  labs(
    title = "PCA of Histological Predictors",
    subtitle = "Black triangles = FGB samples",
    x = "PC1",
    y = "PC2"
  )
######################################################################
#PCA matches rf showing that all samples overlap (Accuracy 0.62 for rf model)
#if we squint: Many FGB points sit within the SCTLD ellipse, some fall toward the WPH region, few sit strongly in the healthy-only region
#FGB lesions look more diseased than healthy, but do not strongly match a single disease type
######################################################################




#Histological predictors show substantial overlap between SCTLD, WPH, and healthy samples. FGB lesions occupy the same feature space and do not form a distinct cluster, suggesting that FGB histology shares characteristics with both disease states rather than representing a clearly separable pathology.