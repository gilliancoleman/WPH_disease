##which disease signs matter?
DRTO <- read_csv("~/Documents/MOTE/DRTO_histoscore.csv")
# Load necessary library
install.packages("randomForest")
library(randomForest)
View(DRTO)
# Convert categorical variables to factors
DRTO$worded_health <- as.factor(DRTO$worded_health)  
DRTO$Individual <- as.factor(DRTO$Individual)  

table(DRTO$species, DRTO$site_name_NPS)
table(DRTO$species, DRTO$worded_health)
table(DRTO$worded_health, DRTO$sample_date)

# Select relevant columns for analysis
DRTO_subset <- DRTO[, c("worded_health",
                    "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                    "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                    "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                    "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure")]

# Remove rows with missing values
DRTO_subset_clean <- na.omit(DRTO_subset)

# Train Random Forest model
set.seed(123)
rf_model <- randomForest(worded_health ~ ., data =DRTO_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Plot feature importance
varImpPlot(rf_model)

# Load the required library
install.packages("gbm")
library(gbm)

# Train a Gradient Boosting model
gbm_model <- gbm(worded_health ~ ., data = DRTO_subset_clean, 
                 distribution = "multinomial", # For classification tasks
                 n.trees = 500, 
                 interaction.depth = 3,  # Depth of the trees
                 shrinkage = 0.1,        # Learning rate
                 cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##we learned that fungus/sponge, blown out gastro, and loss of structure weren't helping, so new subset!

# Select relevant columns for analysis
DRTO_subset2 <- DRTO[, c("worded_health",
                        "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                        "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                        "gastro_sep_I", "degraded_symb_I",
                        "amoebocytes", "loss_of_eosin")]

# Remove rows with missing values
DRTO_subset2_clean <- na.omit(DRTO_subset2)

# Train Random Forest model
set.seed(123)
rf_model <- randomForest(worded_health ~ ., data =DRTO_subset2_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model)

# Plot feature importance
varImpPlot(rf_model)


# Train a Gradient Boosting model
gbm_model <- gbm(worded_health ~ ., data = DRTO_subset2_clean, 
                 distribution = "multinomial", # For classification tasks
                 n.trees = 500, 
                 interaction.depth = 5,  # Depth of the trees
                 shrinkage = 0.1,        # Learning rate
                 cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##loss of eosin is lowest now

# Select relevant columns for analysis
DRTO_subset3 <- DRTO[, c("worded_health",
                         "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                         "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "gastro_sep_I", "degraded_symb_I",
                         "amoebocytes")]

# Remove rows with missing values
DRTO_subset3_clean <- na.omit(DRTO_subset3)

# Train Random Forest model
set.seed(123)
rf_model3 <- randomForest(worded_health ~ ., data = DRTO_subset3_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model3)

# Plot feature importance
varImpPlot(rf_model)


# Train a Gradient Boosting model
gbm_model3 <- gbm(worded_health ~ ., data = DRTO_subset3_clean, 
                 distribution = "multinomial", # For classification tasks
                 n.trees = 500, 
                 interaction.depth = 5,  # Depth of the trees
                 shrinkage = 0.1,        # Learning rate
                 cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model3)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")


##might do this one species at a time?
##MCAV subset
treeCNATsub <- subset(DRTO, species=="CNAT")
View(treeCNATsub)

# Select relevant columns for analysis
CNAT_subset2 <- treeCNATsub[, c("worded_health", "site_name_NPS", "monitored",
                         "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                         "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "gastro_sep_I", "degraded_symb_I",
                         "amoebocytes", "loss_of_eosin")]

# Remove rows with missing values
CNAT_subset2_clean <- na.omit(CNAT_subset2)

# Train Random Forest model
set.seed(123)
rf_model3 <- randomForest(worded_health ~ ., data =CNAT_subset2_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model3)

# Plot feature importance
varImpPlot(rf_model)


# Train a Gradient Boosting model
gbm_model3 <- gbm(worded_health ~ ., data = CNAT_subset2_clean, 
                 distribution = "multinomial", # For classification tasks
                 n.trees = 500, 
                 interaction.depth = 5,  # Depth of the trees
                 shrinkage = 0.1,        # Learning rate
                 cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model3)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##loss of eosin wasn't important
# Select relevant columns for analysis
CNAT_subset3 <- treeCNATsub[, c("worded_health", "site_name_NPS", "monitored",
                                "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                                "gastro_sep_I", "degraded_symb_I",
                                "amoebocytes")]

# Remove rows with missing values
CNAT_subset3_clean <- na.omit(CNAT_subset3)

# Train Random Forest model
set.seed(123)
rf_model4 <- randomForest(worded_health ~ ., data =CNAT_subset3_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model4)

# Plot feature importance
varImpPlot(rf_model4)


# Train a Gradient Boosting model
gbm_model4 <- gbm(worded_health ~ ., data = CNAT_subset3_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 3,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 2)           # Cross-validation

# View the summary of the model
summary(gbm_model4)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")


##what if we try not for worded health but sample health itself, make it binomial
# Select relevant columns for analysis
DRTO_subset4 <- DRTO[, c("sample_health_state",
                         "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                         "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "gastro_sep_I", "degraded_symb_I",
                         "amoebocytes")]

# Remove rows with missing values
DRTO_subset4_clean <- na.omit(DRTO_subset4)

# Train Random Forest model
set.seed(123)
rf_model4 <- randomForest(sample_health_state ~ ., data = DRTO_subset4_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model4)

# Plot feature importance
varImpPlot(rf_mode4)


# Train a Gradient Boosting model
gbm_model4 <- gbm(sample_health_state ~ ., data = DRTO_subset4_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 5,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model4)

# Make predictions
predictions <- predict(gbm_model4, newdata = DRTO_subset4_clean, n.trees = 500, type = "response")

##what if we try not for worded health but sample state, bring it back to HH, HD, DD
# Select relevant columns for analysis
DRTO_subset5 <- DRTO[, c("state",
                         "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                         "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "gastro_sep_I", "degraded_symb_I",
                         "amoebocytes")]

# Remove rows with missing values
DRTO_subset5_clean <- na.omit(DRTO_subset5)
DRTO_subset5_clean$state <- as.factor(DRTO_subset5_clean$state)
# Train Random Forest model
set.seed(123)
rf_model5 <- randomForest(state ~ ., data = DRTO_subset5_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model5)

# Plot feature importance
varImpPlot(rf_model5)


# Train a Gradient Boosting model
gbm_model5 <- gbm(state ~ ., data = DRTO_subset5_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 5,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model5)

# Make predictions
predictions <- predict(gbm_model4, newdata = DRTO_subset4_clean, n.trees = 500, type = "response")

##disease state makes this model run better - theres no difference between naive and H_exposed
# Select relevant columns for analysis
DRTO_subset6 <- DRTO[, c("state",
                        "exocytosis_C", "gastro_sep_C", 
                          "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "gastro_sep_I", "degraded_symb_I",
                         "amoebocytes")]

# Remove rows with missing values
DRTO_subset6_clean <- na.omit(DRTO_subset6)
DRTO_subset6_clean$state <- as.factor(DRTO_subset6_clean$state)
# Train Random Forest model
set.seed(123)
rf_model6 <- randomForest(state ~ ., data = DRTO_subset6_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model6)

# Plot feature importance
varImpPlot(rf_model6)


# Train a Gradient Boosting model
gbm_model6 <- gbm(state ~ ., data = DRTO_subset6_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 5,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model6)

# Make predictions
predictions <- predict(gbm_model6, newdata = DRTO_subset4_clean, n.trees = 500, type = "response")

##removing the bottom three didn't help make a better model but we're gonna try one more time
# Select relevant columns for analysis
DRTO_subset7 <- DRTO[, c("state",
                         "exocytosis_C", "gastro_sep_C", 
                         "necrosis_I", "vacuolization_I", "exocytosis_I", 
                         "amoebocytes")]

# Remove rows with missing values
DRTO_subset7_clean <- na.omit(DRTO_subset7)
DRTO_subset7_clean$state <- as.factor(DRTO_subset7_clean$state)
# Train Random Forest model
set.seed(123)
rf_model7 <- randomForest(state ~ ., data = DRTO_subset7_clean, importance = TRUE, ntree = 500)

# Print model summary
print(rf_model7)

# Plot feature importance
varImpPlot(rf_model7)


# Train a Gradient Boosting model
gbm_model7 <- gbm(state ~ ., data = DRTO_subset7_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 5,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(gbm_model7)

# Make predictions
predictions <- predict(gbm_model7, newdata = DRTO_subset4_clean, n.trees = 500, type = "response")
