##we want to combine the histo measures and histo score into one frame and try the classification on that
##and we're not going to cry

df4$individual <- trimws(df4$individual)
df4$worded_health <- trimws(df4$worded_health)
df4$sample_health_state <- trimws(df4$sample_health_state)
df4$species <- trimws(df4$species)

DRTO$individual <- trimws(DRTO$individual)
DRTO$worded_health <- trimws(DRTO$worded_health)
DRTO$sample_health_state <- trimws(DRTO$sample_health_state)
DRTO$species <- trimws(DRTO$species)

df4$individual <- as.character(df4$individual)
df4$worded_health <- as.character(df4$worded_health)
df4$sample_health_state <- as.character(df4$sample_health_state)
df4$species <- as.character(df4$species)

DRTO$individual <- as.character(DRTO$individual)
DRTO$worded_health <- as.character(DRTO$worded_health)
DRTO$sample_health_state <- as.character(DRTO$sample_health_state)
DRTO$species <- as.character(DRTO$species)


DRTO <- DRTO[rowSums(is.na(DRTO)) != ncol(DRTO), ]
View(DRTO)
merged_df <- merge(df4, DRTO, by = c("individual", "worded_health", "sample_health_state", "species"), all = TRUE)
View(merged_df)

##there are some differences in these datasets so we are making them identical
# Check which individuals in df4 are NOT in DRTO
setdiff(df4$individual, DRTO$individual)

# Check which individuals in DRTO are NOT in df4
setdiff(DRTO$individual, df4$individual)

# Repeat for other merge columns:
setdiff(df4$worded_health, DRTO$worded_health)
setdiff(df4$sample_health_state, DRTO$sample_health_state)
setdiff(df4$species, DRTO$species)

df4$individual <- tolower(df4$individual)
DRTO$individual <- tolower(DRTO$individual)

df4$worded_health <- tolower(df4$worded_health)
DRTO$worded_health <- tolower(DRTO$worded_health)

df4$worded_health <- gsub("[^[:alnum:] ]", "", df4$worded_health) 
DRTO$worded_health <- gsub("[^[:alnum:] ]", "", DRTO$worded_health)

df4$individual <- trimws(df4$individual)
DRTO$individual <- trimws(DRTO$individual)

df4$worded_health <- trimws(df4$worded_health)
DRTO$worded_health <- trimws(DRTO$worded_health)

df4$individual <- as.character(df4$individual)
DRTO$individual <- as.character(DRTO$individual)

df4$worded_health <- as.character(df4$worded_health)
DRTO$worded_health <- as.character(DRTO$worded_health)

df4$worded_health <- gsub("Hexposed", "H_exposed", df4$worded_health)
DRTO$worded_health <- gsub("Hexposed", "H_exposed", DRTO$worded_health)

merged_df <- merge(df4, DRTO, by = c("individual", "worded_health", "sample_health_state", "species"), all = TRUE)
View(merged_df)

# Convert categorical variables to factors
merged_df$worded_health <- as.factor(merged_df$worded_health)  
merged_df$species <- as.factor(merged_df$species)

# Select relevant columns for analysis
merged_subset <- merged_df[, c("worded_health", "species",
                        "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                        "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                        "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                        "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure",
                        "max_gastro", "max_max_vac", "max_prop_exo", "max_symb_vac", 
                        "min_min_symb", "mean_degr_symb")]

# Remove rows with missing values
merged_subset_clean <- na.omit(merged_subset)

# Train Random Forest model
set.seed(123)
mrf_model <- randomForest(worded_health ~ ., data =merged_subset_clean, importance = TRUE, ntree = 500)

# Print model summary
print(mrf_model)

# Plot feature importance
varImpPlot(mrf_model)

# Load the required library
install.packages("gbm")
library(gbm)

# Train a Gradient Boosting model
mgbm_model <- gbm(worded_health ~ ., data = merged_subset_clean, 
                 distribution = "multinomial", # For classification tasks
                 n.trees = 500, 
                 interaction.depth = 3,  # Depth of the trees
                 shrinkage = 0.1,        # Learning rate
                 cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(mgbm_model)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")


##loss of structure, blown out gastro, fungus/sponge, max_gastro and mean_degr_symb contribute less than 1
##does this change if instead of worded_health, we use state

merged_df$state <- as.factor(merged_df$state)
# Select relevant columns for analysis
merged_subset2 <- merged_df[, c("state", "species",
                               "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                               "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                               "gastro_sep_I", "degraded_symb_I", "fungus_sponge",
                               "blown_out_gastro", "amoebocytes", "loss_of_eosin", "loss_of_structure",
                               "max_gastro", "max_max_vac", "max_prop_exo", "max_symb_vac", 
                               "min_min_symb", "mean_degr_symb")]

# Remove rows with missing values
merged_subset2_clean <- na.omit(merged_subset2)

# Train Random Forest model
set.seed(123)
mrf_model2 <- randomForest(state ~ ., data =merged_subset2_clean, importance = TRUE, ntree = 500)

# Print model summary
print(mrf_model2)

# Plot feature importance
varImpPlot(rf_model)

# Train a Gradient Boosting model
mgbm_model2 <- gbm(state ~ ., data = merged_subset2_clean, 
                  distribution = "multinomial", # For classification tasks
                  n.trees = 500, 
                  interaction.depth = 3,  # Depth of the trees
                  shrinkage = 0.1,        # Learning rate
                  cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(mgbm_model2)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##naive and H-exposed can be grouped into one
##42% error so we are removing mean_degr_symb, max_gastro, Fungus_sponge, blown_out_gastro, loss_of_structure, and loss_of_eosin

# Select relevant columns for analysis
merged_subset3 <- merged_df[, c("state", "species",
                                "necrosis_C", "vacuolization_C", "exocytosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                                "gastro_sep_I", "degraded_symb_I",
                                "amoebocytes",
                                "max_max_vac", "max_prop_exo", "max_symb_vac", 
                                "min_min_symb")]

# Remove rows with missing values
merged_subset3_clean <- na.omit(merged_subset3)

# Train Random Forest model
set.seed(123)
mrf_model3 <- randomForest(state ~ ., data =merged_subset3_clean, importance = TRUE, ntree = 500)

# Print model summary
print(mrf_model3)

# Plot feature importance
varImpPlot(mrf_model3)

# Train a Gradient Boosting model
mgbm_model3 <- gbm(state ~ ., data = merged_subset3_clean, 
                   distribution = "multinomial", # For classification tasks
                   n.trees = 500, 
                   interaction.depth = 3,  # Depth of the trees
                   shrinkage = 0.1,        # Learning rate
                   cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(mgbm_model3)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##37% error so the first metric is still great. that's all
## removing amoebocytes, gastro sep intensity, degraded symb intensity, vacuolization consistency


# Select relevant columns for analysis
merged_subset4 <- merged_df[, c("state", "species",
                                "necrosis_C", "exocytosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                                "max_max_vac", "max_prop_exo", "max_symb_vac", 
                                "min_min_symb")]

# Remove rows with missing values
merged_subset4_clean <- na.omit(merged_subset4)

# Train Random Forest model
set.seed(123)
mrf_model4 <- randomForest(state ~ ., data =merged_subset4_clean, importance = TRUE, ntree = 500)

# Print model summary
print(mrf_model4)

# Plot feature importance
varImpPlot(mrf_model4)

# Train a Gradient Boosting model
mgbm_model4 <- gbm(state ~ ., data = merged_subset4_clean, 
                   distribution = "multinomial", # For classification tasks
                   n.trees = 500, 
                   interaction.depth = 3,  # Depth of the trees
                   shrinkage = 0.1,        # Learning rate
                   cv.folds = 5)           # Cross-validation

# View the summary of the model
summary(mgbm_model4)

# Make predictions
predictions <- predict(gbm_model, newdata = DRTO_subset_clean, n.trees = 500, type = "response")

##moved to 38% error
##I think this isn't representing the data well enough now that some are PA and some are measurements/proportions
# Install the e1071 package if not already installed
install.packages("e1071")

# Load the package
library(e1071)
svm_model <- svm(state ~ ., data = merged_subset3_clean)

# Make predictions
predictions <- predict(svm_model, merged_subset3_clean)

# Evaluate the model
table(predictions, merged_subset3_clean$state)

# Calculate Accuracy
accuracy <- sum(diag(table(predictions, merged_subset3_clean$state))) / sum(table(predictions, merged_subset3_clean$state))
print(paste("Accuracy: ", accuracy))

# Precision, Recall, and F1 for each class
install.packages("caret")
library(caret)
conf_matrix <- confusionMatrix(predictions, merged_subset3_clean$state)

# Print Precision, Recall, and F1 Score for each class
print(conf_matrix)

##lets remove some things and see if its better or worse
##making a new subset

merged_subset4 <- merged_df[, c("state", "species",
                                "necrosis_C", "exocytosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "necrosis_I", "vacuolization_I", "exocytosis_I", 
                                "max_max_vac", "max_prop_exo", "max_symb_vac", 
                                "min_min_symb")]
merged_subset4_clean <- na.omit(merged_subset4)
svm_model2 <- svm(state ~ ., data = merged_subset4_clean)

# Make predictions
predictions2 <- predict(svm_model2, merged_subset4_clean)

# Evaluate the model
table(predictions2, merged_subset4_clean$state)

# Calculate Accuracy
accuracy2 <- sum(diag(table(predictions2, merged_subset4_clean$state))) / sum(table(predictions2, merged_subset4_clean$state))
print(paste("Accuracy: ", accuracy2))

# Precision, Recall, and F1 for each class
conf_matrix2 <- confusionMatrix(predictions2, merged_subset4_clean$state)

# Print Precision, Recall, and F1 Score for each class
print(conf_matrix2)

##lets remove some things and see if its better or worse
##when it was a random forest the lowest values were exocytosis consistency and intensity...
##making a new subset

merged_subset5 <- merged_df[, c("state", "species",
                                "necrosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "necrosis_I", "vacuolization_I", 
                                "max_max_vac", "max_prop_exo", "max_symb_vac", 
                                "min_min_symb")]
merged_subset5_clean <- na.omit(merged_subset5)
svm_model3 <- svm(state ~ ., data = merged_subset5_clean)

# Make predictions
predictions3 <- predict(svm_model3, merged_subset5_clean)

# Evaluate the model
table(predictions3, merged_subset5_clean$state)

# Calculate Accuracy
accuracy3 <- sum(diag(table(predictions3, merged_subset5_clean$state))) / sum(table(predictions3, merged_subset5_clean$state))
print(paste("Accuracy: ", accuracy3))

# Precision, Recall, and F1 for each class
conf_matrix3 <- confusionMatrix(predictions3, merged_subset5_clean$state)

# Print Precision, Recall, and F1 Score for each class
print(conf_matrix3)

# Make sure you're passing the formula used to train the model
plot(svm_model2, data = merged_subset4_clean)

plot(svm_model2, formula = state ~ ., data = merged_subset4_clean)

##that made the model less so lets run it back
##does necrosis matter
##making a new subset

merged_subset6 <- merged_df[, c("state", "species",
                                "exocytosis_C", "gastro_sep_C", 
                                "degraded_symb_C", "vacuolization_I", "exocytosis_I", 
                                "max_max_vac", "max_prop_exo", "max_symb_vac", 
                                "min_min_symb")]
merged_subset6_clean <- na.omit(merged_subset6)
svm_model4 <- svm(state ~ ., data = merged_subset6_clean)

# Make predictions
predictions4 <- predict(svm_model4, merged_subset6_clean)

# Evaluate the model
table(predictions4, merged_subset6_clean$state)

# Calculate Accuracy
accuracy4 <- sum(diag(table(predictions4, merged_subset6_clean$state))) / sum(table(predictions4, merged_subset6_clean$state))
print(paste("Accuracy: ", accuracy4))

# Precision, Recall, and F1 for each class
conf_matrix4 <- confusionMatrix(predictions4, merged_subset6_clean$state)

# Print Precision, Recall, and F1 Score for each class
print(conf_matrix4)

##lets give up on improving this model since we've clearly gotten it as good as its gonna get

# Select only numeric columns for the heatmap
numeric_data <- merged_subset4_clean[, sapply(merged_subset4_clean, is.numeric)]

# Create a heatmap for the numeric data
heatmap(as.matrix(numeric_data), scale = "column", col = colorRampPalette(c("blue", "white", "red"))(100), margins = c(5,10))


##lets make that better
# Ensure 'state' is a factor and order the dataset by state
merged_subset4_clean$state <- as.factor(merged_subset4_clean$state)
ordered_data <- merged_subset4_clean[order(merged_subset4_clean$state), ]
# Select only numeric columns
numeric_data <- ordered_data[, sapply(ordered_data, is.numeric)]

# Convert to matrix for heatmap
numeric_matrix <- as.matrix(numeric_data)

# Assign row names as state labels for grouping in the heatmap
rownames(numeric_matrix) <- ordered_data$state
# Load necessary package
install.packages("pheatmap")
library(pheatmap)

# Generate heatmap, clustering only within each 'state' category
pheatmap(numeric_matrix, 
         scale = "column",  # Normalize across columns for better contrast
         annotation_row = data.frame(State = ordered_data$state),  # Add grouping
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Define color scale
         cluster_rows = FALSE,  # Prevents clustering so it's grouped by state
         cluster_cols = TRUE)  # Keeps clustering for features

# Create an annotation data frame
annotation_data <- data.frame(State = ordered_data$state)

# Ensure row names in annotation data match the heatmap data
rownames(annotation_data) <- rownames(numeric_matrix)

# Generate heatmap with annotations
pheatmap(numeric_matrix, 
         scale = "column",
         annotation_row = annotation_data,  # Ensure correct row names
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = FALSE,  # Keep data grouped by state
         cluster_cols = TRUE)
library(pheatmap)

# Ensure 'state' is a factor and order the data by state
ordered_data <- merged_subset4_clean[order(merged_subset4_clean$state), ]
ordered_data$state <- as.factor(ordered_data$state)

# Select only numeric features
numeric_data <- ordered_data[, sapply(ordered_data, is.numeric)]
numeric_matrix <- as.matrix(numeric_data)

# Create annotation data as a column (without setting row names)
annotation_data <- data.frame(State = ordered_data$state)  

# Generate heatmap without row name duplication
pheatmap(numeric_matrix, 
         scale = "column",
         annotation_row = annotation_data,  # Correctly added without setting row names
         color = colorRampPalette(c("red", "white", "blue"))(100),
         cluster_rows = FALSE, 
         cluster_cols = TRUE)
rownames(annotation_data) <- rownames(numeric_matrix)

# Define custom colors for the 'state' annotation
state_colors <- list(State = c("DD" = "darkorange2", "HD" = "goldenrod1", "HH" = "darkolivegreen3"))
pheatmap(numeric_matrix, 
         scale = "column",
         annotation_row = annotation_data,
         annotation_colors = state_colors,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 6,  # Adjust row font size
         height = 10)  # Increase figure height
## Here are the hex colors for this one:
##"#5a7bac","#dee6f6","#bd004a"

##And the full palette:
##  "#5a7bac", "#acc5de", "#dee6f6","#bd004a","#ffcde6"

##And the link to where I generated it:
##  https://pokepalettes.com/#corsola
