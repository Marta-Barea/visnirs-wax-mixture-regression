###############################################################################
######## Machine learning-based approaches on Vis-NIR data for the ############
################### quantification of petroleum wax blends ####################
######################## Marta Barea-Sepúlveda ################################
###############################################################################

###############################################################################
######################## Partial Least Squares (PLS) ##########################
###############################################################################

# Loading Packages

library(readxl)
library(caret)
library(prospectr)
library(MLmetrics)
library(data.table)
library(doParallel)
library(stringr)
library(ggplot2)

# Loading Parallelization

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# Loading data

pw_data <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_PW_Type_Hydrotreating.xlsx", 
                      sheet = "PW_Type_Mix")

ext_val <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_PW_Type_Hydrotreating.xlsx", 
                      sheet = "External_Cal")

# Savitzky Golay Smoothing

sgvec <- savitzkyGolay(X = pw_data[,-c(1,2)], p = 3, w = 11, m = 1)
pw_sg <- cbind.data.frame(Sample = pw_data$ID_Reg, sgvec)

pw_sg$Sample <- as.numeric(pw_sg$Sample)

sgvec_1 <- savitzkyGolay(X = ext_val[,-c(1,2)], p = 3, w = 11, m = 1)
pw_sg_ext <- cbind.data.frame(Sample = ext_val$ID_Reg, sgvec_1)

pw_sg_ext$Sample <- as.numeric(pw_sg_ext$Sample)

# Data slicing

set.seed(1345)

intrain <- createDataPartition(y = pw_sg$Sample, 
                               p = 0.7, 
                               list = FALSE)
pw_train <- pw_sg[intrain,]
pw_test <- pw_sg[-intrain,]

# Hyperparameter tuning and model training

set.seed(1345)

trctrl <- trainControl(method = "cv", 
                       number = 5,
                       returnResamp = "final", 
                       verboseIter = FALSE, 
                       savePredictions = "final")

start_time <- Sys.time()

pls_model <- train(Sample ~., 
                   data = pw_train,
                   trControl = trctrl,
                   method = "pls", 
                   metric = "RMSE",
                   tuneLength = 30) 
pls_model

total_time <- Sys.time() - start_time
total_time

summary(pls_model$finalModel)

ncomp_plot <- plot(pls_model,
                   xlab = "Number of components",
                   ylab = "RMSE (5-Fold CV)",
                   type = "b",
                   pch = 19,
                   col = "#EA6A60")
ncomp_plot

# Train set predictions

train_pred <- predict(pls_model, newdata = pw_train[,-1])

mse_train = MSE(pw_train$Sample, train_pred)
mae_train = MAE(pw_train$Sample, train_pred)
rmse_train = RMSE(pw_train$Sample, train_pred)
r2_train = caret::R2(pw_train$Sample, train_pred)

cat("MAE:", mae_train, "\n", "MSE:", mse_train, "\n", 
    "RMSE:", rmse_train, "\n", "R-squared:", r2_train)

# Test set predictions

test_pred <- predict(pls_model, newdata = pw_test[,-1])

mse_test = MSE(pw_test$Sample, test_pred)
mae_test = MAE(pw_test$Sample, test_pred)
rmse_test = RMSE(pw_test$Sample, test_pred)
r2_test = caret::R2(pw_test$Sample, test_pred)

cat("MAE:", mae_test, "\n", "MSE:", mse_test, "\n", 
    "RMSE:", rmse_test, "\n", "R-squared:", r2_test)

# External validation predictions

ext_pred <- predict(pls_model, newdata = pw_sg_ext[,-1])

mse_ext = MSE(pw_sg_ext$Sample, ext_pred)
mae_ext = MAE(pw_sg_ext$Sample, ext_pred)
rmse_ext = RMSE(pw_sg_ext$Sample, ext_pred)
r2_ext = caret::R2(pw_sg_ext$Sample, ext_pred)

cat("MAE:", mae_ext, "\n", "MSE:", mse_ext, "\n", 
    "RMSE:", rmse_ext, "\n", "R-squared:", r2_ext)

# Predictions' plot

train_pred_matrix <- data.frame(Real = c(pw_train$Sample),
                                Predicted = c(train_pred))

test_pred_matrix <- data.frame(Real = c(pw_test$Sample),
                               Predicted = c(test_pred))

ext_pred_matrix <- data.frame(Real = c(pw_sg_ext$Sample),
                              Predicted = c(ext_pred))

replace_label <- function(mixture) {
  if (mixture == 0) {
    return("Micro_Wax")
  }
  if (mixture == 100) {
    return("Macro_Wax")
  }
  
  return(paste("Mix_", toString(mixture), "%", sep = ""))
}

train_pred_matrix$Sample <- lapply(train_pred_matrix$Real, replace_label)
train_group_labels <- rep("Train set", times = 52)
train_pred_matrix <- cbind(train_pred_matrix, Group = train_group_labels)

test_pred_matrix$Sample <- lapply(test_pred_matrix$Real, replace_label)
test_group_labels <- rep("Test set", times = 20)
test_pred_matrix <- cbind(test_pred_matrix, Group = test_group_labels)

ext_pred_matrix$Sample <- lapply(ext_pred_matrix$Real, replace_label)
ext_group_labels <- rep("External validation", times = 8)
ext_pred_matrix <- cbind(ext_pred_matrix, Group = ext_group_labels)

df <- rbind(train_pred_matrix, test_pred_matrix, ext_pred_matrix)

remove_col<- c("Sample","Group")
labels_id <- as.matrix(df[, !(names(df) %in% remove_col)])
rownames(labels_id) <- df$Sample

group <- as.factor(df$Group)

plot(x = df$Real,
     y = df$Predicted,
     type = "p",
     pch = c(8, 7, 16)[factor(group)],
     col = "black",
     main = "",
     xlab = "Real (%)",
     ylab = "Predicted (%)",
     xlim = c(0,100))
lines(0:100, 0:100, lwd = 1, lty = 1, col = "#EA6A60")
legend(x = 5, y = 95, c("External validation set", "Test set", "Training set"), cex = 0.8, pch = c(8, 7, 16))

# Variable Importance

plot(varImp(object = pls_model), 
     top = 20,
     ylab = "Feature")

var_imp <- varImp(object = pls_model)
var_imp_1 <- var_imp[['importance']]

names <- rownames(var_imp_1)
rownames(var_imp_1) <- NULL

var_imp_2 <- cbind(names,var_imp_1)

sp_feat <- as.data.frame(var_imp_2)
sp_feat_1 <- sp_feat[,-2]

spfeatWithoutQuotes = c()

for (onespfeat in sp_feat_1) {
  onespfeat = str_replace_all(onespfeat, "`", "")
  onespfeat = str_replace_all(onespfeat, r"(\\)", "")
  spfeatWithoutQuotes = append(spfeatWithoutQuotes, onespfeat)
}

print(spfeatWithoutQuotes)

var_imp_3 <- cbind(var_imp_2, Variable = spfeatWithoutQuotes)
var_imp_3 <- var_imp_3[,-1]

var_imp_plot <- ggplot(var_imp_3, aes(x = Variable, y = Overall)) +
  geom_bar(stat = "identity", fill = "#5FA8F1") +
  labs(title = "", x = ("Wavelength (nm)"), y = "Relative Importance (%)") +
  theme_test() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 8, hjust = 1, angle = 90),
        axis.title = element_text(size = 10)) +
  scale_x_discrete(limits = var_imp_3$Variable,
                   breaks = var_imp_3$Variable[seq(1, length(var_imp_3$Variable), by = 50)]) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "blue")

var_imp_plot

# Stop Parallelization

stopCluster(cl)