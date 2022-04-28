###############################################################################
######## Machine learning-based approaches on Vis-NIR data for the ############
################### quantification of petroleum wax blends ####################
######################## Marta Barea-Sepúlveda ################################
###############################################################################

###############################################################################
###################### Support Vector Regression (SVR) ########################
######################### Radial Basis Function (RBF) #########################
###############################################################################

# Loading Packages

library(readxl)
library(caret)
library(prospectr)
library(data.table)
library(MLmetrics)
library(dplyr)
library(graphics)
library(doParallel)

# Loading Parallelization

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# Loading data

pw_data <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_PW_Type_Hydrotreating.xlsx", 
                      sheet = "PW_Type_Mix")

ext_val <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_PW_Type_Hydrotreating.xlsx", 
                      sheet = "External_Cal")

# Savitzky Golay Smoothing

sgvec <- savitzkyGolay(X = pw_data[, -c(1,2)], p = 3, w = 11, m = 1)
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

# Hyperparameters tuning and model training

set.seed(1345)

trctrl_1 <- trainControl(method = "cv", number = 5)
                    
grid_radial_1 <- expand.grid(sigma = c(2^(seq(-10, 10, 0.5))), 
                             C = c(2^(seq(-10, 10, 0.5))))

start_time <- Sys.time()

svr_radial <- train(Sample ~., 
                    data = pw_train,
                    method = "svmRadial",
                    trControl= trctrl_1,
                    tuneGrid = grid_radial_1,
                    metric = "RMSE")

svr_radial

total_time <- Sys.time() - start_time
total_time

svr_radial$finalModel

subset(svr_radial$results, 
       svr_radial$results$C==svr_radial[["bestTune"]][["C"]] & 
         svr_radial$results$sigma==svr_radial[["bestTune"]][["sigma"]])

# Final SVR model

set.seed(1345)

trctrl_2 <- trainControl(method = "cv", number = 5)
grid_radial_2 <- expand.grid(sigma = svr_radial[["bestTune"]][["sigma"]], 
                             C = svr_radial[["bestTune"]][["C"]])

best_svr <- train(Sample ~., 
                  data = pw_train,
                  method = "svmRadial",
                  trControl= trctrl_2,
                  tuneGrid = grid_radial_2,
                  metric = "RMSE")
best_svr

# Train set predictions

train_pred <- predict(best_svr, newdata = pw_train[,-1])

mse_train = MSE(pw_train$Sample, train_pred)
mae_train = MAE(pw_train$Sample, train_pred)
rmse_train = RMSE(pw_train$Sample, train_pred)
r2_train = caret::R2(pw_train$Sample, train_pred)

cat("MAE:", mae_train, "\n", "MSE:", mse_train, "\n", 
    "RMSE:", rmse_train, "\n", "R-squared:", r2_train)

# Test set predictions

test_pred <- predict(best_svr, newdata = pw_test[,-1])

mse_test = MSE(pw_test$Sample, test_pred)
mae_test = MAE(pw_test$Sample, test_pred)
rmse_test = RMSE(pw_test$Sample, test_pred)
r2_test = caret::R2(pw_test$Sample, test_pred)

cat("MAE:", mae_test, "\n", "MSE:", mse_test, "\n", 
    "RMSE:", rmse_test, "\n", "R-squared:", r2_test)

# External validation predictions

ext_pred <- predict(best_svr, newdata = pw_sg_ext[,-1])

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


# Contour plot of the SVR model

rmse_gridsearch <- matrix((svr_radial[["results"]][["RMSE"]]), ncol = 41, nrow = 41)

cost_expression <- expression(log[2] ~ C)
gamma_expression <- expression(log[2] ~ γ)

filled.contour(x = c(seq(-10,10,length.out = 41)),
               y = c(seq(-10,10,length.out = 41)),
               z = as.matrix(rmse_gridsearch),
               color.palette = colorRampPalette(c("red", "orange", "purple", "blue")), 
               plot.title = title(main = "", sub = "", xlab = cost_expression, ylab = gamma_expression),
               plot.axes = {axis(1,seq(-10, 10, 1),cex.axis = 1,las = 2)
                 axis(2,seq(-10, 10, 1),cex.axis = 1,las = 2)})

# Stop Parallelization

stopCluster(cl)