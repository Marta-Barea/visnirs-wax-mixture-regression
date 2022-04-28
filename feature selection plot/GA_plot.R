###############################################################################
######## Machine learning-based approaches on Vis-NIR data for the ############
################### quantification of petroleum wax blends ####################
########################## Marta Barea-Sepúlveda ##############################
###############################################################################

###############################################################################
########################## Genetic algorithm plot #############################
###############################################################################

# Loading Packages

library(doParallel)
library(readxl)
library(prospectr)
library(caret)
library(Boruta)
library(dplyr)
library(stringr)
library(ggplot2)

# Loading Parallelization

cl <- makePSOCKcluster(8)
registerDoParallel(cl)

# Loading data
pw_data <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_PW_Type_Hydrotreating.xlsx", 
                      sheet = "PW_Type_Mix")

# Savitzky Golay Smoothing

sgvec <- savitzkyGolay(X = pw_data[,-c(1,2)], p = 3, w = 11, m = 1)
pw_sg <- cbind.data.frame(Sample = pw_data$ID_Reg, sgvec)

pw_sg$Sample <- as.numeric(pw_sg$Sample)

# Data slicing

set.seed(1345)

intrain <- createDataPartition(y = pw_sg$Sample, 
                               p = 0.7, 
                               list = FALSE)
pw_train <- pw_sg[intrain,]
pw_test <- pw_sg[-intrain,]

# Feature selection by means of the GA

set.seed(1345)

ga_ctrl <- gafsControl(functions = rfGA,
                       method = "cv",
                       number = 5)

rf_ga <- gafs(x = pw_train %>% select(-Sample), 
              y = pw_train$Sample,
              popSize = 10,
              iters = 10,
              gafsControl = ga_ctrl, 
              metric = "RMSE",
              ntree = 100)

rf_ga

ga_feat <- rf_ga[["optVariables"]]

# First derivative spectra along with the selected features

df_2 <- reshape2::melt(pw_sg, "Sample")

ga_plot <- ggplot(data = df_2, aes(x = variable, y = value, color = Sample)) + 
  geom_line() +
  labs(x = "Wavelength (nm)", y = "Absorbance") +
  theme_test() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, hjust = 1, angle = 90),
        axis.title = element_text(size = 8)) +
  scale_x_discrete(limits = df_2$variable,
                   breaks = df_2$variable[seq(1, length(df_2$variable), by = 3000)])+
  geom_vline(xintercept = ga_feat, linetype = 1, colour = "#4286BE", size = 0.25, alpha = 0.20)

ga_plot

# Stop Parallelization

stopCluster(cl)