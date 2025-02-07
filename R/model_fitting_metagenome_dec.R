

# load required libraries
library(tidyverse)
library(SuperLearner)
library(xgboost)
library(ranger)
library(glmnet)
library(vimp)
library(randomForest)
library(rslurm)
library(caret)
library(mboost)

# load data

combined_data_metagenome <- read_csv("data/formatted_data/combined_data_metagenome.csv")
variable_sets_metagenome <- read_rds("data/variable_sets_metagenome.rds")

#------------Preparing Superlearner -------------

# neff = 358
# Use AUC for performance metric
# Use 20 Fold Cross Validation (Practical considerations paper), stratify on Y when assigning units to validation sets
# Use screening, lasso, random forest
# Use up to 70 learners, include parametric learners




# set up data
y_metagenome <- as.numeric(combined_data_metagenome$dec_infection == "Infected")
x_metagenome <- select(combined_data_metagenome, -diarrhea, -sample_id, -alternate_metagenome_id, -dec_infection)  %>% 
  mutate(diversity =(diversity - mean(diversity))/sd(diversity))

SL.caretGAM <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'gamSpline',
           tuneLength = 50,
           trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                            verboseIter = TRUE),
           ...)
}

SL.caretGLMB <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'glmboost',
           tuneLength = 50,
           trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                            verboseIter = TRUE),
           ...)
}

ranger_learners <- create.Learner("SL.ranger", tune = list(max.depth = c(0, 3, 5)))
xgboost_learners <- create.Learner("SL.xgboost", tune = list(shrinkage = c(0.01, 0.1)))


glmnet_learners <- create.Learner("SL.glmnet", tune = list(alpha = c(0, 1)))


# --------------Running VIMP + Superlearner ------------------
run_vimp_holdout_set <- function(indx, y_metagenome, x_metagenome, ...) {
  # set up learners
  
  
  ranger_learners <- create.Learner("SL.ranger", tune = list(max.depth = c(0, 3, 5)))
  xgboost_learners <- create.Learner("SL.xgboost", tune = list(shrinkage = c(0.01, 0.1)))
  
  SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 3,
             trControl =  caret::trainControl(method = "LGOCV", number = 1,
                                              verboseIter = TRUE), ...)
  }
  SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree', tuneLength = 3,
             trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                              verboseIter = TRUE), ...)
  }
  
  SL.caretGAM <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'gamSpline',
             tuneLength = 3,
             trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                              verboseIter = TRUE),
             ...)
  }
  
  SL.caretGLMB <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'glmboost',
             tuneLength = 3,
             trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                              verboseIter = TRUE),
             ...)
  }
  
  glmnet_learners <- create.Learner("SL.glmnet", tune = list(alpha = c(0, 1)))
  
  lasso_screens <- map2(c("SL.ranger", "SL.gam", ranger_learners$names, "SL.glm", glmnet_learners$names), "screen.glmnet", c)
  random_forest_screens <- map2(c("SL.ranger", "SL.gam", ranger_learners$names, "SL.glm", glmnet_learners$names), "screen.randomForest", c)
  
  my_libraries <- c(lasso_screens, random_forest_screens)
  
  
  sl_cvcontrol <- list(V = 30, stratifyCV = TRUE)
  
  set.seed(522305)
  vimp_output <- vimp_auc(y_metagenome, x_metagenome, SL.library = c("SL.mean", my_libraries), 
                          V = 10, indx = indx, cvControl = sl_cvcontrol, ...)
  vimp_output
}

slurm_map(variable_sets_metagenome, run_vimp_holdout_set, 
          jobname = "run_vimp", nodes = 100, cpus_per_node = 1, submit = FALSE, 
          libPaths = ("/home/users/graham18/Desktop/R_lib"),
          y_metagenome = y_metagenome, x_metagenome = x_metagenome,
          method = "method.AUC")

# # 
# start_time <- Sys.time()
# result_9 <- run_vimp_holdout_set(variable_sets_metagenome[[10]], y_metagenome, x_metagenome)
# end_time <- Sys.time()
# 
# end_time - start_time
