


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

combined_data_16s <- read_csv("data/formatted_data/combined_data_16s.csv")
variable_sets_16s <- read_rds("data/variable_sets_16s.rds")

#------------Preparing Superlearner -------------

# neff = 358
# Use AUC for performance metric
# Use 20 Fold Cross Validation (Practical considerations paper), stratify on Y when assigning units to validation sets
# Use screening, lasso, random forest
# Use up to 70 learners, include parametric learners




# set up data
y_16s <- as.numeric(combined_data_16s$dec_infection == "Infected")
x_16s <- select(combined_data_16s, -diarrhea, -sample_id,
                -c(alternate_metagenome_id, 
                   sequenced_16s, 
                   dec_isolate_whole_genome_sequenced, 
                   dec_isolate_id, 
                   shotgun_metagenome_sequenced_from_whole_stool, 
                   dec_infection, 
                   dec_pathotype, 
                   participant_age_years, 
                   geometric_mean))



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

lasso_screens <- map2(c("SL.ranger", "SL.caretGLMB", "SL.caretRF", "SL.gam", ranger_learners$names, "SL.glm", glmnet_learners$names), "screen.glmnet", c)
random_forest_screens <- map2(c("SL.ranger", "SL.caretGLMB", "SL.caretRF", "SL.gam", ranger_learners$names, "SL.glm", glmnet_learners$names), "screen.randomForest", c)

my_libraries <- c(lasso_screens, random_forest_screens)

sl_cvcontrol <- list(V = 25)

set.seed(522305)

test_learner <- SuperLearner(Y = y_16s, X = x_16s, family = binomial(),
             SL.library = c("SL.mean"),
             method = "method.AUC",
             cvControl =sl_cvcontrol)
