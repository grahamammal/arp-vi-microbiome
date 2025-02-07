.libPaths(c(.libPaths(), '/home/users/graham18/Desktop/R_lib'))
library(base, quietly = TRUE)
library(methods, quietly = TRUE)
library(datasets, quietly = TRUE)
library(utils, quietly = TRUE)
library(grDevices, quietly = TRUE)
library(graphics, quietly = TRUE)
library(stats, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tibble, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(forcats, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(nnls, quietly = TRUE)
library(splines, quietly = TRUE)
library(foreach, quietly = TRUE)
library(gam, quietly = TRUE)
library(SuperLearner, quietly = TRUE)
library(xgboost, quietly = TRUE)
library(ranger, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(glmnet, quietly = TRUE)
library(vimp, quietly = TRUE)
library(rslurm, quietly = TRUE)
library(cvAUC, quietly = TRUE)
library(randomForest, quietly = TRUE)



.rslurm_func <- readRDS('f.RDS')
.rslurm_x <- readRDS('x.RDS')
.rslurm_more_args <- readRDS('more_args.RDS')
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
.rslurm_istart <- .rslurm_id * 1 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 1, length(.rslurm_x))

ranger_learners <- create.Learner("SL.ranger", tune = list(max.depth = c(0, 3, 5)))
xgboost_learners <- create.Learner("SL.xgboost", tune = list(shrinkage = c(0.01, 0.1)))

SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 50,
           trControl =  caret::trainControl(method = "LGOCV", number = 1,
                                            verboseIter = TRUE), ...)
}
SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree', tuneLength = 300,
           trControl =  caret::trainControl(method = "LGOCV", number = 1, search = 'random',
                                            verboseIter = TRUE), ...)
}

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

glmnet_learners <- create.Learner("SL.glmnet", tune = list(alpha = c(0, 1)))

lasso_screens <- map2(c("SL.ranger", "SL.caretGLMB", "SL.caretRF", "SL.caretGAM", "SL.caretXGB", "SL.gam", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names), "screen.glmnet", c)
random_forest_screens <- map2(c("SL.ranger", "SL.caretGLMB", "SL.caretRF", "SL.caretGAM", "SL.caretXGB", "SL.gam", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names), "screen.randomForest", c)

my_libraries <- c(lasso_screens, random_forest_screens)

cluster <- parallel::makeCluster(length(.rslurm_x))
parallel::clusterExport(cluster, c("SL.ranger", "SL.caretGLMB", "SL.caretRF", "SL.caretGAM", "SL.caretXGB", "SL.gam", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names))

.rslurm_result <- do.call(parallel::mclapply, c(list(
    X = .rslurm_x[.rslurm_istart:.rslurm_iend],
    FUN = .rslurm_func),
    .rslurm_more_args,
    mc.cores = 1,
    mc.preschedule = TRUE
    ))

saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
parallel::stopCluster()
