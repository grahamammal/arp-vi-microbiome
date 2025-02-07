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

mtry_seq <- floor(sqrt(ncol(.rslurm_more_args$x_16s))* c(0.5, 1, 2))

ranger_learners <- create.Learner("SL.ranger", tune = list(mtry = mtry_seq, max.depth = c(0, 3, 5)))
xgboost_learners <- create.Learner("SL.xgboost", tune = list(max_depth = c(1, 3, 5), shrinkage = c(0.001, 0.01, 0.1)))

glmnet_learners <- create.Learner("SL.glmnet", tune = list(alpha = c(0, 1)))

corP_screens <- map2(c("Sl.ranger", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names), "screen.corP", c)
lasso_screens <- map2(c("Sl.ranger", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names), "screen.glmnet", c)
random_forest_screens <- map2(c("Sl.ranger", ranger_learners$names, xgboost_learners$names, "SL.glm", glmnet_learners$names), "screen.randomForest", c)

my_libraries <- c(corP_screens, lasso_screens, random_forest_screens)

cluster <- parallel::makeCluster(length(.rslurm_x))
parallel::clusterExport(cluster, c(ranger_learners$names, xgboost_learners$names, glmnet_learners$names))

.rslurm_result <- do.call(parallel::mclapply, c(list(
    X = .rslurm_x[.rslurm_istart:.rslurm_iend],
    FUN = .rslurm_func),
    .rslurm_more_args,
    mc.cores = 1,
    mc.preschedule = TRUE
    ))

saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
parallel::stopCluster()
