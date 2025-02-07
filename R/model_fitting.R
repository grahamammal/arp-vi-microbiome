# Setup -------------------------------------------

# Useful calls for server (right click to paste):
# Run simulations:
# qsub -cwd -e /dev/null -o /dev/null -v sim_run='first' run_sim_correlated.sh
# qsub -cwd -e /dev/null -o /dev/null -v sim_run='main' -t 1-1000 -hold_jid 111111 run_sim_correlated.sh
# qsub -cwd -e /dev/null -o /dev/null -v sim_run='last' -hold_jid 111111 run_sim_correlated.sh
# can add -q students.q to run on student queue
# Delete error and output files:
# find . -name 'run_sim.sh.[eo]*' -exec rm {} \;
# --
# To run on the slurm partition:
# Two partitions (queues available)
# "12c-128g-all" (default, for all users)
# "12c-128g-students" (for student use only)
# sbatch --export=sim_run='first' run_sim_correlated_slurm.sh --output=/dev/null
# sbatch --export=sim_run='main' --array=1-20 --depend=afterok:101 run_sim_correlated_slurm.sh --output=/dev/null
# sbatch --export=sim_run='last' --depend=afterok:102 run_sim_correlated_slurm.sh --output=/dev/null
# If error:
# sbatch: error: Batch script contains DOS line breaks (\r\n)
# sbatch: error: instead of expected UNIX line breaks (\n).
# Then run the following:
# sed -i 's/^M//' script_name.sh
# The character ^M is a single special character. To type it press and hold CTRL. Then Press and release v, then still holding CTRL press m



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

# load data

combined_data_16s <- read_csv("data/formatted_data/combined_data_16s.csv")
family_variables <- read_rds("data/family_variables.rds")

#------------Preparing Superlearner -------------

# neff = 358
# Use AUC for performance metric
# Use 20 Fold Cross Validation (Practical considerations paper), stratify on Y when assigning units to validation sets
# Use screening, lasso, random forest
# Use up to 70 learners, include parametric learners




# set up data
y_16s <- as.numeric(combined_data_16s$diarrhea)
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


# --------- Variables to run VIMP on --------------
variable_positions <- data.frame(variable_name = names(x_16s), position = 1:length(names(x_16s)))

species_diversity_variables <- variable_positions %>% 
  filter(variable_name %in% c("species_richness", "shannon_diversity", "simpson_diversity")) %>% 
  pull(position)

index_sets <- family_variables
index_sets[[length(index_sets) + 1]] <- species_diversity_variables
names(index_sets)[length(index_sets)] <- "species_diversity_variables"
# all taxa
taxa_variables <- variable_positions %>% 
  filter(!(variable_name %in% c("species_richness", "shannon_diversity", "simpson_diversity"))) %>% 
  pull(position)
index_sets[[length(index_sets) + 1]] <- taxa_variables
names(index_sets)[length(index_sets)] <- "taxa_variables"
# --------------Running VIMP + Superlearner ------------------


run_vimp_holdout_set <- function(indx, y_16s, x_16s) {
  # set up learners
  

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
  
  sl_cvcontrol <- list(V = 25)
  
  set.seed(522305)
  vimp_output <- vimp_auc(y_16s, x_16s, SL.library = c("SL.mean", my_libraries), 
                          V = 10, indx = indx, cvControl = sl_cvcontrol)
  vimp_output
}

slurm_map(index_sets, run_vimp_holdout_set, 
          jobname = "run_vimp", nodes = 70, cpus_per_node = 1, submit = FALSE, 
          libPaths = ("/home/users/graham18/Desktop/R_lib"),
          y_16s = y_16s, x_16s = x_16s)
