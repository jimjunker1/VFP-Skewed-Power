# estimating counts in each body size
# source("back_estimating_counts.R")

library(sizeSpectra)
library(tidyverse)
library(poweRlaw)
library(parallel)
library(foreach)
library(doParallel)


# source the custom functions written for this simulation study
source("custom_functions.R")
# source variable designations, including setting seed
source("master_variable_designation.R")


rep_df <- expand_grid(
  pr_scenarios,
  n = n,
  lambda = lambda,
  xmin = xmin,
  xmax = xmax,
  vecDiff = vecDiff)

# for working out issues
# rep_df <- rep_df |>
#   sample_n(1)
# 
# n_iter <- 7

# expand df so each parameter set has one row per rep
rep_df <- expand_grid(
  rep_df, 
  rep = 1:n_iter
)

#cores <- detectCores()-1 # when running on its own
cores <- 7 # when running with other simulations
cluster <- makeCluster(cores)
registerDoParallel(cluster)

# old dopar
# error throws out whole rep?
# tictoc::tic()
# results <- foreach(j = 1:n_iter, 
#                    .combine = rbind,
#                    .errorhandling = "remove",
#                    .packages = c("sizeSpectra",
#                                  "tidyverse",
#                                  "poweRlaw")
#                    ) %dopar%{
#   df_out <- parallel_rep_compare_pareto_n(rep_df)
#   df_out$rep <- j
#   #results[[j]] <- df_out
#   df_out}
# 
# results
# tictoc::toc()

# new dopar solution
tictoc::tic()
#set.seed(1130)
results <- foreach(j = 1:nrow(rep_df), 
                   .combine = rbind,
                   .errorhandling = "remove",
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw")) %dopar%{
    set.seed(j+1130)
    df_out <- compare_est_pareto_n(
      n = rep_df[j,]$n, 
      lambda = rep_df[j,]$lambda,
      xmin = rep_df[j,]$xmin, 
      xmax = rep_df[j,]$xmax, 
      h = rep_df[j,]$h, 
      b = rep_df[j,]$b,
      vecDiff = rep_df[j,]$vecDiff)
                                       
    df_out$pr_scenario <- rep_df[j,]$pr_scenario
    df_out$rep <- rep_df[j,]$rep
    df_out
                                 }
#results
tictoc::toc()

stopCluster(cl = cluster)

#results
saveRDS(results, "simulation_results/est_pareto_N_sim_run_test.rds")


