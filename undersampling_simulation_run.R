# run simulation with undersampling probability

# source("undersampling_simulation_run.R")

# This script uses a logistic (type 3) response for sampling probability
# Every body size above ~ 0.095 has a >90% chance of being sampled
# This functional response is similar to that found empirically by Morin et al. 2004  
# the variable values in the sample probabillity function were manually adjusted to recieve qualitative results which matched known empirical samples of freshwater benthic macroinvertebrates  

# load the libraries
library(sizeSpectra)
library(poweRlaw)
library(parallel)
library(foreach)
library(doParallel)
library(tidyverse)

# source the custom functions written for this simulation study
source("custom_functions.R")
# source variable designations, including setting seed
source("master_variable_designation.R")

# make a dataframe with the simulation parameters
sub_lambda_df <- expand_grid(
  n = n,
  lambda = lambda,
  xmin = xmin, 
  xmax = xmax,
  vecDiff = vecDiff,
  pr_scenarios)
# for working out issues
# sub_lambda_df <- sub_lambda_df |>
#   sample_n(10)

# set up parallel processing
#cores <- detectCores()-1 # when running on its own
cores <- 7 # when running with other simulations

cluster <- makeCluster(cores)
registerDoParallel(cluster)

# run simulation 

#n_iter <- 500
sub_lambda_df <- expand_grid(
  sub_lambda_df, 
  rep = 1:n_iter
)

# example data
# used to work out dopar for each row in a df
# df <- expand_grid(
#   n = c(1, NA),
#   mean = c(-2, 2), 
#   rep = 1:3
# )
# y = foreach(j = 1:nrow(df),
#             .combine = rbind,
#             .errorhandling = "remove") %dopar% {
#               x = rnorm(n = df[j,]$n,
#                         mean = df[j,]$mean)
#               x
#               out = data.frame(df[j,], value = x)
#               out}
# y

tictoc::tic()
results <- foreach(j = 1:nrow(sub_lambda_df), 
                   .combine = rbind,
                   .errorhandling = "remove",
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw")) %dopar%{
    # sets a new seed for each parallel instance
    # but makes this reproducible for re-running simulation
    set.seed(j+1130)
    df_out <- plot_sub_lambda(
      n = sub_lambda_df[j,]$n, 
      lambda = sub_lambda_df[j,]$lambda,
      xmin = sub_lambda_df[j,]$xmin, 
      xmax = sub_lambda_df[j,]$xmax, 
      h = sub_lambda_df[j,]$h, 
      b = sub_lambda_df[j,]$b,
      vecDiff = sub_lambda_df[j,]$vecDiff,
      plot = FALSE)
    df_out$pr_scenario <- sub_lambda_df[j,]$pr_scenario
    df_out$rep <- sub_lambda_df[j,]$rep
    df_out
                                 }
#results
tictoc::toc()


# old dopar solution
# I think this was "removing" the whole rep even if only one set of parameter values had an error in it
# i.e., results always had equal counts of reps (1 per parameter set) but not all reps were included

# tictoc::tic()
# results <- foreach(j = 1:n_iter, 
#                    .combine = rbind,
#                    .errorhandling = "remove",
#                    .packages = c("sizeSpectra",
#                                  "tidyverse",
#                                  "poweRlaw")) %dopar%{
#                      df_out <- parallel_rep_sub_lambda(df = sub_lambda_df)
#                      df_out$rep <- j
#                      #results[[j]] <- 
#                      df_out
#                                  }
# results
# tictoc::toc()

stopCluster(cl = cluster)
end <- tictoc::toc()
run <- end$callback_msg
saveRDS(run, paste0("simulation_results/sim_run_time_s_", Sys.Date(), ".rds"))


#result_df <- bind_rows(results)

saveRDS(results, "simulation_results/undersampling_sim_run.rds")

