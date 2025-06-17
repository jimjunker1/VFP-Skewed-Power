# this script looks to determine the effects of changing the origianl N as well as the subsample N

# Objectives
# What N and subsample_N should be chosen?
# what is the variation in the body size which has X% probability of being sampled? i.e., 80%, 90%, 95%
# how does this alter the estimated x_min from poweRlaw?

library(tidyverse)
library(sizeSpectra)
library(poweRlaw)

sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

set.seed(6928)

N_orig = 10000
N_sub = 1000
lambda = -2
xmin = 0.001
xmax = 100

# this doesn't actually need N's. 
# it just summarizes the x_value for different sample probabilities
# could change this to x_pr_summary() or something
# change inputs to h and/or b
# this would tell me what the xmin "should be"
# --> once I firgure out the sampling probability that is ~xmin i..e, 80%, 90%, etc. 
x_pr_subsample <- function(
    N_orig,
    N_sub,
    lambda, 
    xmin, 
    xmax,
    h) {
  x <- rPLB(n = N_orig,
            b = lambda,
            xmin = xmin, 
            xmax = xmax)
  
  x_pr <- data.frame(x = x) |>
    mutate(pr = sample_pr(x, 
                          b = 2,
                          h = h))
  x_90 <- x_pr |>
    filter(pr >=0.89,
            pr <=0.91) |>
    summarize(x_90 = mean(x, na.rm = TRUE))
  x_80 <- x_pr |>
    filter(pr >=0.79,
           pr <=0.81) |>
    summarize(x_80 = mean(x, na.rm = TRUE))
  
  x_95 <- x_pr |>
    filter(pr >=0.94,
           pr <=0.96) |>
    summarize(x_95 = mean(x, na.rm = TRUE))
  
  x_85 <- x_pr |>
    filter(pr >=0.84,
           pr <=0.86) |>
    summarize(x_85 = mean(x, na.rm = TRUE))
  
  
  out <- bind_cols(x_80, 
                   x_85,
                   x_90,
                   x_95)
  out$N_orig = N_orig
  out$N_sub = N_sub
  return(out)
}


x_pr_subsample(10000, 5000, lambda = -2, xmin = 0.001, xmax = 100, h = 0.001)

x_pr_df <- expand_grid(
  N_orig = c(5000, 10000, 50000),
  N_sub = c(100, 500, 1000), 
  lambda = -2, 
  xmin = 0.001, 
  xmax = 100
)

# this probably doesn't need to be here. This just repeats the original x_pr function above, but it doesn't chnage becasue the sample_pr() has set inputs
rep_pr_subsample <- function(x_pr_df){
  out <- NULL
  n_row <- nrow(x_pr_df)
  for (i in 1:n_row){
    sim_in <- x_pr_df[i,]
    out[[i]] <- x_pr_subsample(
      N_orig = sim_in$N_orig,
      N_sub = sim_in$N_sub,
      lambda = sim_in$lambda,
      xmin = sim_in$xmin,
      xmax = sim_in$xmax
    )
  }
  out_df <- bind_rows(out)
  return(out_df)
}

rep_pr_subsample(x_pr_df = x_pr_df)

# xmin for subsamples
xmin_subsample <- function(
    N_orig,
    N_sub,
    lambda, 
    xmin, 
    xmax) {
  x <- rPLB(n = N_orig,
            b = lambda,
            xmin = xmin, 
            xmax = xmax)
  
  x_pr <- data.frame(x = x) |>
    mutate(pr = sample_pr(x, 
                          b = 2,
                          h = 0.001))
  x_sub <- x_pr |>
    sample_n(size = N_sub, weight = pr)
  
  x_power <- conpl$new(x_sub$x)
  est_xmin <- estimate_xmin(x_power)$xmin
  
  xmin_obs <- min(x)
  xmax_obs <- max(x)
  
  out <- data.frame(est_xmin = est_xmin,
                    N_orig = N_orig,
                    N_sub = N_sub,
                    lambda = lambda,
                    xmin_bound = xmin, 
                    xmax_bound = xmax,
                    xmin_obs = xmin_obs,
                    xmax_obs = xmax_obs)
  return(out)
}

xmin_subsample(10000, 500, -2, 0.001, 100)

rep_xmin_subsample <- function(df){
  out <- NULL
  n_row <- nrow(df)
  for (i in 1:n_row){
    sim_in <- df[i,]
    out[[i]] <- xmin_subsample(
      N_orig = sim_in$N_orig,
      N_sub = sim_in$N_sub,
      lambda = sim_in$lambda,
      xmin = sim_in$xmin,
      xmax = sim_in$xmax
    )
  }
  out_df <- bind_rows(out)
  return(out_df)
}


sub_df <- expand_grid(
  N_orig = c(5000, 10000),
  N_sub = c(500, 1000), 
  lambda = -2, 
  xmin = 0.001, 
  xmax = 100
)

est_xmin_df <- rep_xmin_subsample(df = sub_df)

ggplot(est_xmin_df,
       aes(x = N_orig,
           y = est_xmin, 
           color = N_sub)) +
  geom_point() +
  scale_x_log10()

# rep_xmin_subsample with parallel processing
library(parallel)
library(foreach)
library(doParallel)
cores <- detectCores()-1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
n_iter <- 2
results <- list()
tictoc::tic()
results <- foreach(j = 1:n_iter,
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw"),
                   .errorhandling = "pass") %dopar%{
  df_out <- rep_xmin_subsample(df = sub_df)
  df_out$rep <- j
  results[[j]] <- df_out
}
stopCluster(cl = cluster)
tictoc::toc()
results
bind_rows(results)
# repeat the subsample function above for multiple repetitions
# need to parallelize this
repeat_xmin_sim <- function(n, sim_df){
  
  out_list <- list()
  
  for (j in 1:n){
    out_list[[j]] <- rep_xmin_subsample(sim_df)
  }
  out_df <- bind_rows(out_list)
  return(out_df)
}

tictoc::tic()
rep_est_xmin <- repeat_xmin_sim(n = 2, sim_df = sub_df)
tictoc::toc()
# three orig and 3 subs:
# n = 5 ~40s
# n = 50 ~ 402s

# 4 orig and 4 sub
# n = 10 ~134s

ggplot(rep_est_xmin,
       aes(x = N_orig,
           y = est_xmin, 
           color = N_sub,
           group = N_sub)) +
  geom_point() +
  scale_x_log10() +
  stat_smooth(method = "lm")

rep_est_xmin |>
  group_by(N_orig, N_sub) |>
  summarize(m = mean(est_xmin),
            s = sd(est_xmin))

rep_est_xmin |>
  group_by(N_orig, N_sub) |>
  summarize(m = mean(est_xmin),
            s = sd(est_xmin)) |>
  ggplot(aes(x = N_orig, 
             color = N_sub,
             y = m,
             ymin = m-s, 
             ymax = m +s)) +
  geom_pointrange(
    position = position_jitter(width = 10)
  )

sub_df2 <- expand_grid(
  N_orig = c(1000),
  N_sub = c(100, 250, 500, 900), 
  lambda = -2, 
  xmin = 0.001, 
  xmax = 100
)
tictoc::tic()
rep_est_xmin2 <- repeat_xmin_sim(n = 10, sim_df = sub_df2)
tictoc::toc()
rep_est_xmin2
ggplot(rep_est_xmin2,
       aes(x = N_sub,
           y = est_xmin)) +
  geom_point() +
  stat_smooth()




# different subsample idea using rbinom
plot_sub_lambda <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.001, 
    xmax = 100,
    h = 0.00001, 
    b = 2, 
    plot = FALSE){
  # sample from bounded power law
  x <- rPLB(n = n, b = lambda, xmin = xmin, xmax = xmax)
  xmin_obs = min(x)
  xmax_obs = max(x)
  # make a df with the sample probability and not/sampled columns
  x_df <- x |> 
    tibble() |>
    mutate(pr = sample_pr(
      x, 
      h = h, 
      b = b
    ), 
    sampled = rbinom(n(), 1, prob = pr),
    fill = case_when(sampled == 1 ~ "sampled", 
                     .default = "not sampled"))
  
  # under sampled ####
  # filter out the "sampled" data
  # this represents empirical data which has fewer little things than expected
  x_under <- x_df |>
    filter(fill == "sampled")
  # how many body sizes were sampled?
  under_n <- nrow(x_under)
  
  # estimate x_min ####
  # use poweRlaw to estimate where x_min is, i.e., x < x_min is under sampled
  x_power <- conpl$new(x_under$x)
  x_xmin <- estimate_xmin(x_power)$xmin
  
  # add estimated x_min to data frames
  x_df$est_xmin <- x_xmin
  x_under$est_xmin <- x_xmin
  
  # lambdas ####
  # lambda under ####
  # estimate lambda from undersampled data
  x_under_vector <- x_under$x
  lambda_under <- calcLike(negLL.fn = negLL.PLB,
                           x = x_under_vector,
                           xmin = min(x_under_vector), 
                           xmax = max(x_under_vector), 
                           n = length(x_under_vector), 
                           sumlogx = sum(log(x_under_vector)), 
                           p = -1.5,
                           suppress.warnings = TRUE,
                           vecDiff = 2)
  
  # lambda trimmed ####
  # estimate lambda from trimmed data
  x_trimmed_vector <- x_under |>
    filter(x > est_xmin) |>
    pull(x)
  trimmed_n <- length(x_trimmed_vector)
  lambda_trimmed <- calcLike(negLL.fn = negLL.PLB,
                             x = x_trimmed_vector,
                             xmin = min(x_trimmed_vector), 
                             xmax = max(x_trimmed_vector), 
                             n = length(x_trimmed_vector), 
                             sumlogx = sum(log(x_trimmed_vector)), 
                             p = -1.5,
                             suppress.warnings = TRUE,
                             vecDiff = 2)
  # plot ####
    if(plot == TRUE){
    p <- ggplot(x_df,
                aes(x = x, 
                    fill = fill)) +
      geom_histogram(binwidth = 0.01) +
      scale_fill_manual(values = c(NA, "black"))+
      scale_x_log10() +
      geom_vline(aes(xintercept = est_xmin),
                 linetype = "dashed", color = "red") +
      # scale_y_log10() +
      theme_bw() +
      annotate("text",
               x = 0.01, 
               y = 110, 
               label = paste("Under: ",
                             round(lambda_under$MLE, 2))) +
      annotate("text",
               x = 0.01, 
               y = 95, 
               label = paste("Trimmed: ",
                             round(lambda_trimmed$MLE, 2))) +
      annotate("text",
               x = 0.01, 
               y = 125, 
               label = paste("Real \u03bb:",
                             lambda)) +
      annotate("text",
               x = 0.5, 
               y = 125, 
               label = paste("Original N:",
                             n)) +
      annotate("text",
               x = 0.5, 
               y = 110, 
               label = paste("Sub N:",
                             under_n)) +
      annotate("text",
               x = 0.5, 
               y = 70, 
               label = paste("p(sample) variables \nh:",
                             h,
                             "\nb:",
                             b)) +
      annotate("text",
               x = 0.02, 
               y = 60, 
               label = paste("estimated xmin: \n",
                             round(x_xmin,
                                   4))) +
      NULL
  print(p)}
    # end/return ####
    
    
  # return data frame ####
  out_df <- data.frame(
    known_lambda = lambda, 
    lambda_under = lambda_under$MLE,
    lambda_under_lo = lambda_under$conf[1],
    lambda_under_hi = lambda_under$conf[2],
    lambda_trimmed = lambda_trimmed$MLE,
    lambda_trimmed_lo = lambda_trimmed$conf[1],
    lambda_trimmed_hi = lambda_trimmed$conf[2],
    original_n = n,
    under_n = under_n,
    trimmed_n = trimmed_n,
    est_xmin = x_xmin, 
    h = h, 
    b = b,
    xmin_obs = xmin_obs,
    xmax_obs = xmax_obs)
  return(out_df)
  }

plot_sub_lambda(n = 5000, h = 0.00004, b = 2)

###
sub_lambda_df <- expand_grid(
  n = c(5000, 10000),
  lambda = -2,
  xmin = 0.001,
  xmax = 100,
  h = c(0.00004),
  b = 2)
###

parallel_rep_sub_lambda <- function(df){
  out <- NULL
  # number of sim parameter sets
  n_row <- nrow(sub_lambda_df)
  # empty list for reps
    for(i in 1:n_row){
      df_in <- df[i,]
      out[[i]] <- plot_sub_lambda(
        n = df_in$n, 
        lambda = df_in$lambda,
        xmin = df_in$xmin, 
        xmax = df_in$xmax, 
        h = df_in$h, 
        b = df_in$b, 
        plot = FALSE)
    }
  out_df <- bind_rows(out)
  return(out_df)
}
parallel_rep_sub_lambda(df = sub_lambda_df)

cores <- detectCores()-1
cluster <- makeCluster(cores)
registerDoParallel(cluster)
n_iter <- 10
results <- list()
tictoc::tic()
results <- foreach(j = 1:n_iter,
                   .packages = c("sizeSpectra",
                                 "tidyverse",
                                 "poweRlaw"),
                   .errorhandling = "pass") %dopar%{
  df_out <- parallel_rep_sub_lambda(df = sub_lambda_df)
  df_out$rep <- j
  results[[j]] <- df_out
                   }
stopCluster(cl = cluster)
tictoc::toc()
results
bind_rows(results)
# 11 seconds n rep = 2
# 10 reps = 18s

rep_sub_lambda <- function(sub_lambda_df, rep = 2){
  # number of sim parameter sets
  n_row <- nrow(sub_lambda_df)
  # empty list for reps
  final_out <- list()
  param_out <- list()
  for(j in 1:rep){
    for(i in 1:n_row){
      df_in <- sub_lambda_df[i,]
      param_out[[i]] <- plot_sub_lambda(
        n = df_in$n, 
        lambda = df_in$lambda,
        xmin = df_in$xmin, 
        xmax = df_in$xmax, 
        h = df_in$h, 
        b = df_in$b, 
        plot = FALSE)
    }
    final_out[[j]] <- bind_rows(param_out)
    final_out[[j]]$rep <- j
  }
  out_df <- bind_rows(final_out)
  return(out_df)
}

sub_lambda_df <- expand_grid(
  n = c(5000, 10000),
  lambda = -2,
  xmin = 0.001,
  xmax = 100,
  h = c(0.00004),
  b = 2)
tictoc::tic()
rep_sub_lambda(sub_lambda_df = sub_lambda_df, 
               rep = 10)
tictoc::toc()
# sequential with 2 reps = 14.3s
# 10 reps = 72s



#parallel processing for sub sampling and estimating lambdas
