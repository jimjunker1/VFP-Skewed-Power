# custom functions for simulating effects of under sampling small body sizes


# custom functions ####
# sample_pr ####

# this function caluclates the sampling probability of body size x
# the variables to control the response are h and b  
# h approximately controls where the probability starts to increase  
# b controls how rapidly the probability of sampling increases  
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

# plot_sub_lambda ####
# n = original sample size  
# lambda = known value of lambda to sample from
# xmin and xmax = set range of body sizes (x) for the bounded power law
# h and b are variables which are passed to the sample_pr() function internally 
plot_sub_lambda <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.01, 
    xmax = 100,
    h = 0.001, 
    b = 2.5, 
    vecDiff = 2,
    binwidth = 0.1,
    plot = FALSE,
    annotate = FALSE){
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
                           vecDiff = vecDiff)
  
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
    
    plot_x_df <- tibble(x = x, 
                        fill = "Original")
    
    plot_df <- bind_rows(
      x_under |>
        select(x, fill),
      plot_x_df)
    plot_df$est_xmin <- x_xmin
    
    p <- ggplot(plot_df,
                aes(x = x, 
                    fill = fill)) +
      geom_histogram(binwidth = binwidth, 
                     position = "dodge2") +
      scale_fill_manual(values = c("black", "dodgerblue"))+
      scale_x_log10() +
      geom_vline(aes(xintercept = est_xmin),
                 linetype = "dashed", color = "red") +
      # scale_y_log10() +
      theme_bw() 
    if(annotate == TRUE) {
      p +
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
                                     4)))
    }
    print(p)
    }
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


# parallel_rep_sub_lambda ####
# this is a helper function to run the simulations in parallel
parallel_rep_sub_lambda <- function(df){
  out <- NULL
  # number of sim parameter sets
  n_row <- nrow(df)
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
      vecDiff = df_in$vecDiff,
      plot = FALSE)
    out[[i]]$pr_scenario <- df_in$pr_scenario
  }
  out_df <- bind_rows(out)
  return(out_df)
}


# samples x, undersamples, estimate xmin and trim, estimate lambda, back calculate count/distribution of body sizes
plot_unsub_lambda <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.01, 
    xmax = 100,
    h = 0.001, 
    b = 2.5, 
    vecDiff = 2,
    plot = TRUE,
    plot_back_est_n = FALSE,
    bin_width = 0.1){
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
  
  # back calc body size counts
  est_tot_n <- estimate_pareto_N(
    n = trimmed_n,
    lambda = lambda_trimmed$MLE,
    xmin =  x_df$est_xmin[1],
    xmin2 = xmin,
    xmax = xmax)
  
  x_vector <- seq(xmin, xmax, by = 0.01)
  
  x_vector_est <- rPLB(n = est_tot_n,
                       b = lambda_trimmed$MLE,
                       xmin = xmin,
                       xmax = xmax) #* est_tot_n
  
  x_est <- data.frame(x = x_vector_est,
                      fill = "estimated")
  # hi and lo estimates
  x_vector_est_lo <- rPLB(n = est_tot_n,
                          b = lambda_trimmed$conf[1],
                          xmin = xmin,
                          xmax = xmax) #* est_tot_n
  
  x_est_lo <- data.frame(x = x_vector_est_lo,
                         fill = "estimated low")
  
  # hi and lo estimates
  x_vector_est_hi <- rPLB(n = est_tot_n,
                          b = lambda_trimmed$conf[2],
                          xmin = xmin,
                          xmax = xmax) #* est_tot_n
  
  x_est_high <- data.frame(x = x_vector_est_hi,
                           fill = "estimated high")
  
  x_orig_df <- data.frame(x = x,
                          fill = "original")
  
  x_df2 <- bind_rows(x_orig_df, x_est, x_under) #, x_est_high, x_est_lo)
  x_df2$est_xmin <- x_xmin
  
  # plot ####
  if(plot == TRUE & plot_back_est_n == TRUE){
    p <- ggplot(x_df2,
                aes(x = x, 
                    fill = fill)) +
      geom_histogram(binwidth = bin_width,
                     position = "dodge2") +
      scale_fill_manual(values = c("dodgerblue",
                                   #"coral", "maroon",
                                   "black", "pink"))+
      scale_x_log10() +
      geom_vline(aes(xintercept = est_xmin),
                 linetype = "dashed", color = "black") +
      # scale_y_log10() +
      theme_bw() 
    
    print(p)
  }
  
  if(plot == TRUE & plot_back_est_n == FALSE){
    x_df3 <- bind_rows(x_orig_df, x_under)
    p <- ggplot(x_df3,
                aes(x = x, 
                    fill = fill)) +
      geom_histogram(binwidth = bin_width,
                     position = "dodge2") +
      scale_fill_manual(values = c(#"dodgerblue",
                                   #"coral", "maroon",
                                   "black", "pink"))+
      scale_x_log10() +
      # geom_vline(aes(xintercept = est_xmin),
      #            linetype = "dashed", color = "black") +
      # scale_y_log10() +
      theme_bw() 
    
    print(p)
  }
  # end/return ####
  
}

# estimates the total count from a lambda and x-range to a new range of x (xmin2 == "smaller" x)
estimate_pareto_N = function(n, lambda, xmin, xmin2, xmax){
  lambdaPlus = lambda + 1
  n * (xmax^(lambdaPlus) - xmin2^(lambdaPlus)) /
    (xmax^(lambdaPlus) - xmin^(lambdaPlus))
}

# compare total count from under sampled and trimmed data
compare_est_pareto_n <- function(
    n = 5000, 
    lambda = -2,
    xmin = 0.01, 
    xmax = 100,
    h = 0.001, 
    b = 2.5, 
    vecDiff = 2){
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
  
  # lambda under sampled ####
  # estimate lambda from under sampled data
  x_under_vector <- x_under |>
    pull(x)
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
  
  # pareto N ####
  # back calc body size counts
  # trimmed
  est_tot_n_trimmed <- estimate_pareto_N(
    n = trimmed_n,
    lambda = lambda_trimmed$MLE,
    xmin =  x_df$est_xmin[1],
    xmin2 = xmin,
    xmax = xmax)
  
  # under sampled
  est_tot_n_under <- estimate_pareto_N(
    n = under_n,
    lambda = lambda_under$MLE,
# x_min ???####
    xmin =  x_df$est_xmin[1],#should this just be x_min?? 
    xmin2 = xmin,
    xmax = xmax)
  
  out <- data.frame(original_n = n,
                    under_n = under_n, 
                    trimmed_n = trimmed_n,
                    est_under_n = est_tot_n_under,
                    est_trimmed_n = est_tot_n_trimmed,
                    h = h, 
                    b = b, 
                    lambda = lambda,
                    xmin = xmin,
                    xmax = xmax)
  return(out)
}
# function to repeat samples of estimated pareto N
parallel_rep_compare_pareto_n <- function(df){
  out <- NULL
  # number of sim parameter sets
  n_row <- nrow(df)
  # empty list for reps
  for(i in 1:n_row){
    df_in <- df[i,]
    out[[i]] <- compare_est_pareto_n(
      n = df_in$n, 
      lambda = df_in$lambda,
      xmin = df_in$xmin, 
      xmax = df_in$xmax, 
      h = df_in$h, 
      b = df_in$b,
      vecDiff = df_in$vecDiff)
    out[[i]]$pr_scenario <- df_in$pr_scenario
  }
  out_df <- bind_rows(out)
  return(out_df)
}
