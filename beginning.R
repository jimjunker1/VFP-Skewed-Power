library(sizeSpectra)
library(tidyverse)


# x <- rPLB(n = 10000, b = -2, xmin = 0.001, xmax = 100)
# 
mle_and_plot <- function(x, lbn_p = FALSE){

  x_binned <- binData(x = x, binWidth = "2k")
  x_mle <- calcLike(negLL.fn = negLL.PLB,
                    x = x,
                    xmin = min(x), # 0.0026
                    xmax = max(x),
                    n = length(x),
                    sumlogx = sum(log(x)),
                    p = -1.5,
                    suppress.warnings = TRUE)
  x_mle

  x_binvals <- x_binned$binVals

  if(lbn_p == TRUE){
  LBN_bin_plot(
    binValsTibble = x_binvals[complete.cases(x_binvals),],
    b.MLE = x_mle$MLE,
    b.confMin = x_mle$conf[1],
    b.confMax = x_mle$conf[2],
    log.xy = "xy")}
  out <- data.frame(b = x_mle$MLE,
                    ci_low = x_mle$conf[1],
                    ci_high = x_mle$conf[2])
  return(out)
}
# 
# mle_and_plot(x = x, lbn_p = TRUE)
# 
# x_under <- data.frame(x = sort(x), 
#                       theta = NA)
# 
# x_under <- x_under %>%
#   mutate(theta = case_when(
#     x <= 0.03 ~ 0.01, 
#     .default = 1
#   ))
# 
# tail(x_under)
# 
# 
# 
# x_obs <- x_under %>%
#   sample_n(size = 1000, weight = theta)
# 
# x_obs
# 
# mle_and_plot(x = x_obs$x)
# 
# mle_and_plot(x = x_obs$x[x_obs$x>0.035])


# Simulation 1 iteration --------------------------------------------------

set.seed(6928)
run_sim <- function(sim_df,
                    xmin = 0.001,
                    xmax = 100,
                    vecDiff = 2){
  n_row <- nrow(sim_df)
  
  out_df <- data.frame(
    lambda = numeric(n_row),
    theta = numeric(n_row),
    theta_c = numeric(n_row),
    psi = numeric(n_row),
    lambda_est = numeric(n_row),
    lambda_low = numeric(n_row),
    lambda_high = numeric(n_row),
    cut_lambda_est = numeric(n_row),
    cut_lambda_low = numeric(n_row),
    cut_lambda_high = numeric(n_row))
  
  for(i in 1:nrow(sim_df)){
    sim_in <- sim_df[i, ]
    lambda = sim_in$lambda
    theta = sim_in$theta
    theta_c = sim_in$theta_c
    psi = sim_in$psi
    
    x <- rPLB(n = 10000,
              b = lambda,
              xmin = xmin,
              xmax = xmax)
    
    x_under <- data.frame(x = sort(x))
    
    x_under <- x_under %>%
      mutate(weight = case_when(
        x <= theta ~ psi, 
        .default = 1
      ))
    
    x_obs <- x_under %>%
      sample_n(size = 1000, weight = weight)
    
    
    #x_binned <- binData(x = x_obs$x, binWidth = "2k")
    
    x_mle <- calcLike(negLL.fn = negLL.PLB,
                      x = x_obs$x,
                      xmin = min(x_obs$x), 
                      xmax = max(x_obs$x), 
                      n = length(x_obs$x), 
                      sumlogx = sum(log(x_obs$x)), 
                      p = -1.5,
                      suppress.warnings = TRUE,
                      vecDiff = vecDiff)
    
    x_cut <- x_obs$x[x_obs$x>theta_c]
    
    x_cut_mle <- calcLike(negLL.fn = negLL.PLB,
                          x = x_cut,
                          xmin = min(x_cut), 
                          xmax = max(x_cut), 
                          n = length(x_cut), 
                          sumlogx = sum(log(x_cut)), 
                          p = -1.5,
                          suppress.warnings = TRUE,
                          vecDiff = vecDiff)
    
    out_df[i,] <- data.frame(
      sim_in,
      lambda_est = x_mle$MLE,
      lambda_low = x_mle$conf[1],
      lambda_high = x_mle$conf[2],
      cut_lambda_est = x_cut_mle$MLE,
      cut_lambda_low = x_cut_mle$conf[1],
      cut_lambda_high = x_cut_mle$conf[2])
  }
  return(out_df)
}

l_out <- 3
sim_df <- expand_grid(
  lambda = seq(-1, -2, length.out = 5),
  theta = 0.05, #seq(0.01, 0.1, length.out = 10),
  theta_c = c(0.01, 0.0375, 0.025, 0.05, 0.0675, 0.075, 0.09),#seq(0.01, 0.09, length.out = 10),
  psi = c(0.1, 0.5) #seq(0.1, 0.5, length.out = l_out)
)

# nrow(sim_df)
# tictoc::tic()
# result <- run_sim(sim_df)
# tictoc::toc()

# 50 = 2
# 81 rows = 1.86 sec
# 625 rows = 11.97
# 500 rows = 13.5

# names(result)
# head(result)
# result %>%
#   ggplot(aes(x = theta,
#              y = lambda_est)) +
#   geom_hline(aes(yintercept = lambda)) +
#   geom_point() +
#   facet_wrap(~lambda)
# 
# 
# result %>%
#   #filter(theta > 0.075) %>%
#   select(-contains("low"), 
#          -contains("high")) %>%
#   pivot_longer(lambda_est:cut_lambda_est,
#                names_to = "method", 
#                values_to = "value") %>%
#   ggplot(aes(x = theta_c,
#              y = value,
#              #color = psi,
#              shape = method)) +
#   geom_hline(aes(yintercept = lambda)) +
#   geom_vline(aes(xintercept = theta)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~lambda)

repeat_sim <- function(n, sim_df, vecDiff = 1){
  #n_obs = n * nrow(sim_df)
  
  out_list <- list()
  
  for (j in 1:n){
    out_list[[j]] <- run_sim(sim_df, vecDiff = vecDiff)
  }
  out_df <- bind_rows(out_list)
  return(out_df)
}

nrow(sim_df)
n = 100
n*nrow(sim_df)

tictoc::tic()
rep_res <- repeat_sim(n = n,
                      sim_df = sim_df,
                      vecDiff = 5)
tictoc::toc()
# 25 * 2 = 50; vecDiff = 10 --> 8.5 seconds 
# 25 * 10 = 250; vecDiff = 10 --> 32
# 25 * 100 = 2500; vecdiff = 1 --> 57
# 70 * 100 = 7000; vecDiff = 5 --> 411

rep_res %>%
  #filter(theta > 0.075) %>%
  select(-contains("low"), 
         -contains("high")) %>%
  pivot_longer(lambda_est:cut_lambda_est,
               names_to = "method", 
               values_to = "value") %>%
  ggplot(aes(x = theta_c,
             y = value,
             color = method)) +
  geom_point(alpha = 0.6) +
  geom_smooth(alpha = 0.6,
              se = FALSE) +
  geom_hline(aes(yintercept = lambda),
             linetype = 2,
             linewidth = 1,
             color = "red") +
  geom_vline(aes(xintercept = theta), 
             linetype = 2,
             linewidth = 1,
             color = "black") +
  scale_color_viridis_d(option = "inferno",
                        begin = 0.25, 
                        end = 0.75) +
  facet_wrap(
    lambda~psi,
    ncol = 2,
    scales = "free_y",
    labeller = 
      labeller(
        lambda = ~ paste(expression(lambda), ":", .),
        psi = ~ paste(expression(psi), ":", .),
        .multi_line = FALSE
      )) +
  theme_bw() +
  labs(y = "estimated lambda",
       x = "Minimum cut off size",
       caption = "Relationship between lambda estimate and known value (red horizontal dashed line and facet title). Psi is the probability of undersampling body sizes below a threshold. \n Here, the threshold is kept constant at 0.05 and is represented by a vertical black dashed line. The orange line is the lambda estimate when sizes are undersampled but all data is used. \n The purple line is the estimate when data is cutoff at different values (x-axis). The estimates improve when undersampled data is removed, even if the data is cutoff 'above' the undersampled value")

names(rep_res)
rep_res %>%
  mutate(l_diff = lambda - lambda_est,
         cut_diff = lambda - cut_lambda_est) %>%
  group_by(lambda, theta_c) %>%
  summarize(mean(l_diff), sd(l_diff),
            mean(cut_diff), sd(cut_diff)) %>%
  View()

