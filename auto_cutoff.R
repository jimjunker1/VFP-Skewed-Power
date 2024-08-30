# automatic cutoff detection

library(sizeSpectra)
library(tidyverse)

x <- rPLB(n = 10000, b = -2, xmin = 0.001, xmax = 100)

x_binned <- binData(x = x,
                    binWidth = "2k")
x_binned$binVals %>% 
  select(binMin, binMax)
# 6th bin is (0.0312, 0.0625) 
# trying to undersample body sizes smaller than 0.0312

(index_peak <- which.max(x_binned$binVals$binSumNorm))


x_under <- data.frame(x = sort(x), 
                      theta = NA)

# psi_x <- function(x, theta, psi_min = 0.0001, psi_max = 0.1){
#   xmin = min(x)
#   xmax = theta
#   x_trim = x[x<=xmax]
#   psi = (x_trim - xmin) / (xmax - xmin) *
#     (psi_max - psi_min) + psi_min
#   return(psi)
# }
# 
# 
# x_under <- x_under %>%
#   mutate(psi = psi_x(x = x, theta = 0.0312))


x_under <- x_under %>%
  mutate(theta = case_when(
    x <= 0.0312 & x > 0.0156 ~ 0.01,
    x <= 0.0156 & x > 0.00781 ~ 0.01,
    x <= 0.00781 & x > 0.00391 ~ 0.001,
    x < 0.00391 ~ 0.0001,
    .default = 1
  ))

head(x_under)
tail(x_under)
dim(x_under)

x_under %>% group_by(theta) %>% count()


x_obs <- x_under %>%
  sample_n(size = 1000, weight = theta)

x_obs %>% group_by(theta) %>% count()

dim(x_obs)
head(x_obs)
tail(x_obs)

x_obs_bin <- binData(x = x_obs$x,
                     binWidth = "2k")

(index_peak <- which.max(x_obs_bin$binVals$binSumNorm))

#mle_and_plot(x, lbn_p = TRUE)
mle_and_plot(x_obs$x, lbn_p = TRUE)

x_obs_bin$binVals$binMax[index_peak]

mle_and_plot(x = x_obs$x[x_obs$x>x_obs_bin$binVals$binMax[index_peak]], lbn_p = TRUE)

sim_auto_c <- function(sim_df,
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