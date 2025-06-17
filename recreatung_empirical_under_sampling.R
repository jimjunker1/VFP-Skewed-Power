# empirical data to estimate undersampling

library(tidyverse)
library(poweRlaw)
library(sizeSpectra)

# empirical data from project on the Telluride Valley Floor
tvf <- read_csv("empirical_data_examples/tvf_dw.csv") |>
  filter(!is.na(dw))

min(tvf$dw) # 0.000596 mg
# plot raw distribution of body masses  
ggplot(tvf, aes(x = dw)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Original sample of x") +
  facet_wrap(~site_number)

min(tvf$Value) # 0.5 mm 
# plot raw distribution of body lengths  
ggplot(tvf, aes(x = Value)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Original sample of x") +
  facet_wrap(~site_number)

# estimate xmin using the poweRlaw package from Clauset
# 1) first, make a list of data sets  
tvf_list <- tvf %>% 
  group_by(site_number) %>% 
  group_split()

# 2) create empty list to fill with xmins
xmin_list = list()

# 3) get list of xmins for each sample
for(i in 1:length(tvf_list)){
  powerlaw = conpl$new(tvf_list[[i]]$dw) # get power law estimate from poweRlaw package
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin, # extract the xmin from the poweRlaw package
                          site_number = unique(tvf_list[[i]]$site_number))
}

# get the xmin for each site
xmins_clauset = bind_rows(xmin_list)
xmins_clauset

# join tvf data with estimated xmin
tvf <- left_join(tvf, xmins_clauset)

# add a logical column for display and analysis
tvf <- tvf |>
  mutate(undersampled = dw<xmin_clauset)

# plot data showing undersampled values
ggplot(tvf, 
       aes(x = dw,
           fill = undersampled)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Original sample of x") +
  facet_wrap(~site_number)

# estimate lambda with trimmed data for site_number 01
tvf_01_trimmed <- tvf |>
  filter(site_number == "01",
         undersampled == FALSE) |> 
  pull(dw)

# mle estimate
lambda_01_trimmed <- calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = tvf_01_trimmed, # the vector of data
  xmin = min(tvf_01_trimmed), # the minimum body size
  xmax = max(tvf_01_trimmed), # the maximum body size
  n = length(tvf_01_trimmed), # the number of observations
  sumlogx = sum(log(tvf_01_trimmed)), # sum of log-transformed data
  p = -1.5
)
lambda_01_trimmed

# simulate "correct" data set
# first, extract original site 01 data
tvf_01 <- tvf |>
  filter(site_number == "01") |> 
  pull(dw)

est_m <- rPLB(n = 500000, # this consistenly gets similar max values
              b = lambda_01_trimmed$MLE,
              xmin = min(tvf_01),
              xmax = max(tvf_01))
min(est_m)
max(est_m)
max(tvf_01)

sim_correct_orig <- data.frame(est_m = est_m,
                          xmin_c = xmins_clauset[1,1]) |>
  mutate(undersampled = est_m<xmin_clauset)

# simulated data
ggplot(sim_correct_orig, 
       aes(x = est_m,
           fill = undersampled)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Simulated 'correct' body sizes")
# compare to empirical data
ggplot(tvf |> filter(site_number == "01"), 
       aes(x = dw,
           fill = undersampled)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Original sample of x") 

# recreate undersampled distribution

# this is a function which simulates the probability of an individual being sampled based on body length
morin <- function(RL, a = -2.84, b = 5.80, c = -3.18, M = 0.25){
  # a, b, and c are fitted coefficients
  # RL is relative body length (body size / mesh size)
  # M = mesh size. 0.25 mm = 250 microns
  p <- (a + b * log10(RL) + c * log10(RL)*log10(M))
  out <- exp(p) / (1 + exp(1.8))
  return(out)
}

# estimate length from mass
# M = aL^b --> L = (m/a)^(1/b)
# general a and b for all insects are 0.0064 and 2.788 (Table 2 Benke et al. 1999)
sim_correct_orig <- sim_correct_orig |>
  mutate(L = (est_m / 0.0064)^(1/2.788)) 

# sim_correct|>
#   ggplot(aes(x = est_m,
#              y = L)) +
#   geom_point() +
#   scale_x_log10()

sim_correct <- sim_correct_orig |>
  mutate(RL = L / 0.25,
         p = morin(RL = RL, M = 0.05)) |>
  mutate(p = case_when(p >1 ~ 1, 
                       .default = p))

# ggplot(sim_correct,
#        aes(x = est_m, 
#            y = p)) +
#   geom_point() +
#   xlim(0, 0.01)

sim_under <- sim_correct |>
  sample_n(800, replace = FALSE, weight = p)
  
ggplot(sim_under,
       aes(x = est_m,
           fill = undersampled)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "artificial under sampling",
       subtitle = "Color is 'original' undersampling") 
       

sim_power <- conpl$new(sim_under$est_m)
sim_power_xmin <- estimate_xmin(sim_power)$xmin
sim_power_xmin
sim_under$new_xmin <- sim_power_xmin

sim_under |>
  mutate(new_fill = est_m < new_xmin) |>
  ggplot(
       aes(x = est_m,
           fill = new_fill)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "artificial under sampling",
       subtitle = "Color is 'new' undersampling") 
xmins_clauset


calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = sim_under$est_m, 
  xmin = min(sim_under$est_m), 
  xmax = max(sim_under$est_m), 
  n = length(sim_under$est_m), 
  sumlogx = sum(log(sim_under$est_m)), 
  p = -1.5
)


# sinusoidal ####
# function for a sinusoidal sampling probability
# x = body size
# h = where does slope start?
# b = how steep is decline? large b = large decline
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}


tvf_pr <- tibble(x = tvf_01)

ggplot(tvf_pr,
  aes(x = x)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()

ggplot(tvf |> filter(site_number == "01"), 
       aes(x = dw,
           fill = undersampled)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10() +
  labs(title = "Original sample of x") 

tvf_pr |>
  mutate(pr = sample_pr(
    x, 
    h = 0.0004, 
    b = 2
  ),
  fill = pr < 0.95) |>
  ggplot(
         aes(x = x,
             fill = fill)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()

tibble(x = est_m) |>
  mutate(pr = sample_pr(
    x, 
    h = 0.0004, 
    b = 2
  ),
  fill = pr < 0.95) |>
  ggplot(
    aes(x = x,
        fill = fill)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()


sinu_sample <- tibble(x = est_m) |>
  mutate(pr = sample_pr(
    x, 
    h = 0.0004, 
    b = 2
  ),
  fill = pr < 0.95) |>
  sample_n(800, weight = pr) 
sinu_sample|>
  ggplot(
    aes(x = x,
        fill = fill)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()

sinu_power <- conpl$new(sinu_sample$x)
sinu_power_xmin <- estimate_xmin(sinu_power)$xmin
sinu_power_xmin
sinu_sample$new_xmin <- sinu_power_xmin

ggplot(sinu_sample,
  aes(x = x,
      fill = x < new_xmin)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()

xmins_clauset
lambda_01_trimmed$MLE

sinu_trim <- sinu_sample |>
  filter(x > new_xmin) |>
  pull(x)

calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = sinu_trim, 
  xmin = min(sinu_trim), 
  xmax = max(sinu_trim), 
  n = length(sinu_trim), 
  sumlogx = sum(log(sinu_trim)), 
  p = -1.5
)


# how does it work for other sample?
tvf_03 <- tvf |>
  filter(site_number == "03") |> 
  pull(dw)

# estimate lambda with trimmed data for site_number 01
tvf_03_trimmed <- tvf |>
  filter(site_number == "03",
         undersampled == FALSE) |> 
  pull(dw)

# mle estimate
lambda_03_trimmed <- calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = tvf_03_trimmed, # the vector of data
  xmin = min(tvf_03_trimmed), # the minimum body size
  xmax = max(tvf_03_trimmed), # the maximum body size
  n = length(tvf_03_trimmed), # the number of observations
  sumlogx = sum(log(tvf_03_trimmed)), # sum of log-transformed data
  p = -1.5
)
lambda_03_trimmed


est_03 <- rPLB(n = 500000, # this consistenly gets similar max values
              b = lambda_03_trimmed$MLE,
              xmin = min(tvf_03),
              xmax = max(tvf_03))
max(est_03)
max(tvf_03)

ggplot(tibble(x = est_03),
  aes(x = x)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()


sinu_sample_03 <- tibble(x = est_03) |>
  mutate(pr = sample_pr(
    x, 
    h = 0.0004, 
    b = 2
  ),
  fill = pr < 0.85) |>
  sample_n(5000, weight = pr) 

tibble(x = est_03) |>
  mutate(pr = sample_pr(
    x, 
    h = 0.0004, 
    b = 2
  )) |>
  filter(pr > 0.860,
         pr < 0.861) |>
  arrange(pr)

sinu_sample_03|>
  ggplot(
    aes(x = x,
        fill = fill)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()



##
sinu_power <- conpl$new(sinu_sample_03$x)
sinu_power_xmin <- estimate_xmin(sinu_power)$xmin
sinu_power_xmin
sinu_sample_03$new_xmin <- sinu_power_xmin

ggplot(sinu_sample_03,
       aes(x = x,
           fill = x < new_xmin)) +
  geom_histogram(binwidth = 0.01) +
  scale_x_log10()

xmins_clauset
lambda_03_trimmed$MLE

sinu_trim <- sinu_sample_03 |>
  filter(x > new_xmin) |>
  pull(x)

calcLike(
  negLL.fn = negLL.PLB, # continuous estimates of all individuals
  x = sinu_trim, 
  xmin = min(sinu_trim), 
  xmax = max(sinu_trim), 
  n = length(sinu_trim), 
  sumlogx = sum(log(sinu_trim)), 
  p = -1.5
)
