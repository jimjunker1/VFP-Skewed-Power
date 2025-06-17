# function for sinusoidal response of sampling probabilility

# function ####
sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

# SI figures for sampling probability with simulation values
library(tidyverse)
library(sizeSpectra)

source("master_variable_designation.R")
set.seed(2112)
dat <- expand_grid(
  h = h,
  b = b,
  x = c(exp(seq(log(0.001), log(10), length.out = 1000)))) |>
mutate(pr = sample_pr(x = x, h = h, b = b))

dat |>
  ggplot(aes(x = x, 
             y = pr, 
             color = as.factor(h))) +
  geom_point() +
  geom_vline(aes(xintercept = 0.005), linetype = "dashed") +
  geom_hline(aes(yintercept = 0.9), linetype = "dashed") +
  scale_x_log10() +
  facet_wrap(b~., labeller = label_both) +
  theme_bw() +
  guides(color = guide_legend(title= "h")) +
  labs(title = "Sampling probability as a function of body mass",
       x =expression(Log[10]~dry~mass),
       y = "Sampling probability",
       caption = "The horizontal dashed line represents a 90% sampling probability for reference. The vertical line shows an individual of 0.1 mg dry mass. \nDashed lines are chosen arbritarily for illustrative purposes.")

ggsave("plots/sampling_probability.png", scale = 2)
# Example of undersampling ####
# make SI figures of "real" data and undersampling results

# simulate x from bounded power law
original_dat <- tibble(
  x = rPLB(5000, b = -2, xmin = xmin, xmax = xmax))

# expand grid to include each variable combination with each body size
sample_dat <- expand_grid(
  original_dat, 
  h = h, 
  b = b
)

# add sample probabilities as a function of x (body mass)
sampled_dat <- sample_dat |> 
  mutate(pr = sample_pr(
    x, 
    h = h, 
    b = b
  ), 
  sampled = rbinom(n(), 1, prob = pr),
  fill = case_when(sampled == 1 ~ "sampled", 
                   .default = "not sampled"))

# add "fill = original" to simulated data
original_x <- original_dat |>
  mutate(fill = "original")

# combine "original and sampled data
both_dats <- left_join(sampled_dat, original_x, by = "x") |>
  pivot_longer(fill.x:fill.y)

both_dats |>
  filter(value!= "not sampled") |>
  ggplot(aes(x = x, 
             fill = value)) +
  geom_histogram(binwidth = 0.1, 
                 position = "dodge2") +
  scale_x_log10() +
  facet_wrap(b~h, labeller = label_both) +
  scale_fill_manual(values = c("dodgerblue","black")) + 
  theme_bw() +
  guides(color = guide_legend(title= "Observation")) +
  labs(title = "Undersampling a power law distribution",
       x =expression(Log[10]~dry~mass),
       caption = "Examples of how samples from a power law distribution (blue bars) are affected by different combinations of variables (facet titles) \nin the sampling probability equation, and how under sampled data (black bars) appears. \nUnder sampling is minor in top left and increases to extreme in the bottom right. Original N = 5000, \u03bb = -2")
ggsave("plots/original_under_sample_example.png", scale = 2)



x_pr_subsample <- function(
    df
    #N_orig,
    #N_sub,
    #lambda, 
    #xmin, 
    #xmax,
    #h,
    #b
    #x_pr
    ) {
  # x <- rPLB(n = N_orig,
  #           b = lambda,
  #           xmin = xmin, 
  #           xmax = xmax)
  
  # x_pr <- data.frame(x = x) |>
  #   mutate(pr = sample_pr(x, 
  #                         b = 2,
  #                         h = h))
  x_90 <- df |>
    filter(pr >=0.89,
           pr <=0.91) |>
    group_by(h, b) |>
    summarize(x_90 = mean(x, na.rm = TRUE)) |>
    ungroup() |>
    select(-h, -b)
  
  x_80 <- df |>
    filter(pr >=0.79,
           pr <=0.81) |>
    group_by(h, b) |>
    summarize(x_80 = mean(x, na.rm = TRUE))
  
  x_95 <- df |>
    filter(pr >=0.94,
           pr <=0.96) |>
    group_by(h, b) |>
    summarize(x_95 = mean(x, na.rm = TRUE))|>
    ungroup() |>
    select(-h, -b)
  
  x_85 <- df |>
    filter(pr >=0.84,
           pr <=0.86) |>
    group_by(h, b) |>
    summarize(x_85 = mean(x, na.rm = TRUE))|>
    ungroup() |>
    select(-h, -b)
  
  out <- bind_cols(x_80, 
                   x_85,
                   x_90,
                   x_95)
  #out$N_orig = N_orig
  #out$N_sub = N_sub
  return(out)
}

dat %>%
  #group_by(h, b) %>%
  reframe(x_pr_subsample(.)) |>
  arrange(b, h) |>
  mutate(across(starts_with("x_"), ~ round(.x, 3))) |>
  write_csv("simulation_summaries/body_size_sampling_probabillities.csv")

dat %>%
  #group_by(h, b) %>%
  reframe(x_pr_subsample(.)) |>
  pivot_longer(x_80:x_95) |>
  separate(name, into = c("x", "probability")) |>
  mutate(probability = as.numeric(probability)) |>
  ggplot(aes(x = value, 
             y = probability,
             color = interaction(h, b))) +
  geom_point()
  