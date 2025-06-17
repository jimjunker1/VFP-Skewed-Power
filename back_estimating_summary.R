# back_estimating_summary

library(tidyverse)
library(tidybayes)
library(poweRlaw)
library(sizeSpectra)

source("custom_functions.R")

# show example of back estimating counts
set.seed(999)
plot_unsub_lambda(n = 5000, 
                  xmin = 0.001, 
                  h = 1e-05, 
                  b = 2)
ggsave("plots/back_est_n_example.png", width = 9.35, height = 4.1)

# read in simulaiton
# note that I forgot to remove the "test" from the name when I saved the full run. This is the most recent and complete data. 
back_dat <- readRDS("simulation_results/est_pareto_N_sim_run_test.rds") |>
  as_tibble() 

# how many replicates were successful?
# should be 36 parameter sets * 500 reps = 18000
back_dat |>
  group_by(rep) |>
  count() |>
  ungroup() |>
  group_by(n) |>
  count()
# 20 of the replicates are missing 1 parameter set
# 480 reps have all parameter sets
# 17,980 / 18,000 = 99.89% of simulations successful

# unique sampling probabillity sets:
back_dat |>
  distinct(pr_scenario, h, b)

# rename the pr_scenarios so 01 = least under sampling and 04 = most undersampling

back_dat <- back_dat |>
  mutate(pr_scenario = 
           case_when(pr_scenario == "03" ~ 1,
                     pr_scenario == "01" ~ 2,
                     pr_scenario == "04" ~ 3,
                     pr_scenario == "02" ~ 4))

# summarize data ####

# just the n = 5000, lambda = -2 
back_dat |>
  filter(original_n == 5000,
         lambda == -2) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05,
                       end = 0.75,
                       labels = c("Trimmed", "Under sampled")) +
  labs(y = "density",
       x = "Estimated total N",
       fill = "Data",
       title = "Estimating total N based on \u03bb estimates.",
       caption = "Distribution of estimates of the total N in a community when using \u03bb estimates from under sampled (green) and trimmed (pruple) data. \nPanels showing increasing effect of undersampling from top left to bottom right. \nI.e., the top left panel has very little under sampling occuring and the bottom right panel is extremely undersampled. \nSolid vertical line shows the original N = 5000. This plot shows results when \u03bb = -2")
ggsave("plots/back_est_n_l_2_N5000.png", scale = 2)

# just the n = 5000, lambda = -1.5 
back_dat |>
  filter(original_n == 5000,
         lambda == -1.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75) +
  labs(title = "\u03bb = -1.5; N = 5000")

# just the n = 5000, lambda = -2.5 
back_dat |>
  filter(original_n == 5000,
         lambda == -2.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75)+
  labs(title = "\u03bb = -2.5; N = 5000")

# the other n's ####
# N = 1000 ####
back_dat |>
  filter(original_n == 1000,
         lambda == -2) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05,
                       end = 0.75,
                       labels = c("Trimmed", "Under sampled")) +
  labs(y = "density",
       x = "Estimated total N",
       fill = "Data",
       title = "\u03bb = -2; N = 1000")
# just the n = 5000, lambda = -1.5 
back_dat |>
  filter(original_n == 1000,
         lambda == -1.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75) +
  labs(title = "\u03bb = -1.5; N = 1000")

# just the n = 1000, lambda = -2.5 
back_dat |>
  filter(original_n == 1000,
         lambda == -2.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75)+
  labs(title = "\u03bb = -2.5; N = 1000")


# N = 10000 ####
back_dat |>
  filter(original_n == 10000,
         lambda == -2) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05,
                       end = 0.75,
                       labels = c("Trimmed", "Under sampled")) +
  labs(y = "density",
       x = "Estimated total N",
       fill = "Data",
       title = "\u03bb = -2; N = 10000")
# just the n = 10000, lambda = -1.5 
back_dat |>
  filter(original_n == 10000,
         lambda == -1.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75) +
  labs(title = "\u03bb = -1.5; N = 10000")

# just the n = 10000, lambda = -2.5 
back_dat |>
  filter(original_n == 10000,
         lambda == -2.5) |>
  pivot_longer(est_under_n:est_trimmed_n) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  geom_vline(aes(xintercept = original_n)) +
  facet_wrap(.~pr_scenario, 
             scales = "free", 
             ncol = 2) +
  theme_bw() +
  scale_fill_viridis_d(begin = 0.05, end = 0.75)+
  labs(title = "\u03bb = -2.5; N = 10000")
