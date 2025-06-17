# summarize undersampling simulation

library(tidyverse)
library(tidybayes)

sim <- as_tibble(readRDS("simulation_results/undersampling_sim_run.rds"))

sim |>
  group_by(rep) |>
  count() |>
  ungroup() |>
  group_by(n) |>
  count()
# 36 original parameter sets and 500 reps = 18,000
# 26 parameter sets are missing a single rep
# 17,974 / 18,000 = 99.86% of reps successful

sim |>
  summarize(min_xmax_obs = min(xmax_obs),
            mean_xmax_obs = mean(xmax_obs),
            med_xmax_obs = median(xmax_obs),
            max_xmax_obs = max(xmax_obs))
# xmax_obs has a lot of variation

sim |>
  summarize(min_xmin_obs = min(xmin_obs),
            mean_xmin_obs = mean(xmin_obs),
            med_xmin_obs = median(xmin_obs),
            max_xmin_obs = max(xmin_obs))
# xmin_obs is very consistent at the set boundary of 0.001
  
ggplot(sim, 
       aes(x = original_n,
           y = under_n)) +
  geom_point(position = position_jitter(width = 500))

# distribution of lambda estimates (not including CIs)
sim |>
  as_tibble() |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  select(b, h, known_lambda, lambda_under, lambda_trimmed, original_n) |>
  rename(`under sampled` = lambda_under,
         `trimmed data` = lambda_trimmed) |>
  pivot_longer(`under sampled`:`trimmed data`) |>
  ggplot(aes(x = value, 
             fill = name)) +
  stat_halfeye() +
  geom_vline(aes(xintercept = known_lambda),
             linetype = "dashed") +
  scale_fill_manual(values = c(c("#FF1984",
                                 "#019AFF"))) +
  facet_wrap(b ~ h,
             #scales = "free_x",
             labeller = label_both) +
  theme_bw() +
  labs(x = "\u03bb estimate",
       y = "density",
       caption = "Distribution of \u03bb estimates for under sampled (blue) and data which has been trimmed at the estimated x_min (pink). \nFacet titles show variables for undersampling function and the number of data points in the under sampled data \nincreases from left to right and top to bottom. Dashed line shows the known value of \u03bb")
ggsave("plots/under_trimmed_estimates_lambda_2.png", 
       scale = 2)
  

# proportion of samples
sim |>
  select(h, b, original_n, under_n) |>
  mutate(prop = under_n / original_n) |>
  group_by(h, b) |>
  summarize(med_prop = median(prop),
            sd_prop = sd(prop))
# h == 0.0001 
  # b == 1.5 ~ 89% +- 2.7
  # b == 2 ~ 40% +- 12.2
# h == 0.001 
  # b == 1.5 ~53% +- 9.6
  # b == 2 ~15% +- 11.5
sim |>
  select(h, b, original_n, under_n) |>
  mutate(prop = under_n / original_n) |>
  group_by(h, b, original_n) |>
  summarize(med_prop = median(prop),
            mean_prop = mean(prop), 
            sd_prop = sd(prop))|>
  write_csv("simulation_summaries/proportion_samples_pr_scenarios.csv")

set.seed(111)
sim |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  group_by(pr_scenario) |>
  sample_n(50) |>
  arrange(lambda_under_lo) |>
  group_by(known_lambda, lambda_under_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_under_lo < known_lambda & lambda_under_hi > known_lambda) |>
  ggplot(aes(y = id,
             x = lambda_under,
             xmin = lambda_under_lo,
             xmax = lambda_under_hi,
             color = color)) +
  geom_pointrange(
    position = position_jitter(height = 0.05)) +
  #geom_vline(aes(xintercept = known_lambda)) +
  facet_wrap(b~h,
             scales = "free",
             labeller = label_both) +
  theme_bw() +
  scale_color_viridis_d(option = "plasma") +
  guides(
    color = guide_legend(
      title= "\u03bb in 95% CI?")) +
  labs(x = "95% CI for \u03bb estimate",
       #y = "",
       title = "95% CIs for undersampled data",
       caption = "95% CI's for under sampled data. Color corresponds to whether or not the CI has the true \u03bb in it. This plot is for \u03bb = -2. and original N = 5000 \n Only 50 CIs are shown for visualization") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("plots/ci_under_lambda_2.png", width = 9, height = 4.5)

# CI plot color if lambda is in or out
set.seed(111)
sim |>
  filter(original_n == 5000,
         known_lambda == -2) |>
  group_by(pr_scenario) |>
  sample_n(50) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ggplot(aes(y = id,
             x = lambda_trimmed,
             xmin = lambda_trimmed_lo,
             xmax = lambda_trimmed_hi,
             color = color)) +
  geom_pointrange(position = position_jitter(height = 0.05)) +
  geom_vline(aes(xintercept = known_lambda)) +
  facet_wrap(b~h,
             scales = "free",
             labeller = label_both)+
  theme_bw() +
  scale_color_viridis_d(option = "plasma", end = 0.75) +
  guides(
    color = guide_legend(
      title= "\u03bb in 95% CI?")) +
  labs(x = "95% CI for \u03bb estimate",
       title = "95% CIs for Trimmed data",
       caption = "95% CI's for data which has been trimmed at the estimated xmin vlaue. Color corresponds to whether or not the CI has the true \u03bb in it. \nThis plot is for \u03bb = -2 and original N = 5000. Only 50 CIs are shown for visualization") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("plots/ci_trimmed_lambda_2.png",
       width = 9, height = 4.5)

# proportion of trimmed CIs with known lambda
sim |>
  filter(original_n == 5000) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ungroup() |>
  group_by(known_lambda, h, b, original_n) |>
  summarize(ci_true = sum(color) / n()) |>
  print(n = 36)

# plot of proportion of ci with lambda
sim |>
  #filter(original_n == 5000) |>
  arrange(lambda_trimmed_lo) |>
  group_by(known_lambda, h, lambda_trimmed_lo) |>
  mutate(id = cur_group_id(),
         color = lambda_trimmed_lo < known_lambda & lambda_trimmed_hi > known_lambda) |>
  ungroup() |>
  group_by(known_lambda, h, b, original_n) |>
  summarize(ci_true = mean(color),
            ci_sd = sd(color)) |>
  ggplot(aes(x = original_n,
             y = ci_true, 
             ymin = ci_true - ci_sd, 
             ymax = ci_true + ci_sd,
             color = interaction(h, b))) +
  geom_pointrange(
    position = position_dodge(width = 1000)
  )+
  facet_wrap(.~known_lambda) +
  scale_color_viridis_d(option = "plasma")


