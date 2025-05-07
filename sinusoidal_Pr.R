# function for sinusoidal response of sampling probabilility


sample_pr <- function(x, h, b){
  pr = 1 / (1 + (h / x**b))
  return(pr)
}

x = seq(0.01, 1, length.out = 100)
h = seq(0.01, 1, length.out = 10)
b = seq(1, 20, by = 1)
pr = sample_pr(x = x,
          h = 0.01, 
          b = 5)
dat <- data.frame(x = x, pr = pr)
plot(pr~x, data = dat)

library(tidyverse)

dat2 <- expand_grid(h = h,
                    b = b,
                    x = x)
dat2 |>
  mutate(pr = sample_pr(x = x, h = h, b = b)) |>
  ggplot(aes(x = x, 
             y = pr, 
             color = (h))) +
  geom_point() +
  facet_wrap(~b)
