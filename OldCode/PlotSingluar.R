library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# Load
fit <- readRDS("Outputs/fit_phyloHemFlyEOO_1AmOnly_1000.rds")
summary(fit)
boot_mat <- fit$bootstrap[,1:30]

# Label columns with coefficient names
beta_hat <- coef(fit)
pnames <- names(beta_hat)

# Build a tidy table of bootstrap draws
boot_long <- as.data.frame(boot_mat) |>
  mutate(.draw = row_number()) |>
  pivot_longer(-.draw, names_to = "term", values_to = "value")

# Compute CIs and medians
ci_tbl <- boot_long |>
  group_by(term) |>
  summarise(
    median = median(value, na.rm = TRUE),
    lo95   = quantile(value, 0.025, na.rm = TRUE),
    hi95   = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot coefficient bootstrap distributions with 95% CI and median
p_coefs <- ggplot(boot_long, aes(x = value, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_point(data = ci_tbl, aes(x = median, y = term), size = 1.8) +
  geom_errorbarh(data = ci_tbl, aes(xmin = lo95, xmax = hi95, y = term), height = 0.15, linewidth = 0.7, inherit.aes = FALSE) +
  theme_minimal(base_size = 12)

p_coefs

# Plot bootstrap for alpha
alpha_vals <- fit$bootstrap[, 31]
alpha_df <- tibble(alpha = as.numeric(alpha_vals))
alpha_ci <- alpha_df |>
    summarise(
      median = median(alpha), lo95 = quantile(alpha, 0.025), hi95 = quantile(alpha, 0.975)
    )
  
p_alpha <- ggplot(alpha_df, aes(x = alpha, y = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.9) +
  geom_point(data = alpha_ci, aes(y = 0, x = median), size = 2, inherit.aes = FALSE) +
  geom_errorbar(data = alpha_ci, aes(y = 0, xmin = lo95, xmax = hi95), width = 0.12, inherit.aes = FALSE) +
  scale_y_continuous(breaks = NULL) +
  labs(y = NULL, x = expression(alpha), title = "Bootstrap distribution and 95% CI for alpha") +
  theme_minimal(base_size = 12)

p_alpha
