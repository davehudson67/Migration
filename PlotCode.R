library(ape)
library(dplyr)
library(phylolm)
library(phytools)
library(ape)
library(nimble)
library(MCMCglmm)
library(Matrix)
library(tidyverse)
library(coda)
library(mcmcplots)
library(patchwork)


# 1. Read and extract everything, including SEs
files <- list.files(
  path       = "Outputs",
  pattern    = "^FitPhyloHemFlyEOO_10_ref_.*\\.rds$",
  full.names = TRUE
)

results <- lapply(files, function(f) {
#  browser()
  ref <- sub(".*_ref_(.*)\\.rds$", "\\1", basename(f))
  fit <- readRDS(f)
  summ <- summary(fit)
  
  # your existing coefficient bits
  s       <- summ$coefficients
  migr    <- s["migratory1", "Estimate"]
  se_m    <- s["migratory1", 2]
  #migr_lo <- s["migratory1", 4]
  migr_lo <- migr - 1.96 * se_m
  #migr_hi <- s["migratory1", 5]
  migr_hi <- migr + 1.96 * se_m
  
  eoo     <- s["EOO_log_cent", "Estimate"]
  se_e    <- s["EOO_log_cent", 2]
  #eoo_lo  <- s["EOO_log_cent", 4]
  eoo_lo <- eoo - 1.96 * se_e
  #eoo_hi  <- s["EOO_log_cent", 5]
  eoo_hi <- eoo + 1.96 * se_e
  
  # alpha
  #alpha    <- summ$alpha 
  #alpha_lo <- summ$bootconfint95[1, 32]   # the lower 95% CI
  #alpha_hi <- summ$bootconfint95[2, 32]   # the upper 95% CI
  
  hemisphere <- sub("^(\\d+)_.*$", "\\1", ref)
  flyway     <- sub("^\\d+_(.*)$", "\\1", ref)
  
  data_frame(
    ref, hemisphere, flyway,
    migr,   se_m,   migr_lo,   migr_hi,
    eoo,    se_e,   eoo_lo,    eoo_hi
    #alpha, alpha_lo, alpha_hi
  )
})

df <- do.call(rbind, results)
df$hemisphere <- factor(df$hemisphere, levels = unique(df$hemisphere))
df$flyway     <- factor(df$flyway,     levels = unique(df$flyway))
df$ref        <- factor(df$ref,        levels = unique(df$ref))

# 2. Plot migratory1 with 95% CI
m <- ggplot(df, aes(x = ref, y = migr, color = flyway)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = migr_lo, ymax = migr_hi),  
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ hemisphere, scales = "free_x") +
  labs(
    x      = "Reference HemFly level",
    y      = "Migratory effect (coef ± 95% CI)",
    colour = "Flyway",
    title  = "Migratory effect by HemFly (with 95% CI)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Plot EOO_log_cent effect with 95% CI
e <- ggplot(df, aes(x = ref, y = eoo, color = flyway)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = eoo_lo, ymax = eoo_hi),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ hemisphere, scales = "free_x") +
  labs(
    x      = "Reference HemFly level",
    y      = "EOO effect (coef ± 95% CI)",
    colour = "Flyway",
    title  = "Effect of EOO_log_cent by HemFly (with 95% CI)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Plot alpha effect with 95% CI
a <- ggplot(df, aes(x = ref, y = alpha, color = flyway)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = alpha_lo, ymax = alpha_hi),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ hemisphere, scales = "free_x") +
  labs(
    x      = "Reference HemFly level",
    y      = "alpha (coef ± 95% CI)",
    colour = "Flyway",
    title  = "alpha by HemFly (with 95% CI)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

m
e
a

m/e + plot_layout(guides = 'collect')
m


