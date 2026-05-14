# ======================================================================
# RIDGE PLOTS FROM SAVED PHYLOGLM RDS OBJECTS
# ======================================================================

rm(list = ls())

# ── Packages ────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(ggridges)
library(patchwork)
library(scales)

# ── Paths ───────────────────────────────────────────────────────────────
setwd("~/Migration")

out_dir <- "Outputs"

# ── Load saved objects ──────────────────────────────────────────────────
data <- readRDS(file.path(out_dir, "Data_for_plots_phyloglm.rds"))
fit_phylo_mig <- readRDS(file.path(out_dir, "fit_phylo_mig.rds"))
fit_phyloHemFly <- readRDS(file.path(out_dir, "fit_phyloHemFly.rds"))
fit_phyloHemFlyEOO <- readRDS(file.path(out_dir, "fit_phyloHemFlyEOO.rds"))

# ── Helpers ─────────────────────────────────────────────────────────────
invlogit <- function(x) 1 / (1 + exp(-x))

`%||%` <- function(x, y) if (is.null(x)) y else x

get_terms_rhs <- function(fit, dat) {
  stopifnot(!is.null(fit$call$formula))
  fml   <- as.formula(fit$call$formula)
  Terms <- terms(fml, data = dat)
  delete.response(Terms)
}

get_boot_betas <- function(fit) {
  beta_names <- names(coef(fit))
  B_raw <- as.matrix(fit$bootstrap)
  
  if (all(beta_names %in% colnames(B_raw))) {
    B <- B_raw
  } else if (all(beta_names %in% rownames(B_raw))) {
    B <- t(B_raw)
  } else {
    stop(
      "Could not align bootstrap matrix with coefficient names.\n",
      "coef names: ", paste(beta_names, collapse = ", "), "\n",
      "bootstrap colnames: ", paste(colnames(B_raw), collapse = ", "), "\n",
      "bootstrap rownames: ", paste(rownames(B_raw), collapse = ", ")
    )
  }
  
  if ("alpha" %in% colnames(B)) {
    B <- B[, colnames(B) != "alpha", drop = FALSE]
  }
  
  B[, beta_names, drop = FALSE]
}

get_contrasts <- function(Terms_rhs, dat) {
  mm_dat <- model.matrix(Terms_rhs, data = dat)
  attr(mm_dat, "contrasts")
}

build_mm <- function(Terms_rhs, newdata, contrs, beta_names) {
  mm <- model.matrix(Terms_rhs, data = newdata, contrasts.arg = contrs)
  mm <- mm[, colnames(mm) %in% beta_names, drop = FALSE]
  
  missing <- setdiff(beta_names, colnames(mm))
  if (length(missing)) {
    mm <- cbind(
      mm,
      matrix(
        0,
        nrow = nrow(mm),
        ncol = length(missing),
        dimnames = list(NULL, missing)
      )
    )
  }
  
  mm[, beta_names, drop = FALSE]
}

make_grid <- function(Terms_rhs, dat, eoo = 0) {
  rhs_vars <- all.vars(Terms_rhs)
  
  if (!"HemFly" %in% rhs_vars) {
    stop("RHS must contain 'HemFly' for these region plots.")
  }
  
  hf_levels <- levels(dat$HemFly)
  comps <- list(HemFly = factor(hf_levels, levels = hf_levels))
  
  if ("migratory" %in% rhs_vars) {
    if (is.factor(dat$migratory)) {
      comps$migratory <- factor(levels(dat$migratory), levels = levels(dat$migratory))
    } else {
      comps$migratory <- c(0, 1)
    }
  }
  
  if ("EOO_log_cent" %in% rhs_vars) {
    comps$EOO_log_cent <- eoo
  }
  
  do.call(tidyr::expand_grid, comps)
}

predict_draws_grid <- function(fit, dat, eoo = 0, scale = c("prob", "logit")) {
  scale <- match.arg(scale)
  
  Terms_rhs  <- get_terms_rhs(fit, dat)
  B          <- get_boot_betas(fit)
  beta_names <- colnames(B)
  contrs     <- get_contrasts(Terms_rhs, dat)
  
  grid <- make_grid(Terms_rhs, dat, eoo = eoo)
  mm   <- build_mm(Terms_rhs, grid, contrs, beta_names)
  
  eta <- B %*% t(mm)
  D   <- if (scale == "prob") invlogit(eta) else eta
  
  list(draws = D, grid = grid, scale = scale)
}

get_boot_col <- function(fit, colname) {
  bm <- as.matrix(fit$bootstrap)
  
  if (colname %in% colnames(bm)) {
    B <- bm
  } else if (colname %in% rownames(bm)) {
    B <- t(bm)
  } else {
    stop(sprintf("Column '%s' not found in fit$bootstrap.", colname))
  }
  
  as.numeric(B[, colname])
}

draws_to_long <- function(draws, group_names, value_name = "value", group_var = "group") {
  as.data.frame(draws) %>%
    setNames(group_names) %>%
    tibble::rowid_to_column("draw") %>%
    tidyr::pivot_longer(
      cols = -draw,
      names_to = group_var,
      values_to = value_name
    )
}

predict_draws_long <- function(fit, dat, eoo = 0, scale = c("prob", "logit")) {
  pred <- predict_draws_grid(fit, dat, eoo = eoo, scale = scale)
  D <- pred$draws
  grid <- pred$grid
  scale <- pred$scale
  
  long_df <- as.data.frame(t(D))
  names(long_df) <- paste0("draw_", seq_len(ncol(long_df)))
  
  bind_cols(grid, long_df) %>%
    tidyr::pivot_longer(
      cols = starts_with("draw_"),
      names_to = "draw",
      values_to = "value"
    ) %>%
    mutate(
      Region = HemFly,
      Hemisphere = substr(as.character(HemFly), 1, 1),
      Flyway = substr(as.character(HemFly), 3, nchar(as.character(HemFly))),
      Mig = ifelse(as.character(migratory) == levels(dat$migratory)[2], "Migrant", "Resident"),
      Mig = factor(Mig, levels = c("Resident", "Migrant")),
      .scale = scale
    )
}

predict_mig_diff_long <- function(fit, dat, eoo = 0, scale = c("prob", "logit")) {
  pred <- predict_draws_grid(fit, dat, eoo = eoo, scale = scale)
  D <- pred$draws
  grid <- pred$grid
  scale <- pred$scale
  
  if (is.factor(grid$migratory)) {
    is_res <- as.character(grid$migratory) == levels(grid$migratory)[1]
    is_mig <- as.character(grid$migratory) == levels(grid$migratory)[2]
  } else {
    is_res <- grid$migratory == 0
    is_mig <- grid$migratory == 1
  }
  
  stopifnot(all(grid$HemFly[is_res] == grid$HemFly[is_mig]))
  
  d <- D[, is_mig, drop = FALSE] - D[, is_res, drop = FALSE]
  
  diff_df <- as.data.frame(d)
  names(diff_df) <- as.character(grid$HemFly[is_res])
  
  diff_df %>%
    tibble::rowid_to_column("draw") %>%
    tidyr::pivot_longer(
      cols = -draw,
      names_to = "Region",
      values_to = "value"
    ) %>%
    mutate(
      Region = factor(Region, levels = levels(dat$HemFly)),
      .scale = scale
    )
}

# ── Figure 1A: violin for migratory effect ─────────────────────────────
boot_mig <- as.matrix(fit_phylo_mig$bootstrap)
if (all(names(coef(fit_phylo_mig)) %in% rownames(boot_mig))) {
  boot_mig <- t(boot_mig)
}

mig_coef_name <- grep("^migratory", colnames(boot_mig), value = TRUE)[1]
if (is.na(mig_coef_name)) stop("Could not find migratory coefficient in bootstrap matrix.")

ridge_F1a <- tibble(
  group = factor("Migratory effect", levels = "Migratory effect"),
  value = boot_mig[, mig_coef_name]
)

ci_F1a <- ridge_F1a %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE)
  )

p_F1a <- ggplot(ridge_F1a, aes(x = group, y = value)) +
  geom_violin(fill = "#40B0A6", alpha = 0.7, colour = NA, width = 0.8) +
  geom_errorbar(
    data = ci_F1a,
    aes(x = 1, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    data = ci_F1a,
    aes(x = 1, y = med),
    inherit.aes = FALSE,
    size = 2,
    colour = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
  coord_flip(ylim = c(min(ridge_F1a$value) - 0.05, 0.05)) +
  labs(
    x = NULL,
    y = "Log-odds parameter estimate",
    title = "Migratory effect (phyloglm bootstrap)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

p_F1a

# ── Figure 1B: violins for predicted probabilities ─────────────────────
newdata_simple <- data.frame(
  migratory = factor(levels(data$migratory), levels = levels(data$migratory))
)

Terms_simple <- get_terms_rhs(fit_phylo_mig, data)
B_simple <- get_boot_betas(fit_phylo_mig)
contrs_simple <- get_contrasts(Terms_simple, data)
mm_simple <- build_mm(Terms_simple, newdata_simple, contrs_simple, colnames(B_simple))
eta_simple <- B_simple %*% t(mm_simple)
prob_simple <- invlogit(eta_simple)

plot_F1b_df <- draws_to_long(
  draws = prob_simple,
  group_names = c("Resident", "Migrant"),
  value_name = "value",
  group_var = "group"
) %>%
  mutate(group = factor(group, levels = c("Resident", "Migrant")))

ci_F1b <- plot_F1b_df %>%
  group_by(group) %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_F1b <- ggplot(plot_F1b_df, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.75, colour = NA, width = 0.8) +
  geom_errorbar(
    data = ci_F1b,
    aes(x = group, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    data = ci_F1b,
    aes(x = group, y = med),
    inherit.aes = FALSE,
    size = 2,
    colour = "black"
  ) +
  scale_fill_manual(values = c("Resident" = "#E1BE6A", "Migrant" = "#40B0A6")) +
  coord_flip(ylim = c(0.4, 0.8)) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Probability of decline",
    title = "Predicted probability of decline"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

fig1 <- p_F1a | p_F1b
fig1

# ── Plot 2A: ridge plot of HemFly × migratory probabilities ────────────
ridge_region_prob <- predict_draws_long(
  fit_phyloHemFlyEOO,
  data,
  eoo = 0,
  scale = "prob"
) %>%
  mutate(
    Region = fct_rev(factor(Region, levels = levels(data$HemFly))),
    Mig = factor(Mig, levels = c("Resident", "Migrant"))
  )

ci_region_prob <- ridge_region_prob %>%
  group_by(Region, Mig) %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y = as.numeric(Region) - 0.10   # both on same line
  )

p_region_prob <- ggplot(
  ridge_region_prob,
  aes(x = value, y = Region, fill = Mig)
) +
  geom_density_ridges(
    alpha = 0.4,
    colour = NA,
    scale = 0.9,
    rel_min_height = 0.001,
    position = "identity"
  ) +
  geom_segment(
    data = ci_region_prob,
    aes(x = lo95, xend = hi95, y = y, yend = y, colour = Mig),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_point(
    data = ci_region_prob,
    aes(x = med, y = y, colour = Mig),
    inherit.aes = FALSE,
    size = 1.5
  ) +
  scale_fill_manual(values = c("Resident" = "#F8766D", "Migrant" = "#00BFC4")) +
  scale_colour_manual(values = c("Resident" = "#F8766D", "Migrant" = "#00BFC4")) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Pr(decline)",
    y = "Region (HemFly)",
    title = "Pr(decline) by region and migratory status",
    subtitle = "Bootstrap distributions"
  ) +
  theme_bw(base_size = 12)

p_region_prob

# ── Plot 2B: ridge plot of migrant − resident differences ──────────────
ridge_region_diff <- predict_mig_diff_long(
  fit_phyloHemFlyEOO,
  data,
  eoo = 0,
  scale = "prob"
) %>%
  mutate(Region = fct_rev(factor(Region, levels = levels(data$HemFly))))

ci_region_diff <- ridge_region_diff %>%
  group_by(Region) %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_region_diff <- ggplot(
  ridge_region_diff,
  aes(x = value, y = Region)
) +
  geom_density_ridges(
    fill = "grey70",
    alpha = 0.8,
    colour = NA,
    scale = 0.9,
    rel_min_height = 0.001
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(
    data = ci_region_diff,
    aes(x = lo95, xend = hi95, y = Region, yend = Region),
    inherit.aes = FALSE,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    data = ci_region_diff,
    aes(x = med, y = Region),
    inherit.aes = FALSE,
    size = 1.5,
    colour = "black"
  ) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Migrant − Resident (probability)",
    y = "Region (HemFly)",
    title = "Migration effect on Pr(decline) by region",
    subtitle = "Bootstrap distributions"
  ) +
  theme_bw(base_size = 12)

p_region_diff

# ── Plot 3: violins for EOO slope and alpha ────────────────────────────

boot_eoo <- tibble(
  group = factor("EOO slope", levels = "EOO slope"),
  value = get_boot_col(fit_phyloHemFlyEOO, "EOO_log_cent")
)

boot_alpha <- tibble(
  group = factor("Alpha (phylo signal)", levels = "Alpha (phylo signal)"),
  value = get_boot_col(fit_phyloHemFlyEOO, "alpha")
)

ci_eoo <- boot_eoo %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE)
  )

ci_alpha <- boot_alpha %>%
  summarise(
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    med  = median(value, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE)
  )

p_eoo <- ggplot(boot_eoo, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.75, colour = NA, width = 0.8) +
  geom_errorbar(
    data = ci_eoo,
    aes(x = 1, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    data = ci_eoo,
    aes(x = 1, y = med),
    inherit.aes = FALSE,
    size = 2,
    colour = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("EOO slope" = "#40B0A6")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Effect Size"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

p_alpha <- ggplot(boot_alpha, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.75, colour = NA, width = 0.8) +
  geom_errorbar(
    data = ci_alpha,
    aes(x = 1, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    data = ci_alpha,
    aes(x = 1, y = med),
    inherit.aes = FALSE,
    size = 2,
    colour = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("Alpha (phylo signal)" = "#40B0A6")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Estimate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

p_boot <- p_eoo / p_alpha

# ── Combine figures ─────────────────────────────────────────────────────
fig1 <- p_F1a | p_F1b
fig2 <- p_region_prob | p_region_diff

# ── Print ───────────────────────────────────────────────────────────────
print(fig1)
print(fig2)
print(p_boot)

# ── Save ────────────────────────────────────────────────────────────────
ggsave(file.path(out_dir, "Fig1_phyloglm_migratory_ridges.png"),
       fig1, width = 12, height = 6, dpi = 300)

ggsave(file.path(out_dir, "Fig2_phyloglm_HemFly_ridges.png"),
       fig2, width = 14, height = 8, dpi = 300)

ggsave(file.path(out_dir, "Fig3_phyloglm_bootstrap_ridges.png"),
       p_boot, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------
