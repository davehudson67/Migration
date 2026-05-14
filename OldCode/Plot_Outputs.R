# ── Packages ───────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(purrr)
library(scales)

# ── Small helpers ──────────────────────────────────────────────────────────────
invlogit <- function(x) 1/(1+exp(-x))

get_terms_rhs <- function(fit, dat) {
  stopifnot(!is.null(fit$call$formula))
  fml   <- as.formula(fit$call$formula)
  Terms <- terms(fml, data = dat)
  delete.response(Terms)
}

get_boot_betas <- function(fit) {
  beta_names <- names(coef(fit))               # fixed effects only
  B_raw <- as.matrix(fit$bootstrap)
  if (nrow(B_raw) < ncol(B_raw)) B_raw <- t(B_raw)
  if (!all(beta_names %in% colnames(B_raw))) {
    stop("Bootstrap columns don't contain all coefficient names in coef(fit).")
  }
  B_raw[, beta_names, drop = FALSE]            # n_boot x p (drops 'alpha')
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
    mm <- cbind(mm, matrix(0, nrow = nrow(mm), ncol = length(missing),
                           dimnames = list(NULL, missing)))
  }
  mm[, beta_names, drop = FALSE]
}

# Build a grid only for variables that exist on the RHS (EOO included if present)
make_grid <- function(Terms_rhs, dat, eoo = 0) {
  rhs_vars <- all.vars(Terms_rhs)
  
  if (!"HemFly" %in% rhs_vars) {
    stop("RHS must contain 'HemFly' for these plots.")
  }
  if (!is.factor(dat$HemFly)) dat$HemFly <- factor(dat$HemFly)
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

# Get *draws* for every Region × Migratory at an EOO value ──
# scale = "prob" (default) or "logit"
predict_draws_grid <- function(fit, dat, eoo = 0, scale = c("prob","logit")) {
  scale <- match.arg(scale)
  Terms_rhs <- get_terms_rhs(fit, dat)
  B         <- get_boot_betas(fit)
  beta_names <- colnames(B)
  contrs    <- get_contrasts(Terms_rhs, dat)
  
  grid <- make_grid(Terms_rhs, dat, eoo = eoo)
  mm   <- build_mm(Terms_rhs, grid, contrs, beta_names)
  
  eta <- B %*% t(mm)                  # n_boot x n_rows
  D   <- if (scale == "prob") invlogit(eta) else eta
  
  list(draws = D, grid = grid, scale = scale)
}

# Summaries for Region × Migratory (paired medians & 95% CIs) ───────────────
summarise_pairs <- function(pred, nice_labels = TRUE) {
  D <- pred$draws; grid <- pred$grid; scale <- pred$scale
  
  qs <- function(v) quantile(v, c(.025, .5, .975), na.rm = TRUE)
  S  <- apply(D, 2, qs)  # 3 x n_rows
  
  out <- cbind(grid, t(S)) %>%
    as.data.frame() %>%
    rename(lo95 = `2.5%`, med = `50%`, hi95 = `97.5%`) %>%
    mutate(
      Region     = HemFly,
      Hemisphere = substr(as.character(HemFly), 1, 1),
      Flyway = substr(as.character(HemFly), 3, nchar(as.character(HemFly)))
    )
  
  # Nice labels for migratory
  if ("migratory" %in% names(out) && nice_labels) {
    out <- out %>%
      mutate(
        Mig = ifelse(migratory == 1, "Migrant", "Resident"),
        Mig = factor(Mig, levels = c("Resident","Migrant"))
      )
  }
  
  attr(out, "scale") <- scale
  out
}

# ── Paired bootstrap difference (Migrant − Resident) on chosen scale ──────────
summarise_mig_effect_paired <- function(pred) {
  D <- pred$draws; grid <- pred$grid; scale <- pred$scale
  if (!"migratory" %in% names(grid)) {
    stop("Model has no 'migratory' term; cannot form Migrant − Resident contrast.")
  }
  
  # Map resident/migrant rows per Region (assumes same HemFly order)
  if (is.factor(grid$migratory)) {
    is_res <- as.character(grid$migratory) == levels(grid$migratory)[1]
    is_mig <- as.character(grid$migratory) == levels(grid$migratory)[2]
  } else {
    is_res <- grid$migratory == 0
    is_mig <- grid$migratory == 1
  }
  
  # Ensure we have one resident and one migrant row per Region in the same order
  stopifnot(all(grid$HemFly[is_res] == grid$HemFly[is_mig]))
  
  # Paired differences per bootstrap draw and region
  d <- D[, is_mig, drop = FALSE] - D[, is_res, drop = FALSE]   # n_boot x n_regions
  
  qs <- function(v) stats::quantile(v, c(.025, .5, .975), na.rm = TRUE)
  S  <- apply(d, 2, qs)
  
  tibble(
    Region   = grid$HemFly[is_res],
    diff_lo95 = S[1,],
    diff_med  = S[2,],
    diff_hi95 = S[3,],
    .scale    = scale
  )
}

# ── Plots (auto-format axes based on scale) ────────────────────────────────────
plot_paired <- function(summ_df) {
  scale <- attr(summ_df, "scale") %||% "prob"
  x_lab <- if (scale == "prob") "Pr(decline)" else "Linear predictor (log-odds)"
  x_fmt <- if (scale == "prob") label_percent(accuracy = 1) else identity
  
  ggplot(summ_df, aes(x = med, y = Region, colour = Mig)) +
    geom_linerange(aes(xmin = lo95, xmax = hi95),
                   position = position_dodge(width = 0.6),
                   linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.6), size = 2.7) +
    scale_x_continuous(labels = x_fmt) +
    labs(x = x_lab, y = "Region (HemFly)",
         title = if (scale == "prob") "Pr(decline) by region and migratory"
         else "Log-odds by region and migratory",
         subtitle = "Points = medians; bars = 95% paired bootstrap CIs") +
    theme_bw(base_size = 12)
}

plot_mig_effect_paired <- function(diff_df) {
  scale <- unique(diff_df$.scale)
  x_lab <- if (scale == "prob") "Migrant − Resident (probability)"
  else "Migrant − Resident (log-odds)"
  x_fmt <- if (scale == "prob") label_percent(accuracy = 1) else identity
  
  ggplot(diff_df, aes(x = diff_med, y = Region)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_linerange(aes(xmin = diff_lo95, xmax = diff_hi95), linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_x_continuous(labels = x_fmt) +
    labs(x = x_lab, y = "Region (HemFly)",
         title = if (scale == "prob")
           "Migration effect on Pr(decline) by region"
         else
           "Migration effect on log-odds by region",
         subtitle = "Paired bootstrap medians and 95% CIs") +
    theme_minimal(base_size = 12)
}

# ── Usage ─────────────────────────────────────────────────────────────────────
# Load your objects
fit <- readRDS("Outputs/fit_phyloHemFly_1AmOnly_1000_scaledT.rds")
summary(fit)
fit <- readRDS("Outputs/fit_phyloHemFly_1AmOnly_1000.rds")
summary(fit)
data <- readRDS("Data_for_plots.rds")
levels(data$HemFly)
data$HemFly <- factor(data$HemFly, levels = c("1_Am_only", 
                                              "1_Af_only",
                                              "1_As_only",
                                              "1_Af_As",
                                              "1_Am_As",
                                              "1_Am_Af_As",
                                              "2_Am_only",
                                              "2_Af_only",
                                              "2_As_only",
                                              "3_Am_only",
                                              "3_Af_only",
                                              "3_As_only",
                                              "3_Af_As",
                                              "3_Am_As",
                                              "3_Am_Af_As"))

## A) Probability scale (default)
pred_prob <- predict_draws_grid(fit, data, eoo = 0, scale = "prob")
summ_prob <- summarise_pairs(pred_prob) %>% 
  mutate(Region = fct_inorder(Region))
diff_prob <- summarise_mig_effect_paired(pred_prob)

p1_prob <- plot_paired(summ_prob)
p2_prob <- plot_mig_effect_paired(diff_prob)
p1_prob
p2_prob

## B) Log-odds (linear predictor) scale
pred_logit <- predict_draws_grid(fit, data, eoo = 0, scale = "logit")
summ_logit <- summarise_pairs(pred_logit) %>% mutate(Region = fct_inorder(Region))
diff_logit <- summarise_mig_effect_paired(pred_logit)

p1_logit <- plot_paired(summ_logit)
p2_logit <- plot_mig_effect_paired(diff_logit)
p1_logit
p2_logit

#-------------------------------------------------------------------------------
# Plot EOO and alpha

# Helper: pull a named bootstrap column safely (drops to numeric vector)
get_boot_col <- function(fit, colname) {
  bm <- as.matrix(fit$bootstrap)
  if (nrow(bm) < ncol(bm)) bm <- t(bm)   # ensure rows = draws
  stopifnot(!is.null(colnames(bm)))
  if (!colname %in% colnames(bm)) {
    stop(sprintf("Column '%s' not found in fit$bootstrap. Available: %s",
                 colname, paste(colnames(bm), collapse = ", ")))
  }
  as.numeric(bm[, colname])
}

# Extract draws
eoo_draws   <- get_boot_col(fit, "EOO_log_cent")  # log-odds slope per unit of centered log EOO
alpha_draws <- get_boot_col(fit, "alpha")         # phylogenetic signal

# Tidy long format
boot_df <- bind_rows(
#  tibble(Parameter = "EOO slope (β_EOO)", value = eoo_draws),
  tibble(Parameter = "Alpha (phylo signal)", value = alpha_draws)
)

# Plot
summ_df <- boot_df %>%
  group_by(Parameter) %>%
  summarise(
    med  = median(value, na.rm = TRUE),
    lo95 = quantile(value, 0.025, na.rm = TRUE),
    hi95 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# --- Violin plot with median point and 95% bar ---
p <- ggplot(boot_df, aes(x = Parameter, y = value)) +
  geom_violin(fill = "grey80", colour = "grey40", alpha = 0.6, trim = FALSE) +
  geom_point(data = summ_df, aes(x = Parameter, y = med), colour = "red", size = 2, inherit.aes = FALSE) +
  geom_errorbar(
    data = summ_df,
    aes(x = Parameter, y = med, ymin = lo95, ymax = hi95),
    width = 0.2, colour = "red",
    inherit.aes = FALSE
  )+
  coord_flip() +
  labs(x = NULL, y = "Value",
       title = "Bootstrap distributions",
       subtitle = "Violin = density; red point = median; red bar = 95% interval") +
  theme_bw(base_size = 13) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(~Parameter, scales = "free")

p

