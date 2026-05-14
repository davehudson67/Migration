# ======================================================================
# COMBINED ANALYSIS + PLOTS FOR PHYLOGLM MODELS
# Replaces old MCMCglmm-dependent plotting workflow
# ======================================================================

rm(list = ls())

# ── Packages ────────────────────────────────────────────────────────────
library(ape)
library(phytools)
library(phylolm)
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)
library(future)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)

# ── Paths ───────────────────────────────────────────────────────────────
# Edit these if needed
setwd("~/Migration")

data_file <- "Data/IUCNdata_updated_Taxonomy_030925.csv"
tree_file <- "Data/NewFixedTree_030925.rds"
out_dir   <- "Outputs"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Small helpers ──────────────────────────────────────────────────────
invlogit <- function(x) 1 / (1 + exp(-x))

`%||%` <- function(x, y) if (is.null(x)) y else x

# Build quantile "block" data for interval-rectangle plots
make_block_df <- function(draws, group_labels) {
  out <- lapply(seq_along(group_labels), function(i) {
    x <- draws[, i]
    qs <- quantile(x, probs = c(0.025, 0.25, 0.45, 0.55, 0.75, 0.975), na.rm = TRUE)
    tibble(
      group = group_labels[i],
      q025 = qs[1],
      q25  = qs[2],
      q45  = qs[3],
      q55  = qs[4],
      q75  = qs[5],
      q975 = qs[6],
      med  = median(x, na.rm = TRUE)
    )
  })
  bind_rows(out)
}

plot_interval_blocks <- function(block_df, fills, xlab, xintercept = NULL,
                                 title = NULL, percent_scale = FALSE) {
  block_df <- block_df %>%
    mutate(y = rev(seq_len(n())))
  
  p <- ggplot(block_df) +
    geom_rect(aes(xmin = q025, xmax = q975, ymin = y - 0.45, ymax = y + 0.45, fill = group),
              alpha = 0.25, colour = NA) +
    geom_rect(aes(xmin = q25, xmax = q75, ymin = y - 0.45, ymax = y + 0.45, fill = group),
              alpha = 0.50, colour = NA) +
    geom_rect(aes(xmin = q45, xmax = q55, ymin = y - 0.45, ymax = y + 0.45, fill = group),
              alpha = 1.00, colour = NA) +
    geom_segment(aes(x = med, xend = med, y = y - 0.45, yend = y + 0.45), linewidth = 1.2) +
    scale_y_continuous(breaks = block_df$y, labels = block_df$group) +
    scale_fill_manual(values = fills) +
    labs(x = xlab, y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(xintercept)) {
    p <- p + geom_vline(xintercept = xintercept, linetype = "dashed", alpha = 0.7)
  }
  
  if (percent_scale) {
    p <- p + scale_x_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1))
  }
  
  p
}

get_terms_rhs <- function(fit, dat) {
  stopifnot(!is.null(fit$call$formula))
  fml   <- as.formula(fit$call$formula)
  Terms <- terms(fml, data = dat)
  delete.response(Terms)
}

get_boot_betas <- function(fit) {
  beta_names <- names(coef(fit))
  B_raw <- as.matrix(fit$bootstrap)
  
  # Case 1: already correct orientation
  if (all(beta_names %in% colnames(B_raw))) {
    B <- B_raw
    # Case 2: transposed
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
  
  # drop alpha if present
  if ("alpha" %in% colnames(B)) {
    B <- B[, colnames(B) != "alpha", drop = FALSE]
  }
  
  # reorder to exactly match coef(fit)
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
    mm <- cbind(mm,
                matrix(0, nrow = nrow(mm), ncol = length(missing),
                       dimnames = list(NULL, missing)))
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

summarise_pairs <- function(pred, nice_labels = TRUE) {
  D <- pred$draws
  grid <- pred$grid
  scale <- pred$scale
  
  qs <- function(v) quantile(v, c(.025, .5, .975), na.rm = TRUE)
  S  <- apply(D, 2, qs)
  
  out <- cbind(grid, t(S)) %>%
    as.data.frame() %>%
    rename(lo95 = `2.5%`, med = `50%`, hi95 = `97.5%`) %>%
    mutate(
      Region     = HemFly,
      Hemisphere = substr(as.character(HemFly), 1, 1),
      Flyway     = substr(as.character(HemFly), 3, nchar(as.character(HemFly)))
    )
  
  if ("migratory" %in% names(out) && nice_labels) {
    out <- out %>%
      mutate(
        Mig = ifelse(as.character(migratory) == levels(migratory)[2], "Migrant", "Resident"),
        Mig = factor(Mig, levels = c("Resident", "Migrant"))
      )
  }
  
  attr(out, "scale") <- scale
  out
}

summarise_mig_effect_paired <- function(pred) {
  D <- pred$draws
  grid <- pred$grid
  scale <- pred$scale
  
  if (!"migratory" %in% names(grid)) {
    stop("Model has no migratory term.")
  }
  
  if (is.factor(grid$migratory)) {
    is_res <- as.character(grid$migratory) == levels(grid$migratory)[1]
    is_mig <- as.character(grid$migratory) == levels(grid$migratory)[2]
  } else {
    is_res <- grid$migratory == 0
    is_mig <- grid$migratory == 1
  }
  
  stopifnot(all(grid$HemFly[is_res] == grid$HemFly[is_mig]))
  
  d <- D[, is_mig, drop = FALSE] - D[, is_res, drop = FALSE]
  
  qs <- function(v) quantile(v, c(.025, .5, .975), na.rm = TRUE)
  S  <- apply(d, 2, qs)
  
  tibble(
    Region    = grid$HemFly[is_res],
    diff_lo95 = S[1, ],
    diff_med  = S[2, ],
    diff_hi95 = S[3, ],
    .scale    = scale
  )
}

plot_paired <- function(summ_df) {
  scale <- attr(summ_df, "scale") %||% "prob"
  x_lab <- if (scale == "prob") "Pr(decline)" else "Linear predictor (log-odds)"
  x_fmt <- if (scale == "prob") label_percent(accuracy = 1) else identity
  
  ggplot(summ_df, aes(x = med, y = Region, colour = Mig)) +
    geom_linerange(aes(xmin = lo95, xmax = hi95),
                   position = position_dodge(width = 0.6),
                   linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.6), size = 2.6) +
    scale_x_continuous(labels = x_fmt) +
    labs(
      x = x_lab,
      y = "Region (HemFly)",
      title = if (scale == "prob") "Pr(decline) by region and migratory status"
      else "Log-odds by region and migratory status",
      subtitle = "Points = medians; bars = 95% bootstrap intervals"
    ) +
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
    labs(
      x = x_lab,
      y = "Region (HemFly)",
      title = if (scale == "prob")
        "Migration effect on Pr(decline) by region"
      else
        "Migration effect on log-odds by region",
      subtitle = "Median and 95% bootstrap intervals"
    ) +
    theme_bw(base_size = 12)
}

get_boot_col <- function(fit, colname) {
  bm <- as.matrix(fit$bootstrap)
  
  # orient using names, not dimensions
  if (colname %in% colnames(bm)) {
    B <- bm
  } else if (colname %in% rownames(bm)) {
    B <- t(bm)
  } else {
    stop(sprintf("Column '%s' not found in fit$bootstrap.", colname))
  }
  
  as.numeric(B[, colname])
}

# ── Load + prepare data ─────────────────────────────────────────────────
data <- read.csv(
  data_file,
  header = TRUE,
  na.strings = "",
  stringsAsFactors = FALSE
)

data <- data %>%
  mutate(
    Hem3 = case_when(
      Hemisphere %in% c(1, 2) ~ "Breed_North_only",
      Hemisphere == 3         ~ "Present_South_only",
      Hemisphere %in% c(4, 5) ~ "Breed_Both",
      TRUE                    ~ NA_character_
    ),
    Hem3 = factor(Hem3, levels = c("Breed_North_only", "Present_South_only", "Breed_Both"))
  ) %>%
  filter(!is.na(Hem3))

data$Hem3n <- as.numeric(data$Hem3)

data <- data %>%
  mutate(
    flyway_combo = case_when(
      American == 1 & AfroPal == 0 & Asian == 0 ~ "Am_only",
      American == 0 & AfroPal == 1 & Asian == 0 ~ "Af_only",
      American == 0 & AfroPal == 0 & Asian == 1 ~ "As_only",
      American == 1 & AfroPal == 1 & Asian == 0 ~ "Am_Af",
      American == 1 & AfroPal == 0 & Asian == 1 ~ "Am_As",
      American == 0 & AfroPal == 1 & Asian == 1 ~ "Af_As",
      American == 1 & AfroPal == 1 & Asian == 1 ~ "Am_Af_As",
      TRUE                                       ~ NA_character_
    ),
    flyway_combo = factor(flyway_combo)
  ) %>%
  filter(
    flyway_combo != "Am_Af",
    !(Hem3n == 2 & flyway_combo %in% c("Am_As", "Af_As", "Am_Af_As"))
  ) %>%
  droplevels()

# keep your original column removal if still needed
data <- select(data, -c(19:25))

data$populationTrend <- as.factor(data$populationTrend)
data <- filter(data, populationTrend %in% c("Decreasing", "Increasing", "Stable"))
data$decline <- ifelse(data$populationTrend == "Decreasing", 1, 0)

allFactors <- c("decline", "migratory", "Hemisphere", "Hem3n",
                "American", "AfroPal", "Asian", "flyway_combo")
data[allFactors] <- lapply(data[allFactors], factor)

data$y <- as.integer(as.character(data$decline))

data$EOO <- as.numeric(data$EOO)
data <- data[!is.na(data$EOO), ]

data$EOO_log <- log(data$EOO + 1)
data$EOO_log_cent <- as.numeric(scale(data$EOO_log, center = TRUE, scale = FALSE))

data$rowID <- paste0("sp_", seq_len(nrow(data)))
rownames(data) <- data$rowID

data <- data %>%
  mutate(HemFly = interaction(Hem3n, flyway_combo, sep = "_", drop = TRUE))

desired_levels <- c(
  "1_Am_only", "1_Af_only", "1_As_only", "1_Af_As", "1_Am_As", "1_Am_Af_As",
  "2_Am_only", "2_Af_only", "2_As_only",
  "3_Am_only", "3_Af_only", "3_As_only", "3_Af_As", "3_Am_As", "3_Am_Af_As"
)

desired_levels <- desired_levels[desired_levels %in% unique(as.character(data$HemFly))]
data$HemFly <- factor(data$HemFly, levels = desired_levels)
data$HemFly <- relevel(data$HemFly, ref = "1_Am_only")

# Store cleaned data for reuse
saveRDS(data, file.path(out_dir, "Data_for_plots_phyloglm.rds"))

# ── Load + scale tree ───────────────────────────────────────────────────
tree <- readRDS(tree_file)
tree$edge.length <- tree$edge.length / max(nodeHeights(tree))

# keep only species present in both data and tree
keep_tips <- intersect(tree$tip.label, rownames(data))
tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))
data <- data[keep_tips, , drop = FALSE]

# ── Descriptive summary plot (like old Summary_S1 style) ───────────────
desc_df <- data %>%
  mutate(
    Migratory = ifelse(as.character(migratory) == levels(migratory)[2], "Migratory", "Non-migratory"),
    Decline = ifelse(as.character(decline) == "1", "Declining", "Not declining")
  ) %>%
  count(Migratory, Decline) %>%
  group_by(Migratory) %>%
  mutate(Prop = 100 * n / sum(n)) %>%
  ungroup()

p_desc <- ggplot(desc_df, aes(x = Migratory, y = Prop, fill = interaction(Migratory, Decline))) +
  geom_col() +
  scale_fill_manual(values = c("#f2e3c0", "#a1ddd8", "#E1BE6A", "#40B0A6")) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = NULL, y = NULL, title = "Observed proportions by migratory status") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# ── Fit models ──────────────────────────────────────────────────────────
# Basic GLMs
m1  <- glm(decline ~ migratory, family = binomial, data = data)
m1a <- glm(decline ~ migratory + EOO_log_cent, family = binomial, data = data)
m1b <- glm(decline ~ migratory * EOO_log_cent, family = binomial, data = data)

m2  <- glm(decline ~ HemFly:migratory - 1, family = binomial, data = data)
m2a <- glm(decline ~ HemFly:migratory + EOO_log_cent - 1, family = binomial, data = data)
m3  <- glm(decline ~ migratory * HemFly, family = binomial, data = data)
m4  <- glm(decline ~ migratory * HemFly + EOO_log_cent, family = binomial, data = data)

# Phylogenetic models with bootstrap retained for plotting
fit_phylo_mig <- phyloglm(
  decline ~ migratory,
  data = data,
  phy = tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4,
  boot = 1000,
  save = TRUE,
  start.beta = coef(m1),
  start.alpha = 0.5
)

fit_phyloHemFly <- phyloglm(
  decline ~ migratory * HemFly,
  data = data,
  phy = tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4,
  boot = 1000,
  save = TRUE,
  start.beta = coef(m3),
  start.alpha = 0.5
)

fit_phyloHemFlyEOO <- phyloglm(
  decline ~ migratory * HemFly + EOO_log_cent,
  data = data,
  phy = tree,
  method = "logistic_MPLE",
  btol = 50,
  log.alpha.bound = 4,
  boot = 1000,
  save = TRUE,
  start.beta = coef(m4),
  start.alpha = 0.5
)

saveRDS(fit_phylo_mig,       file.path(out_dir, "fit_phylo_mig.rds"))
saveRDS(fit_phyloHemFly,     file.path(out_dir, "fit_phyloHemFly.rds"))
saveRDS(fit_phyloHemFlyEOO,  file.path(out_dir, "fit_phyloHemFlyEOO.rds"))