# ======================================================================
# SINGLE-HEMFLY DENSITY PLOTS
# - Resident and Migrant densities on the same axis
# - slightly overlapping
# - with median dot + 95% interval at the bottom of each density
# - and Migrant - Resident difference as a separate plot
# ======================================================================

rm(list = ls())

# ── Packages ────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

# ── Paths ───────────────────────────────────────────────────────────────
setwd("~/Migration")
out_dir <- "Outputs"

# ── Load saved objects ──────────────────────────────────────────────────
data <- readRDS(file.path(out_dir, "Data_for_plots_phyloglm.rds"))
fit_phyloHemFlyEOO <- readRDS(file.path(out_dir, "fit_phyloHemFlyEOO.rds"))

# ── Helpers ─────────────────────────────────────────────────────────────
invlogit <- function(x) 1 / (1 + exp(-x))

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
    stop("RHS must contain 'HemFly' for these plots.")
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

predict_draws_long <- function(fit, dat, eoo = 0, scale = c("prob", "logit")) {
  pred <- predict_draws_grid(fit, dat, eoo = eoo, scale = scale)
  D <- pred$draws
  grid <- pred$grid
  
  long_df <- as.data.frame(t(D))
  names(long_df) <- paste0("draw_", seq_len(ncol(long_df)))
  
  bind_cols(grid, long_df) %>%
    pivot_longer(
      cols = starts_with("draw_"),
      names_to = "draw",
      values_to = "value"
    ) %>%
    mutate(
      HemFly = factor(HemFly, levels = levels(dat$HemFly)),
      Mig = ifelse(as.character(migratory) == levels(dat$migratory)[2], "Migrant", "Resident"),
      Mig = factor(Mig, levels = c("Resident", "Migrant"))
    )
}

predict_mig_diff_long <- function(fit, dat, eoo = 0, scale = c("prob", "logit")) {
  pred <- predict_draws_grid(fit, dat, eoo = eoo, scale = scale)
  D <- pred$draws
  grid <- pred$grid
  
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
    rowid_to_column("draw") %>%
    pivot_longer(
      cols = -draw,
      names_to = "HemFly",
      values_to = "value"
    ) %>%
    mutate(HemFly = factor(HemFly, levels = levels(dat$HemFly)))
}

summarise_dist <- function(df) {
  df %>%
    summarise(
      lo95 = quantile(value, 0.025, na.rm = TRUE),
      med  = median(value, na.rm = TRUE),
      hi95 = quantile(value, 0.975, na.rm = TRUE)
    )
}

# ── Build long data once ────────────────────────────────────────────────
paired_long <- predict_draws_long(
  fit_phyloHemFlyEOO,
  data,
  eoo = 0,
  scale = "prob"
)

diff_long <- predict_mig_diff_long(
  fit_phyloHemFlyEOO,
  data,
  eoo = 0,
  scale = "prob"
)

paired_xlim <- c(0, 1)

max_abs_diff <- max(abs(diff_long$value), na.rm = TRUE)
pad_diff <- max_abs_diff * 0.08
diff_xlim <- c(-max_abs_diff - pad_diff, max_abs_diff + pad_diff)

paired_breaks <- seq(0, 1, by = 0.25)
diff_breaks <- pretty(diff_xlim, n = 5)
diff_breaks <- diff_breaks[diff_breaks >= diff_xlim[1] & diff_breaks <= diff_xlim[2]]

# ── One HemFly: resident + migrant slightly overlapping ────────────────
plot_hemfly_paired_density <- function(hf,
                                       paired_df = paired_long,
                                       xlim = paired_xlim,
                                       breaks = paired_breaks,
                                       resident_fill = "#F4B6AE",
                                       migrant_fill = "#73D2D6") {
  
  d_sub <- paired_df %>% filter(HemFly == hf)
  
  s_sub <- d_sub %>%
    group_by(Mig) %>%
    summarise(
      lo95 = quantile(value, 0.025, na.rm = TRUE),
      med  = median(value, na.rm = TRUE),
      hi95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_base = c(0.03, 0.14),   # closer together = overlap
      y_dot  = y_base,
      y_lab  = c(0.09, 0.20)
    )
  
  dens_res <- density(
    d_sub$value[d_sub$Mig == "Resident"],
    from = xlim[1], to = xlim[2], na.rm = TRUE
  )
  dens_mig <- density(
    d_sub$value[d_sub$Mig == "Migrant"],
    from = xlim[1], to = xlim[2], na.rm = TRUE
  )
  
  max_y <- max(c(dens_res$y, dens_mig$y))
  height_scale <- 0.22
  
  dens_df <- bind_rows(
    tibble(
      x = dens_res$x,
      density = dens_res$y / max_y * height_scale,
      Mig = "Resident",
      y0 = 0.05
    ),
    tibble(
      x = dens_mig$x,
      density = dens_mig$y / max_y * height_scale,
      Mig = "Migrant",
      y0 = 0.16
    )
  ) %>%
    mutate(y = y0 + density)
  
  ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(x = x, ymin = y0, ymax = y, fill = Mig),
      alpha = 0.78,
      colour = NA
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y, colour = Mig),
      linewidth = 0.45
    ) +
    geom_segment(
      data = s_sub,
      aes(x = lo95, xend = hi95, y = y_base, yend = y_base),
      linewidth = 0.6,
      colour = "black"
    ) +
    geom_point(
      data = s_sub,
      aes(x = med, y = y_dot),
      size = 1.8,
      colour = "black"
    ) +
    geom_text(
      data = s_sub,
      aes(x = xlim[1] + 0.02 * diff(xlim), y = y_lab, label = Mig),
      hjust = 0,
      size = 3
    ) +
    scale_fill_manual(values = c("Resident" = resident_fill, "Migrant" = migrant_fill)) +
    scale_colour_manual(values = c("Resident" = resident_fill, "Migrant" = migrant_fill)) +
    coord_cartesian(xlim = xlim, ylim = c(0, 0.42), expand = FALSE) +
    scale_x_continuous(
      breaks = breaks,
      labels = label_percent(accuracy = 1)
    ) +
    labs(
      title = as.character(hf),
      #subtitle = "Resident and Migrant",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.ticks.x = element_line(linewidth = 0.3),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ── One HemFly: difference density ─────────────────────────────────────
plot_hemfly_difference_density <- function(hf,
                                           diff_df = diff_long,
                                           xlim = diff_xlim,
                                           breaks = diff_breaks,
                                           fill = "grey72") {
  
  d_sub <- diff_df %>% filter(HemFly == hf)
  s_sub <- summarise_dist(d_sub)
  
  ggplot(d_sub, aes(x = value)) +
    geom_density(
      fill = fill,
      alpha = 0.9,
      colour = "grey20",
      linewidth = 0.5,
      adjust = 1
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_segment(
      data = s_sub,
      aes(x = lo95, xend = hi95, y = 0.02, yend = 0.02),
      inherit.aes = FALSE,
      linewidth = 0.6,
      colour = "black"
    ) +
    geom_point(
      data = s_sub,
      aes(x = med, y = 0.02),
      inherit.aes = FALSE,
      size = 1.8,
      colour = "black"
    ) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    scale_x_continuous(
      breaks = breaks,
      labels = label_percent(accuracy = 1)
    ) +
    labs(
      title = as.character(hf),
      subtitle = "Migrant − Resident",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.ticks.x = element_line(linewidth = 0.3),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ── Choose one HemFly to plot ───────────────────────────────────────────
levels(data$HemFly)
#levels(data$HemFly)
# "1_Am_only"  "1_Af_only"  "1_As_only"  "1_Af_As"    "1_Am_As"    "1_Am_Af_As" "2_Am_only"  "2_Af_only"  "2_As_only" 
# "3_Am_only"  "3_Af_only"  "3_As_only"  "3_Af_As"    "3_Am_As"    "3_Am_Af_As"

chosen_hemfly <- "1_Af_only"

# ── Make plots ──────────────────────────────────────────────────────────
p_paired <- plot_hemfly_paired_density(chosen_hemfly)
#p_diff   <- plot_hemfly_difference_density(chosen_hemfly)

print(p_paired)
#print(p_diff)

# ── Save if wanted ──────────────────────────────────────────────────────
ggsave(
  filename = file.path(out_dir, paste0("HemFly_", chosen_hemfly, "_paired_density_overlap.png")),
  plot = p_paired,
  width = 2.8,
  height = 2.0,
  dpi = 300,
  bg = "transparent"
)

ggsave(
  filename = file.path(out_dir, paste0("HemFly_", chosen_hemfly, "_difference_density.png")),
  plot = p_diff,
  width = 2.8,
  height = 1.9,
  dpi = 300,
  bg = "transparent"
)


#=========================================================================

hemfly_labels <- c(
  "1_Am_only"   = "American Flyway only",
  "1_Af_only"   = "African Flyway only",
  "1_As_only"   = "Asian Flyway only",
  "1_Af_As"     = "African + Asian Flyways",
  "1_Am_As"     = "American + Asian Flyways",
  "1_Am_Af_As"  = "American + African + Asian Flyways",
  "2_Am_only"   = "American Flyway only",
  "2_Af_only"   = "African Flyway only",
  "2_As_only"   = "Asian Flyway only",
  "3_Am_only"   = "American Flyway only",
  "3_Af_only"   = "African Flyway only",
  "3_As_only"   = "Asian Flyway only",
  "3_Af_As"     = "African + Asian Flyways",
  "3_Am_As"     = "American + Asian Flyways",
  "3_Am_Af_As"  = "American + African + Asian Flyways"
)

pretty_hemfly_label <- function(hf) {
  out <- hemfly_labels[[hf]]
  if (is.null(out)) out <- hf
  out
}

safe_file_label <- function(hf) {
  hemi <- substr(hf, 1, 1)
  hemi_lab <- c("1" = "North", "2" = "South", "3" = "Both")[[hemi]]
  fly_lab <- pretty_hemfly_label(hf)
  
  fly_lab <- gsub("\\+", "plus", fly_lab)
  fly_lab <- gsub("[[:space:]]+", "_", fly_lab)
  fly_lab <- gsub("[^A-Za-z0-9_]", "", fly_lab)
  
  paste0(hemi_lab, "_", fly_lab)
}

plot_hemfly_paired_density <- function(hf,
                                       paired_df = paired_long,
                                       xlim = paired_xlim,
                                       breaks = paired_breaks,
                                       resident_fill = "#F4B6AE",
                                       migrant_fill = "#73D2D6") {
  
  d_sub <- paired_df %>% filter(HemFly == hf)
  
  if (nrow(d_sub) == 0) {
    stop(paste0(
      "No rows found for HemFly = '", hf, "'. Valid levels are:\n",
      paste(levels(paired_df$HemFly), collapse = ", ")
    ))
  }
  
  s_sub <- d_sub %>%
    group_by(Mig) %>%
    summarise(
      lo95 = quantile(value, 0.025, na.rm = TRUE),
      med  = median(value, na.rm = TRUE),
      hi95 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_base = ifelse(Mig == "Resident", 0.03, 0.14),
      y_dot  = y_base,
      y_lab  = ifelse(Mig == "Resident", 0.09, 0.20)
    )
  
  dens_list <- list()
  
  if (any(d_sub$Mig == "Resident")) {
    dens_res <- density(
      d_sub$value[d_sub$Mig == "Resident"],
      from = xlim[1], to = xlim[2], na.rm = TRUE
    )
    dens_list[["Resident"]] <- tibble(
      x = dens_res$x,
      density_raw = dens_res$y,
      Mig = "Resident",
      y0 = 0.05
    )
  }
  
  if (any(d_sub$Mig == "Migrant")) {
    dens_mig <- density(
      d_sub$value[d_sub$Mig == "Migrant"],
      from = xlim[1], to = xlim[2], na.rm = TRUE
    )
    dens_list[["Migrant"]] <- tibble(
      x = dens_mig$x,
      density_raw = dens_mig$y,
      Mig = "Migrant",
      y0 = 0.16
    )
  }
  
  dens_df <- bind_rows(dens_list)
  
  max_y <- max(dens_df$density_raw)
  height_scale <- 0.22
  
  dens_df <- dens_df %>%
    mutate(
      density = density_raw / max_y * height_scale,
      y = y0 + density
    )
  
  ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(x = x, ymin = y0, ymax = y, fill = Mig),
      alpha = 0.78,
      colour = NA
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y, colour = Mig),
      linewidth = 0.45
    ) +
    geom_segment(
      data = s_sub,
      aes(x = lo95, xend = hi95, y = y_base, yend = y_base),
      linewidth = 0.6,
      colour = "black"
    ) +
    geom_point(
      data = s_sub,
      aes(x = med, y = y_dot),
      size = 1.8,
      colour = "black"
    ) +
    geom_text(
      data = s_sub,
      aes(x = xlim[1] + 0.02 * diff(xlim), y = y_lab, label = Mig),
      hjust = 0,
      size = 3
    ) +
    scale_fill_manual(values = c("Resident" = resident_fill, "Migrant" = migrant_fill)) +
    scale_colour_manual(values = c("Resident" = resident_fill, "Migrant" = migrant_fill)) +
    coord_cartesian(xlim = xlim, ylim = c(0, 0.42), expand = TRUE) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      title = pretty_hemfly_label(hf),
      x = NULL,
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.ticks.x = element_line(linewidth = 0.3),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
      #plot.margin = margin(3, 3, 3, 3)
    )
}

plot_hemfly_difference_density <- function(hf,
                                           diff_df = diff_long,
                                           xlim = diff_xlim,
                                           breaks = diff_breaks,
                                           fill = "grey72") {
  
  d_sub <- diff_df %>% filter(HemFly == hf)
  
  if (nrow(d_sub) == 0) {
    stop(paste0(
      "No rows found for HemFly = '", hf, "'. Valid levels are:\n",
      paste(levels(diff_df$HemFly), collapse = ", ")
    ))
  }
  
  s_sub <- summarise_dist(d_sub)
  
  ggplot(d_sub, aes(x = value)) +
    geom_density(
      fill = fill,
      alpha = 0.9,
      colour = "grey20",
      linewidth = 0.5,
      adjust = 1
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_segment(
      data = s_sub,
      aes(x = lo95, xend = hi95, y = 0.02, yend = 0.02),
      inherit.aes = FALSE,
      linewidth = 0.6,
      colour = "black"
    ) +
    geom_point(
      data = s_sub,
      aes(x = med, y = 0.02),
      inherit.aes = FALSE,
      size = 1.8,
      colour = "black"
    ) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    scale_x_continuous(
      breaks = breaks,
      labels = label_percent(accuracy = 1)
    ) +
    labs(
      title = pretty_hemfly_label(hf),
      subtitle = "Migrant − Resident",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.ticks.x = element_line(linewidth = 0.3),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(3, 3, 3, 3)
    )
}

#levels(data$HemFly)
# "1_Am_only"  "1_Af_only"  "1_As_only"  "1_Af_As"    "1_Am_As"    "1_Am_Af_As" "2_Am_only"  "2_Af_only"  "2_As_only" 
# "3_Am_only"  "3_Af_only"  "3_As_only"  "3_Af_As"    "3_Am_As"    "3_Am_Af_As"

chosen_hemfly <- "3_Am_Af_As"

p_paired <- plot_hemfly_paired_density(chosen_hemfly)
#p_diff   <- plot_hemfly_difference_density(chosen_hemfly)

print(p_paired)
#print(p_diff)

ggsave(
  filename = file.path(out_dir, paste0("HemFly_", safe_file_label(chosen_hemfly), "_paired_density_overlap.png")),
  plot = p_paired,
  width = 2.8,
  height = 2.0,
  dpi = 300,
  bg = "transparent"
)

ggsave(
  filename = file.path(out_dir, paste0("HemFly_", safe_file_label(chosen_hemfly), "_difference_density.png")),
  plot = p_diff,
  width = 2.8,
  height = 1.9,
  dpi = 300,
  bg = "transparent"
)
