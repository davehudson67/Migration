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
library(future)

setwd("~/Migration")
data <- read.csv("Data/IUCNdata_updated_July2025.csv", header=TRUE,
                 na.strings="", stringsAsFactors=FALSE)

# Recode flyways
data <- data %>% 
  mutate(
    flyway_combo = case_when(
      American==1 & AfroPal==0 & Asian==0 ~ "Am_only",
      American==0 & AfroPal==1 & Asian==0 ~ "Af_only",
      American==0 & AfroPal==0 & Asian==1 ~ "As_only",
      American==1 & AfroPal==1 & Asian==0 ~ "Am_Af",
      American==1 & AfroPal==0 & Asian==1 ~ "Am_As",
      American==0 & AfroPal==1 & Asian==1 ~ "Af_As",
      American==1 & AfroPal==1 & Asian==1 ~ "Am_Af_As",
      TRUE                                ~ NA_character_
    ),
    flyway_combo = factor(flyway_combo)
  ) %>%
  filter(flyway_combo != "Am_Af",
         !(Hemisphere == 2 & flyway_combo %in% c("Am_As","Af_As","Am_Af_As"))) %>%
  droplevels()
data <- select(data, -c(9:16, 20:27))

data %>% 
  count(Hemisphere, flyway_combo) %>% 
  tidyr::pivot_wider(names_from = flyway_combo, values_from = n, values_fill = 0)

# factor categorical predictors
allFactors <- c("decline","migratory","Hemisphere","Flyway",
                "American","AfroPal","Asian", "flyway_combo")#, "HemFly6", "HemFly4")

data[allFactors] <- lapply(data[allFactors], factor)

# Numeric response
data$y <- as.integer(as.character(data$decline))

# Add HemFLy
data <- data %>%
  mutate(HemFly = interaction(Hemisphere, flyway_combo, sep = "_", drop = TRUE))
# Change ref category
data$HemFly <- relevel(data$HemFly, ref = "2_Am_only")

# Log‐transform & center EOO
data$EOO       <- as.numeric(data$EOO)
data$EOO_log   <- log(data$EOO + 1)
data$EOO_log_cent  <- scale(data$EOO_log, center=TRUE, scale=FALSE)
remove <- which(is.na(data$EOO))
data <- data[-remove,]

# Label rows
data$rowID <- paste0("sp_", seq_len(nrow(data)))
names(data)

# Read tree
#tree <- read.tree("Data/ult_5k_tree.trees")
tree <- readRDS("Data/NewFixedTree.rds")

rownames(data) <- data$rowID

#------------------------------------------------------------------------------
# Load phylo object
fit_phylo <- readRDS("Outputs/FitPhyloHemFlyEOO_H2AmOnly_1000.rds")

#------------------------------------------------------------------------------
#
# Plot coefficicents....

# Bootstrap matrix (rows = bootstraps)
bmat <- fit_phylo$bootstrap

# Fixed-effect names and draws (drop alpha/logalpha if present)
beta_names <- names(coef(fit_phylo))
keep_cols  <- intersect(colnames(bmat), beta_names)
betas      <- as.data.frame(bmat[, keep_cols, drop = FALSE])
betas$.boot <- seq_len(nrow(betas))

# alpha as a scalar per bootstrap
a_col  <- grep("alpha", colnames(bmat), value = TRUE)[1]
alpha  <- if (!is.na(a_col)) bmat[, a_col] else rep(NA_real_, nrow(bmat))

# Long format for ggplot
post_long <- betas %>%
  pivot_longer(-.boot, names_to = "term", values_to = "value")

#------------------------------------------------------------------------------
# Faceted 'density' plots with glm lines

# Narrow, unimodal densities → stable estimates.
# Wide or multi-modal densities → instability; different resamples pull the coefficient to different regions.
# Blue far from the density’s center → MPLE point estimate is atypical relative to the bootstrap draws.
# Big gap vs GLM (red) → phylogeny meaningfully changes that coefficient (phylogenetic signal affecting that contrast).

glm_fit <- glm(decline ~ migratory * HemFly + EOO_log_cent, family = binomial, data)

gb <- coef(glm_fit)
glm_lines <- tibble(term = intersect(names(gb), unique(post_long$term)),
         glm  = unname(gb[term]))

phylo_lines <- tibble(
  term = names(coef(fit_phylo)[keep_cols]),
  phy  = unname(coef(fit_phylo)[keep_cols])
)

ggplot(post_long, aes(x = value)) +
  geom_density(fill = "grey90") +
  geom_vline(data = phylo_lines, aes(xintercept = phy), colour = "steelblue", linewidth = 0.6) +
  geom_vline(data = glm_lines, aes(xintercept = glm), colour = "firebrick", linetype = 2, linewidth = 0.6)+
  facet_wrap(~ term, scales = "free", ncol = 4) +
  labs(title = "Bootstrap coefficient distributions",
       subtitle = "Blue = phyloglm point estimate; Red dashed = GLM (if provided)",
       x = "Coefficient (log-odds)", y = "Density") +
  theme_minimal(base_size = 12)

#------------------------------------------------------------------------------
# Median and 95% CI ordered by SD with GLM...

# Width = uncertainty/instability under resampling.
# CI crossing 0 = not clearly different from 0 on the log-odds scale.
# Large shift between red × and the median → phylogenetic model materially differs from GLM.
# Very wide intervals for interactions often reflect sparse cells/near-separation.

summ <- post_long %>%
  group_by(term) %>%
  summarise(med = median(value),
            lo  = quantile(value, 0.025),
            hi  = quantile(value, 0.975),
            sd  = sd(value),
            .groups = "drop") %>%
  arrange(desc(sd))  # show the wobbliest first

# Add GLM estimates
summ <- left_join(summ, glm_lines, by = "term")

ggplot(summ, aes(x = med, y = reorder(term, sd))) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(size = 2) +
  geom_point(aes(x = glm), shape = 4, size = 2, color = "firebrick") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  labs(title = "Bootstrap 95% intervals (ordered by SD)",
       subtitle = "Point = median; × = GLM (if provided)",
       x = "Coefficient (log-odds)", y = NULL) +
  theme_minimal(base_size = 12)

#------------------------------------------------------------------------------
# Index plots (coefficient value by bootstrap index)
#
# Each point is a single bootstrap draw for that coefficient. 
# The x-axis is the bootstrap replicate index; points are colored red when α is very low (e.g., bottom 5%).

# Random scatter: good—no special subset of bootstraps drives the coefficient.
# Red points forming a band far from the rest: when α is low (weak phylogenetic signal), 
# this coefficient systematically shifts (often changes sign or becomes extreme).
library(ggforce)
post_idx <- post_long %>%
  mutate(alpha = alpha[.boot],
         low_alpha = ifelse(is.finite(alpha) & alpha <= quantile(alpha, 0.005, na.rm = TRUE),
                            "low α (bottom 0.5%)", "other"))

nrow <- 1; ncol <- 4
p_base <- ggplot(post_idx, aes(x = .boot, y = value)) +
  geom_point(aes(color = low_alpha), alpha = 0.5, size = 0.7) +
  scale_color_manual(values = c("low α (bottom 0.5%)" = "firebrick", "other" = "grey50")) +
  labs(x = "Bootstrap index", y = "Coefficient") +
  theme_minimal(base_size = 12)

npages <- ggforce::n_pages(p_base + facet_wrap_paginate(~ term, nrow = nrow, ncol = ncol))

# draw page 1
p_base +
  facet_wrap_paginate(~ term, nrow = nrow, ncol = ncol, page = 8, scales = "free") +
  labs(title = "Coefficient draws across bootstraps",
       subtitle = paste("Page 1 of", npages))




#------------------------------------------------------------------------------
# Alpha: beta relationship...

# For each coefficient, a scatter of α (x-axis) against the bootstrap 
# coefficient value (y-axis), with a smooth loess curve.

# Flat loess: coefficient is robust to the strength of phylogenetic signal.
# Strong positive/negative trend: coefficient depends on α; e.g., as α→0 
# (less phylogenetic structure), the slope moves to a different regime.

ggplot(post_long %>% mutate(alpha = alpha[.boot]),
         aes(x = alpha, y = value)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
    facet_wrap(~ term, scales = "free", ncol = 4) +
    labs(title = "Association between α and coefficients across bootstraps",
         x = "alpha", y = "Coefficient") +
    theme_minimal(base_size = 12)

#------------------------------------------------------------------------------
# define weird by alpha and/or by coefficient outliers
weird_alpha <- which(alpha <= quantile(alpha, 0.05, na.rm = TRUE))

weird_coef <- post_long %>%
  group_by(term) %>%
  mutate(lo = quantile(value, 0.005), hi = quantile(value, 0.995)) %>%
  ungroup() %>%
  filter(value < lo | value > hi) %>%
  distinct(.boot) %>%
  pull(.boot)
# vector of bootstrap indices flagged by coefficient outliers (for any term)

# Union of the two criteria (low-α OR outlier coefficients).
weird_boots <- sort(unique(c(weird_alpha, weird_coef)))
length(weird_boots); head(weird_boots)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Extract alpha draws from the bootstrap matrix
bmat <- fit_phylo$bootstrap
alpha_col <- grep("alpha", colnames(bmat), value = TRUE)[1]
stopifnot(!is.na(alpha_col))
alpha <- bmat[, alpha_col]

# 0.2 species-by-bootstrap inclusion indicator (YOU may need to point to the right slot)
# Assume inc[b, s] = 1 if the species s is INCLUDED in bootstrap b (flip if it's "excluded").
inc <- fit_phylo$bootdata
stopifnot(ncol(inc) == nrow(bmat))
sp_names <- rownames(inc)

# 3) Define "weird" bootstraps (e.g., lowest 5% alpha);
#    you can OR in slope criteria if you wish.
lowCut <- quantile(alpha, 0.05)
weird  <- alpha <= lowCut                     # logical length K (bootstraps)

# 4) Species-level influence scores (WORKS WITH SPECIES x BOOT inc)
# Inclusion rates in weird vs ok boots
rate_weird <- rowMeans(inc[, weird,  drop = FALSE])   # P(included | weird)
rate_ok    <- rowMeans(inc[, !weird, drop = FALSE])   # P(included | ok)
delta_rate <- rate_weird - rate_ok                    # negative => often excluded in weird boots

# Correlation of inclusion with alpha (across boots)
cor_alpha  <- apply(inc, 1, function(v) cor(v, alpha, use = "complete.obs"))

# Enrichment test: are EXCLUSIONS enriched among weird boots?
exc <- 1 - inc  # 1 = excluded
p_excl <- sapply(seq_len(nrow(exc)), function(i){
  a <- sum(exc[i, weird])     # excluded & weird
  b <- sum(1 - exc[i, weird]) # included & weird
  c <- sum(exc[i, !weird])    # excluded & ok
  d <- sum(1 - exc[i, !weird])# included & ok
  fisher.test(matrix(c(a,b,c,d), 2, 2, byrow = TRUE))$p.value
})
p_excl <- p.adjust(p_excl, method = "BH")

infl_tbl <- data.frame(
  species       = sp_names,
  delta_rate    = delta_rate,     # big negative -> suspect (absence linked to low alpha)
  cor_alpha     = cor_alpha,      # positive -> presence raises alpha
  p_excl_enrich = p_excl,
  stringsAsFactors = FALSE
)
infl_tbl <- infl_tbl[order(infl_tbl$delta_rate), ]
head(infl_tbl, 20)      # most suspicious (often excluded in weird boots)
tail(infl_tbl, 20)      # often included in weird boots (less suspicious)

ggplot(infl_tbl, aes(x = delta_rate, y = reorder(species, delta_rate))) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  labs(x = "Δ inclusion rate  [weird − ok]   (negative ⇒ excluded in weird boots)",
       y = NULL, title = "Species linked to low-α bootstraps") +
  theme_minimal(base_size = 12)

#-------------------------------------------------------------------------------
#===============================================================================

bmat <- fit_phylo$bootstrap              # K × parameters (rows = boots)
inc  <- fit_phylo$bootdata          # S × K (rows = species, cols = boots), 1 = included
stopifnot(ncol(inc) == nrow(bmat))
alpha_col <- grep("alpha", colnames(bmat), value = TRUE)[1]
alpha <- if (length(alpha_col)) bmat[, alpha_col] else rep(NA_real_, nrow(bmat))
if (length(alpha_col) && grepl("^log", alpha_col)) alpha <- exp(alpha)  # back-transform

# (Optional) grab specific slopes if you want to also flag slope-weirdness
coef_names <- names(coef(fit_phylo))
keep <- intersect(colnames(bmat), coef_names)
b_betas <- bmat[, keep, drop = FALSE]
# e.g., tweak these names to your model terms:
eoo_draw <- b_betas[, "EOO_log_cent"]
mig_draw <- b_betas[, "migratory1"]

lowA <- is.finite(alpha) & alpha <= quantile(alpha, 0.05, na.rm = TRUE)
weird_coef <- (eoo_draw <= quantile(eoo_draw, 0.005) | eoo_draw >= quantile(eoo_draw, 0.995) |
                 mig_draw <= quantile(mig_draw, 0.005) | mig_draw >= quantile(mig_draw, 0.995))
weird <- lowA | weird_coef
table(weird)

#------------------------------------
W  <- which(weird)
OK <- which(!weird)
S  <- nrow(inc)

ex_weird <- rowSums(1 - inc[, W,  drop = FALSE])  # exclusions among weird boots
ex_ok    <- rowSums(1 - inc[, OK, drop = FALSE])  # exclusions among normal boots
nW <- length(W); nO <- length(OK)

rate_exc_weird <- ex_weird / pmax(nW, 1)
rate_exc_ok    <- ex_ok    / pmax(nO, 1)
delta_exc      <- rate_exc_weird - rate_exc_ok     # > 0 ⇒ excluded more often in weird

# odds ratio of exclusion (add 0.5 Haldane prior to avoid zeros)
OR <- ((ex_weird + .5) / (nW - ex_weird + .5)) / ((ex_ok + .5) / (nO - ex_ok + .5))

# per-species Fisher tests
p_excl <- vapply(seq_len(S), function(i){
  a <- ex_weird[i]; b <- nW - ex_weird[i]
  c <- ex_ok[i];    d <- nO - ex_ok[i]
  fisher.test(matrix(c(a,b,c,d), 2, 2, byrow = TRUE))$p.value
}, numeric(1))
padj <- p.adjust(p_excl, "BH")

infl_tbl <- tibble(
  species          = rownames(inc),
  rate_exc_weird,
  rate_exc_ok,
  delta_exc,
  OR_exclusion = OR,
  p_excl_enrich = padj
) %>%
  arrange(desc(delta_exc), desc(OR_exclusion))

head(infl_tbl, 15)

topN <- 25
plot_tbl <- infl_tbl %>% slice_head(n = topN) %>%
  tidyr::pivot_longer(c(rate_exc_weird, rate_exc_ok),
                      names_to = "set", values_to = "rate")

ggplot(plot_tbl,
       aes(y = reorder(species, infl_tbl$delta_exc[match(species, infl_tbl$species)]),
           x = rate, fill = set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("rate_exc_weird" = "firebrick", "rate_exc_ok" = "grey60"),
                    labels = c("weird boots", "other boots")) +
  labs(x = "Exclusion rate", y = NULL,
       title = "Species excluded more often in 'weird' bootstraps",
       fill = NULL) +
  theme_minimal(base_size = 12)

#===================================================================================

## Included (drivers of wierdness)

# Inputs
bmat <- fit_phylo$bootstrap                 # K x params
inc  <- fit_phylo$bootdata            # S x K, 1 = included
stopifnot(ncol(inc) == nrow(bmat))

# alpha per bootstrap (back-transform if stored on log-scale)
alpha_col <- grep("alpha", colnames(bmat), value = TRUE)[1]
alpha <- if (length(alpha_col)) bmat[, alpha_col] else rep(NA_real_, nrow(bmat))
if (length(alpha_col) && grepl("^log", alpha_col)) alpha <- exp(alpha)

# Define "weird" boots (e.g., lowest 5% alpha). You can OR in slope outliers too.
weird <- is.finite(alpha) & alpha <= quantile(alpha, 0.05, na.rm = TRUE)

W  <- which(weird)
OK <- which(!weird)
nW <- length(W); nO <- length(OK)

# --- Inclusion diagnostics ---
# Inclusion rates when weird vs ok
inc_weird <- rowSums(inc[, W,  drop = FALSE])          # counts of inclusion among weird boots
inc_ok    <- rowSums(inc[, OK, drop = FALSE])          # counts of inclusion among normal boots

rate_inc_weird <- inc_weird / pmax(nW, 1)
rate_inc_ok    <- inc_ok    / pmax(nO, 1)
delta_inc      <- rate_inc_weird - rate_inc_ok         # > 0 ⇒ included more often in weird boots

# Odds ratio for inclusion in weird vs ok (Haldane +0.5 to avoid zeros)
OR_inclusion <- ((inc_weird + .5) / (nW - inc_weird + .5)) / ((inc_ok + .5) / (nO - inc_ok + .5))

# Fisher test: is inclusion enriched among weird boots?
p_incl <- vapply(seq_len(nrow(inc)), function(i){
  a <- inc_weird[i]; b <- nW - inc_weird[i]   # included/excluded within weird
  c <- inc_ok[i];    d <- nO - inc_ok[i]      # included/excluded within ok
  fisher.test(matrix(c(a,b,c,d), 2, 2, byrow = TRUE))$p.value
}, numeric(1))
p_incl <- p.adjust(p_incl, method = "BH")

# Correlation with weirdness indicator (1 = weird)
is_weird <- rep(0L, ncol(inc)); is_weird[W] <- 1L
cor_weird <- apply(inc, 1, function(v) suppressWarnings(cor(v, is_weird, use = "complete.obs")))

incl_tbl <- data.frame(
  species        = rownames(inc),
  rate_inc_weird = rate_inc_weird,
  rate_inc_ok    = rate_inc_ok,
  delta_inc      = delta_inc,        # primary ranking: higher ⇒ present more often in weird runs
  OR_inclusion   = OR_inclusion,
  cor_weird      = cor_weird,        # positive ⇒ presence aligns with weird boots
  p_incl_enrich  = p_incl,
  stringsAsFactors = FALSE
)
incl_tbl <- incl_tbl[order(-incl_tbl$delta_inc, -incl_tbl$OR_inclusion), ]
head(incl_tbl, 15)

topN <- 25
plot_tbl <- incl_tbl %>%
  slice_head(n = topN) %>%
  pivot_longer(c(rate_inc_weird, rate_inc_ok),
               names_to = "set", values_to = "rate")

ggplot(plot_tbl,
       aes(y = reorder(species, incl_tbl$delta_inc[match(species, incl_tbl$species)]),
           x = rate, fill = set)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.6) +
  scale_fill_manual(values = c("rate_inc_weird" = "firebrick", "rate_inc_ok" = "grey65"),
                    labels = c("weird boots", "other boots")) +
  labs(x = "Inclusion rate", y = NULL,
       title = "Species included more often in 'weird' (low-α) bootstraps",
       fill = NULL) +
  theme_minimal(base_size = 12)





