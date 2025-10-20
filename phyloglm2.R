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
data <- read.csv("Data/IUCNdata_updated_Taxonomy_030925.csv", header=TRUE,
                 na.strings="", stringsAsFactors=FALSE)

data <- data %>%
  mutate(
    Hem3 = case_when(
      Hemisphere %in% c(1, 2) ~ "Breed_North_only",   # combine 1 & 2
      Hemisphere == 3         ~ "Present_South_only", # keep 3
      Hemisphere %in% c(4, 5) ~ "Breed_Both",         # combine 4 & 5
      TRUE ~ NA_character_
    ),
    Hem3 = factor(Hem3, levels = c("Breed_North_only", "Present_South_only", "Breed_Both"))
  ) %>%
  filter(!is.na(Hem3))

data$Hem3n <- as.numeric(data$Hem3)

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
         !(Hem3n == 2 & flyway_combo %in% c("Am_As","Af_As","Am_Af_As"))) %>%
  droplevels()
data <- select(data, -c(19:25))

data %>% 
  count(Hem3, flyway_combo) %>% 
  tidyr::pivot_wider(names_from = flyway_combo, values_from = n, values_fill = 0)

# Seems sensible to collapse --------------------------------------------------#

#data <- data %>% 
#  mutate(
#    flyway2 = fct_collapse(
#      flyway_combo,
#      # keep these three
#      Pure_American = "Am_only",
#      Pure_AfroPal   = "Af_only",
#      Pure_Asian     = "As_only",
#      # lump all the rest
#      Other = c("Am_Af", "Am_As", "Af_As", "Am_Af_As")
#    )
#  )

#data %>% 
#  count(Hemisphere, flyway2) %>% 
#  pivot_wider(names_from = flyway2, values_from = n, values_fill = 0)

#------------------------------------------------------------------------------#

## Sort decline
data$populationTrend <- as.factor(data$populationTrend)
summary(data$populationTrend)
data <- filter(data, populationTrend == "Decreasing" | populationTrend == "Increasing" | populationTrend == "Stable")
data$decline <- ifelse(data$populationTrend == "Decreasing", 1, 0)

# factor categorical predictors
allFactors <- c("decline","migratory","Hemisphere", "Hem3n",
                "American","AfroPal","Asian", "flyway_combo")#, "HemFly6", "HemFly4")

data[allFactors] <- lapply(data[allFactors], factor)

# Numeric response
data$y <- as.integer(as.character(data$decline))

# Log‐transform & center EOO
data$EOO       <- as.numeric(data$EOO)
data$EOO_log   <- log(data$EOO + 1)
data$EOO_log_cent  <- scale(data$EOO_log, center=TRUE, scale=FALSE)
remove <- which(is.na(data$EOO))
data <- data[-remove,]

data %>% 
  count(Hem3, flyway_combo) %>% 
  tidyr::pivot_wider(names_from = flyway_combo, values_from = n, values_fill = 0)

# Label rows
data$rowID <- paste0("sp_", seq_len(nrow(data)))
names(data)

# Read tree
#tree <- read.tree("Data/ult_5k_tree.trees")
tree <- readRDS("Data/NewFixedTree_030925.rds")
tree$edge.length <- tree$edge.length / max(nodeHeights(tree))

summary(tree)
#eps  <- 1e-8

#for(i in seq_len(nrow(data))) {
#  newName <- data$rowID[i]
#  anchor  <- data$animal[i]
#  tree <- bind.tip(tree,
#                   tip.label = newName,
#                   where     = which(tree$tip.label == anchor),
#                   position  = 0)
#}

# now jitter so no zero‐length edges remain
#tree$edge.length[tree$edge.length <= 0] <- eps

#keep_sp <- intersect(data$rowID, tree$tip.label)
#tree    <- drop.tip(tree, setdiff(tree$tip.label, keep_sp))
#saveRDS(tree, "Data/NewFixedTree_030925.rds")
rownames(data) <- data$rowID

# allow up to 2 GiB of globals
options(future.globals.maxSize = +Inf)

#------------------------------------------------------------------------------#
m1 <- glm(decline ~ migratory, family = binomial, data)
m1$coefficients
summary(m1)

m1a <- glm(decline ~ migratory + EOO_log_cent, family = binomial, data)
summary(m1a)

m1b <- glm(decline ~ migratory * EOO_log_cent, family = binomial, data)
summary(m1b)


#-----------------------------------------------------------------------------#
# Try with combined Hem/Flyway
data <- data %>%
  mutate(HemFly = interaction(Hem3n, flyway_combo, sep = "_", drop = TRUE))

# Check numbers...
with(data, table(migratory, HemFly))

# Change ref category
data$HemFly <- relevel(data$HemFly, ref = "1_Am_only")
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

m2 <- glm(decline ~ HemFly:migratory - 1, family = binomial, data)
m2$coefficients
summary(m2)

coefs <- coef(m2)
vc    <- vcov(m2)
se    <- sqrt(diag(vc))

wald <- tibble(
  Cell = names(coefs),
  logit = unname(coefs),
  se    = se
) %>%
  mutate(
    lo_logit = logit - 1.96*se,
    hi_logit = logit + 1.96*se,
    p_med    = plogis(logit),
    p_lo95   = plogis(lo_logit),
    p_hi95   = plogis(hi_logit)
  )

wald <- wald %>%
  tidyr::separate(Cell, into = c("HemFly","Migratory"), sep = ":", remove = FALSE) %>%
  mutate(
    HemFly = sub("^HemFly", "", HemFly),  # remove "HemFly" prefix
    HemFly = factor(HemFly, levels = levels(data$HemFly)), # keep ordering
    Migratory = sub("^migratory", "", Migratory),
    Migratory = ifelse(Migratory %in% c("0","1"),
                       ifelse(Migratory=="1","Migrant","Resident"),
                       Migratory),
    Migratory = factor(Migratory, levels = c("Resident", "Migrant"))
  )

ggplot(wald, aes(x = p_med, y = HemFly, colour = Migratory)) +
  geom_linerange(aes(xmin = p_lo95, xmax = p_hi95),
                 position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.6) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Pr(decline)", y = "Region (HemFly)",
       title = "GLM estimates with Wald 95% CIs") +
  theme_bw(base_size = 12)

##------------------------------------------------------------------------------
# decline ~ HemFly:migratory + EOO - 1

m2a <- glm(decline ~ HemFly:migratory + EOO_log_cent - 1, family = binomial, data)
m2a$coefficients
summary(m2a)

coefs <- coef(m2a)[-1]
vc    <- vcov(m2a)
se    <- sqrt(diag(vc))[-1]

wald <- tibble(
  Cell = names(coefs),
  logit = unname(coefs),
  se    = se
) %>%
  mutate(
    lo_logit = logit - 1.96*se,
    hi_logit = logit + 1.96*se,
    p_med    = plogis(logit),
    p_lo95   = plogis(lo_logit),
    p_hi95   = plogis(hi_logit)
  )

wald <- wald %>%
  tidyr::separate(Cell, into = c("HemFly","Migratory"), sep = ":", remove = FALSE) %>%
  mutate(
    HemFly = sub("^HemFly", "", HemFly),  # remove "HemFly" prefix
    HemFly = factor(HemFly, levels = levels(data$HemFly)), # keep ordering
    Migratory = sub("^migratory", "", Migratory),
    Migratory = ifelse(Migratory %in% c("0","1"),
                       ifelse(Migratory=="1","Migrant","Resident"),
                       Migratory),
    Migratory = factor(Migratory, levels = c("Resident", "Migrant"))
  )

ggplot(wald, aes(x = p_med, y = HemFly, colour = Migratory)) +
  geom_linerange(aes(xmin = p_lo95, xmax = p_hi95),
                 position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.6) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Pr(decline)", y = "Region (HemFly)",
       title = "GLM estimates with Wald 95% CIs",
       subtitle = "decline ~ HemFly:migratory + EOO - 1") +
  theme_bw(base_size = 12)

#-------------------------------------------------------------------------------
# decline ~ migratory * Hemfly with phylo signal
#Fit in parallel...
plan(multisession, workers = 50)

m3 <- glm(decline ~ migratory * HemFly, family = binomial, data)
m3$coefficients
startB <- m3$coefficients
startA <- 0.5

# Wrap the call to reduce what gets exported to workers
fit_phyloHemFly <- phyloglm(
  decline ~ migratory * HemFly,
  data = data, phy = tree,
  method = "logistic_MPLE",
  btol = 50, 
  log.alpha.bound = 4,
  boot = 1000,
  save = TRUE,
  #full.matrix = TRUE,
  start.beta = startB,
  start.alpha = startA
  )

saveRDS(fit_phyloHemFly , "Outputs/fit_phyloHemFly_1AmOnly_1000_scaledT.rds")

#-------------------------------------------------------------------------------
# decline ~ migratory * Hemfly + EOO with phylo signal
#Fit in parallel...
plan(multisession, workers = 50)

m4 <- glm(decline ~ migratory * HemFly + EOO_log_cent, family = binomial, data)
m4$coefficients
startB <- m4$coefficients
startA <- 0.5

# Wrap the call to reduce what gets exported to workers
fit_phyloHemFlyEOO <- phyloglm(
  decline ~ migratory * HemFly + EOO_log_cent,
  data = data, phy = tree,
  method = "logistic_MPLE",
  btol = 50, 
  log.alpha.bound = 4,
  boot = 1000,
  save = TRUE,
  #full.matrix = TRUE,
  start.beta = startB,
  start.alpha = startA
)

saveRDS(fit_phyloHemFlyEOO , "Outputs/fit_phyloHemFlyEOO_1AmOnly_1000_scaledT.rds")
summary(fit_phyloHemFlyEOO)
fit2 <- readRDS("Outputs/fit_phyloHemFlyEOO_1AmOnly_1000_scaledT.rds")
summary(fit2)

fit_phyloHemFlyEOO$coefficients
fit2$coefficients
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------#
# Loop over different factor levels...

# Get all possible reference levels
refs <- levels(data$HemFly)

for (ref in refs) {
  # Relevel HemFly
  data$HemFly <- relevel(data$HemFly, ref = ref)
  
  # Fit the standard GLM
  m1 <- glm(decline ~ migratory * HemFly + EOO_log_cent,
            family = binomial,
            data   = data)
  coefs <- m1$coefficients
  
  # Fit the phylogenetic GLM
  fit_phylo <- phyloglm(decline ~ migratory * HemFly + EOO_log_cent,
                        data            = data,
                        phy             = tree,
                        method          = "logistic_MPLE",
                        btol            = 50,
                        log.alpha.bound = 4,
                        boot            = 10)
  
  # Save model output
  saveRDS(fit_phylo,
          file = paste0("Outputs/FitPhyloHemFlyEOO_10_ref_", ref, ".rds"))
  
  # print progress
  message("Finished ref = ", ref)
}

#------------------------------------------------------------------------------
# relevel to Hem 2 - Am Only
data$HemFly <- relevel(data$HemFly, ref = "2_Am_only")

# Fit with EOO and bootstrap
m1 <- glm(decline ~ EOO_log_cent, family = binomial, data)
m1$coefficients
startB <- m1$coefficients
startA <- 0.5

fit_phyloHemFlyEOO <- phyloglm(
  decline ~ migratory * HemFly + EOO_log_cent,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 50,
  log.alpha.bound = 4,
  boot            = 0
)

summary(fit_phyloHemFlyEOO)
saveRDS(fit_phyloHemFlyEOO, "Outputs/FitPhyloHemFlyEOO_boot100.rds")

# relevel to Hem 2 - Am Only
data$HemFly <- relevel(data$HemFly, ref = "2_Am_only")

# Fit with EOO and bootstrap
m1 <- glm(decline ~ migratory * HemFly + EOO_log_cent, family = binomial, data)
m1$coefficients
startB <- m1$coefficients
startA <- 0.5

fit_phyloHemFlyEOO <- phyloglm(
  decline ~ migratory * HemFly + EOO_log_cent,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 50,
  log.alpha.bound = 4,
  boot            = 100
)

summary(fit_phyloHemFly)
saveRDS(fit_phyloHemFlyEOO, "Outputs/FitPhyloHemFlyEOO_boot100.rds")
