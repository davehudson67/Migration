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

# factor categorical predictors
allFactors <- c("decline","migratory","Hemisphere","Flyway",
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

# Label rows
data$rowID <- paste0("sp_", seq_len(nrow(data)))
names(data)

# Read tree
#tree <- read.tree("Data/ult_5k_tree.trees")
tree <- readRDS("Data/NewFixedTree.rds")
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
#saveRDS(tree, "Data/NewFixedTree.rds")
rownames(data) <- data$rowID

# allow up to 2 GiB of globals
options(future.globals.maxSize = +Inf)

#------------------------------------------------------------------------------#

m1 <- glm(decline ~ migratory * Hemisphere * flyway_combo + EOO_log_cent, family = binomial, data)
m1$coefficients
m1$coefficients[is.na(m1$coefficients)] <- 1
startB <- m1$coefficients
startA <- 0.5

# Fit phylo logistic regression
fit_phylo <- phyloglm(
  decline ~ migratory * Hemisphere * flyway_combo + EOO_log_cent,
  data   = data,
  phy    = tree,
  method = "logistic_MPLE",
  btol   = 40,
  log.alpha.bound = 4,
  start.beta = startB,
  start.alpha = startA,
  boot = 1
)

summary(fit_phylo)
#-----------------------------------------------------------------------------#
# Drop the 3 way interaction term

#fit_phylo2 <- phyloglm(
#  decline ~ migratory * Hemisphere + migratory * flyway2 + Hemisphere * flyway2 + EOO_log_cent,
#  data            = data,
#  phy             = tree,
#  method          = "logistic_MPLE",
#  btol            = 50,
#  log.alpha.bound = 4,
#  boot            = 100
#)
#summary(fit_phylo2)
#saveRDS(fit_phylo2, "FitPhylo2.rds")

#-----------------------------------------------------------------------------#
# Try with combined Hem/Flyway

data <- data %>%
  mutate(HemFly = interaction(Hemisphere, flyway_combo, sep = "_", drop = TRUE))

# Check numbers...
with(data, table(migratory, HemFly))

# Change ref category
data$HemFly <- relevel(data$HemFly, ref = "2_Am_only")

# Fit with no EOO and no bootstrap
m1 <- glm(decline ~ migratory * HemFly, family = binomial, data)
m1$coefficients
startB <- m1$coefficients
startA <- 0.5

fit_phyloHemFly <- phyloglm(
  decline ~ migratory * HemFly,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 50,
  log.alpha.bound = 4,
  boot            = 5
)

summary(fit_phyloHemFly)
saveRDS(fit_phyloHemFly, "Outputs/FitPhyloHemFly.rds")
#------------------------------------------------------------------------------#
#Fit in parallel...
library(future)

plan(multisession, workers = 50)
#options(future.globals.maxSize = 8 * 1024^3)

# Wrap the call to reduce what gets exported to workers
run_phylo <- function(dat, tr, K){
  phylolm::phyloglm(
    decline ~ migratory * HemFly + EOO_log_cent,
    data = dat, phy = tr,
    method = "logistic_MPLE",
    btol = 50, log.alpha.bound = 4,
    boot = K,
    save = TRUE,
    full.matrix = TRUE
  )
}

fit_phyloHemFly <- run_phylo(data, tree, 1000)
saveRDS(fit_phyloHemFly , "Outputs/fit_phyloHemFlyEOO_2AmOnly_1000.rds")

#------------------------------------------------------------------------------#
# Fit with EOO
m1 <- glm(decline ~ migratory * HemFly + EOO_log_cent, family = binomial, data)
m1$coefficients
startB <- m1$coefficients
startA <- 0.5

plan(multisession, workers = 20)
plan()

fit_phyloHemFlyEOO <- phyloglm(
  decline ~ migratory * HemFly + EOO_log_cent,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 20,
  log.alpha.bound = 4,
  boot            = 5000,
  start.beta = startB,
  start.alpha = startA,
)

summary(fit_phyloHemFlyEOO)
saveRDS(fit_phyloHemFlyEOO, "Outputs/FitPhyloHemFlyEOO.rds")

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
