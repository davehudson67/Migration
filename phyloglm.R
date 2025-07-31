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

setwd("~/FraserCode")
data <- read.csv("Data/IUCN_data_B.csv", header=TRUE,
                 na.strings="", stringsAsFactors=FALSE)

data$rowID <- paste0("sp_", seq_len(nrow(data)))

names(data)

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
  )

data %>% 
  count(Hemisphere, flyway_combo) %>% 
  tidyr::pivot_wider(names_from = flyway_combo, values_from = n, values_fill = 0)

# Seems sensible to collapse --------------------------------------------------#

data <- data %>% 
  mutate(
    flyway2 = fct_collapse(
      flyway_combo,
      # keep these three
      Pure_American = "Am_only",
      Pure_AfroPal   = "Af_only",
      Pure_Asian     = "As_only",
      # lump all the rest
      Other = c("Am_Af", "Am_As", "Af_As", "Am_Af_As")
    )
  )

data %>% 
  count(Hemisphere, flyway2) %>% 
  pivot_wider(names_from = flyway2, values_from = n, values_fill = 0)

#------------------------------------------------------------------------------#

# factor categorical predictors
allFactors <- c("decline","migratory","Hemisphere","Artificial","Forest",
                "Grassland","Savanna","Shrubland","Wetlands",
                "American","AfroPal","Asian", "flyway2")#, "HemFly6", "HemFly4")

data[allFactors] <- lapply(data[allFactors], factor)

# Numeric response
data$y <- as.integer(as.character(data$decline))

# Log‐transform & center EOO
data$EOO       <- as.numeric(data$EOO)
data$EOO_log   <- log(data$EOO + 1)
data$EOO_log_cent  <- scale(data$EOO_log, center=TRUE, scale=FALSE)

# Read tree
#tree <- read.tree("Data/ult_5k_tree.trees")
tree <- readRDS("FixedTree.rds")
eps  <- 1e-8

#for(i in seq_len(nrow(data))) {
#  newName <- data$rowID[i]
#  anchor  <- data$animal[i]
#  tree <- bind.tip(tree,
#                   tip.label = newName,
#                   where     = which(tree$tip.label == anchor),
#                   position  = 0)
#}

# now jitter so no zero‐length edges remain
tree$edge.length[tree$edge.length <= 0] <- eps

keep_sp <- intersect(data$rowID, tree$tip.label)
tree    <- drop.tip(tree, setdiff(tree$tip.label, keep_sp))
rownames(data) <- data$rowID

# allow up to 2 GiB of globals
options(future.globals.maxSize = 2 * 1024^3)

m1 <- glm(decline ~ migratory * Hemisphere * flyway2, family = binomial, data)
m1$coefficients[23] <- 1
m1$coefficients[25] <- 1
startB <- m1$coefficients
startA <- 0.5

# Fit phylo logistic regression
fit_phylo <- phyloglm(
  decline ~ migratory * Hemisphere * flyway2 + EOO_log_cent,
  data   = data,
  phy    = tree,
  method = "logistic_IG10",
  btol   = 40,
  log.alpha.bound = 4,
  start.beta = startB,
  start.alpha = startA,
  boot = 1
)

#-----------------------------------------------------------------------------#
# Drop the 3 way interaction term

fit_phylo2 <- phyloglm(
  decline ~ migratory * Hemisphere + migratory * flyway2 + Hemisphere * flyway2 + EOO_log_cent,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 50,
  log.alpha.bound = 4,
  boot            = 100
)
summary(fit_phylo2)

saveRDS(fit_phylo2, "FitPhylo2.rds")
#-----------------------------------------------------------------------------#
# Try with combined Hem/Flyway

data <- data %>%
  mutate(HemFly = interaction(Hemisphere, flyway2, sep = "_", drop = TRUE))

with(data, table(migratory, HemFly))


bad <- names(which(tbl["1", ] == 0))
bad


m1 <- glm(decline ~ migratory * HemFly, family = binomial, data)
m1$coefficients[23] <- 1
m1$coefficients[25] <- 1
startB <- m1$coefficients
startA <- 0.5

fit_phyloHemFly <- phyloglm(
  decline ~ migratory * Hemisphere + migratory * HemFly + EOO_log_cent,
  data            = data,
  phy             = tree,
  method          = "logistic_MPLE",
  btol            = 60,
  log.alpha.bound = 4,
  boot            = 1
)
summary(fit_phylo2)

saveRDS(fit_phylo2, "FitPhylo2.rds")