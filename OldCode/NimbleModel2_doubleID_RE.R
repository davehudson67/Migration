library(ape)
library(nimble)
library(MCMCglmm)
library(Matrix)
library(tidyverse)
library(coda)
library(mcmcplots)
set.seed(123)  

setwd("~/FraserCode")
data <- read.csv("Data/IUCN_data_B.csv", header=TRUE,
                 na.strings="", stringsAsFactors=FALSE)

# Remove exact duplicates??
data <- distinct(data, .keep_all = TRUE)
data <- distinct(data, animal, .keep_all = TRUE)

# Pick 100 species at random
sample_spp <- sample(unique(data$animal), size = 1000)

# subset your data to only those species
data <- data %>% filter(animal %in% sample_spp)

# Recode flyways
data$flyway <- data$American + 2 * data$AfroPal + 4 * data$Asian
data$flyway2 <- data$flyway
data$flyway2[data$flyway2 %in% c(3, 5, 6, 7)] <- 5
data$flyway2[data$flyway2 == 4] <- 3
data$flyway2[data$flyway2 == 5] <- 4

# factor categorical predictors
allFactors <- c("decline","migratory","Hemisphere","Artificial","Forest",
                "Grassland","Savanna","Shrubland","Wetlands",
                "American","AfroPal","Asian","flyway","flyway2")
data[allFactors] <- lapply(data[allFactors], factor)

# Numeric response
data$y <- as.integer(as.character(data$decline))

# Log‐transform & center EOO
data$EOO       <- as.numeric(data$EOO)
data$EOO_log   <- log(data$EOO + 1)
data$EOO_cent  <- scale(data$EOO_log, center=TRUE, scale=FALSE)


# Load phlogeny tree
fullTree   <- read.tree("Data/ult_5k_tree.trees")

# Exclude species not in the dataset
keep_sp    <- intersect(data$animal, fullTree$tip.label)
prunedTree <- drop.tip(fullTree, setdiff(fullTree$tip.label, keep_sp))
rm(keep_sp)

# Jitter any zero-length edges
eps <- 1e-8
prunedTree$edge.length <- prunedTree$edge.length + eps
Ainv <- inverseA(prunedTree, nodes = "TIPS")$Ainv

data$animalID <- as.numeric(factor(data$animal, levels = prunedTree$tip.label))
nAnimal <- length(unique(data$animal))

# Design matrix
X_mat  <- model.matrix(~ migratory * Hemisphere * flyway2 + EOO_cent, data)

# Constants, data, inits
constants <- list(
  nObs = nrow(data),
  p = ncol(X_mat),
  nAnimal = nAnimal,
  animalID = data$animalID,
  Ainv = as.matrix(Ainv),
  zeroMean = rep(0, nAnimal))

nim_data <- list(
  X = X_mat,
  y = data$y)

inits <- list(
  beta = rep(0, constants$p),
  tau_id = 1,
  tau_phy = 1,
  animal_raw = rep(0, nAnimal),
  animal_id = rep(0, nAnimal))

data$animalID <- as.numeric(factor(data$animal, levels = prunedTree$tip.label))
nAnimal <- length(unique(data$animal))

# NIMBLE code with direct Bernoulli–logit ––

phyloCode <- nimbleCode({
  # fixed effects
  for(j in 1:p) {
    beta[j] ~ dnorm(0, sd=10)
  }
  
  # two variance‐components (precisions)
  tau_phy   ~ dgamma(0.001, 0.001)
  sigma_phy  <- 1/sqrt(tau_phy)
  tau_id    ~ dgamma(0.001, 0.001)
  sigma_id   <- 1/sqrt(tau_id)
  
  # phylogenetic random effect
  animal_raw[1:nAnimal] ~ dmnorm(mean = zeroMean[1:nAnimal], prec = Ainv[1:nAnimal,1:nAnimal])
  
  # independent species effect
  for(k in 1:nAnimal) {
    animal_id[k] ~ dnorm(0, sd = sigma_id)
  }
  
  # scale the phylo effect
  for(k in 1:nAnimal) {
    animal[k] <- animal_raw[k] * sigma_phy
  }
  
  # likelihood
  for(i in 1:nObs) {
    logit(prob[i]) <- inprod(X[i,1:p], beta[1:p]) + animal[animalID[i]] + animal_id[animalID[i]]
    y[i] ~ dbern(prob[i])
  }
})

# Compile & run ––
Rmodel <- nimbleModel(
  code      = phyloCode,
  data      = nim_data,
  constants = constants,
  inits     = inits
)

Cmodel <- compileNimble(Rmodel)

conf  <- configureMCMC(Rmodel, monitors = c("beta","tau_phy", "tau_id"))
conf$removeSampler(c("animal_raw", "animal_id", "tau_phy", "tau_id", "beta"))
conf$addSampler(target = "animal_raw", type = "AF_slice", control = list(adapt = TRUE))
conf$addSampler(target = "animal_id", type = "AF_slice", control = list(adapt = TRUE))
conf$addSampler(target = "tau_phy", type = "slice")
conf$addSampler(target = "tau_id", type = "slice")
conf$addSampler(target = "beta", type = "AF_slice", control = list(adapt = TRUE))
conf$printSamplers()                    

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

samples <- runMCMC(
  Cmcmc,
  niter    = 50000,
  thin     = 100,
  nburnin  = 1000,
  nchains  = 2,
  setSeed  = TRUE
)

saveRDS(samples, "samples_NimbleModel_LATEST.rds")

# Inspect posterior
mcmcplot(samples, parms = c("beta", "lambda"))
