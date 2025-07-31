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

# Remove exact duplicates??
data <- distinct(data, .keep_all = TRUE)
data <- distinct(data, animal, .keep_all = TRUE)

#dupes <- data %>%
#  group_by(animal) %>%
#  filter(n() > 1) %>%
#  ungroup()
  
#length(unique(dupes$animal))

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

# or do it using ape package
#tree_M <- vcv(prunedTree) # get the variance/covariances
#A <- tree_M/max(tree_M) # scale
#Ainv <- solve(A)
#dim(Ainv)

# Jitter any zero-length edges
#eps <- 1e-8
#prunedTree$edge.length <- prunedTree$edge.length + eps

# Get the full node‐wise precision from inverseA()
#inv_res <- inverseA(prunedTree)
#Ainv <- as.matrix(inverseA(prunedTree)$Ainv)

# Get rid of the Node entries
#isTip <- !grepl("^Node", rownames(Ainv))
#Ainv <- as.matrix(Ainv[isTip, isTip])
#Ainv <- as(Ainv, "dgCMatrix")
#rm(isTip)

# Turn data$animal into a factor with levels exactly in the tree order
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
  beta       = rep(0, constants$p),
  tau        = 1,
  animal_raw = rep(0, nAnimal))

#–– 4) NIMBLE code with direct Bernoulli–logit ––

phyloCode <- nimbleCode({
  # fixed‐effect priors
  for(j in 1:p) {
    beta[j] ~ dnorm(0, sd=10)
  }
  
  # phylo‐variance prior
  tau ~ dgamma(0.001, 0.001)
  sigma_animal <- 1 / sqrt(tau)
  
  # raw phylogenetic effects
  animal_raw[1:nAnimal] ~ dmnorm(mean = zeroMean[1:nAnimal], prec = Ainv[1:nAnimal, 1:nAnimal])
  
  for(k in 1:nAnimal) {
    animal[k] <- animal_raw[k] * sigma_animal
  }
  
  # Bernoulli–logit likelihood
  for(i in 1:nObs) {
    logit(prob[i]) <- inprod(X[i, 1:p], beta[1:p]) + animal[animalID[i]]
    y[i] ~ dbern(prob[i])
  }
})

#–– 5) Compile & run ––

Rmodel <- nimbleModel(
  code      = phyloCode,
  data      = nim_data,
  constants = constants,
  inits     = inits
)

Cmodel <- compileNimble(Rmodel)

conf  <- configureMCMC(Rmodel,
                       monitors = c("beta","tau","sigma_animal","animal"))

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
mcmcplot(samples, parms = c("beta", "tau"))
