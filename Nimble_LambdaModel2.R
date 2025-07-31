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

#data[allFactors] <- lapply(data[allFactors], as.numeric)

# Numeric response
data$y <- as.integer(as.character(data$decline))

# Log‐transform & center EOO
data$EOO       <- as.numeric(data$EOO)
data$EOO_log   <- log(data$EOO + 1)
data$EOO_cent  <- scale(data$EOO_log, center=TRUE, scale=FALSE)

# Design matrix
X_mat  <- model.matrix(~ migratory * Hemisphere * flyway2 + EOO_cent, data)

#------------------------------------------------------------------------------#
# Read tree
#tree <- read.tree("Data/ult_5k_tree.trees")
tree <- readRDS("FixedTree.rds")
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
#rownames(data) <- data$rowID

# save out fixed tree
saveRDS(tree, "FixedTree.rds")

keep_sp <- intersect(data$rowID, tree$tip.label)
tree    <- drop.tip(tree, setdiff(tree$tip.label, keep_sp))

# Get the Brownian‐motion correlation matrix among tips:
R <- vcv.phylo(tree, corr = TRUE)

# Eigen‐decompose it once:
E  <- eigen(R, symmetric = TRUE)
U  <- E$vectors      # orthonormal matrix (n×n)
d  <- E$values       # eigenvalues (length n)

# 4. Drop tiny negative eigenvalues (numerical noise)
d[d < 0] <- 0

data$animalID <- as.numeric(factor(data$rowID, levels = tree$tip.label))
nAnimal <- length(unique(data$rowID))

# pack data & constants
constants <- list(
  nObs     = nrow(data),
  p        = ncol(X_mat),
  nAnimal  = nAnimal,
  animalID = data$animalID,
  U        = U,
  d        = d)

nim_data <- list(X = X_mat, y = data$y)

inits <- list(beta = rep(0, ncol(X_mat)),
              lambda = 0.5,
              tau_phy = 1,
              z = rep(0, nAnimal))

phyloCode <- nimbleCode({
  # Fixed effects
  for(j in 1:p) beta[j] ~ dnorm(0, sd = 10)
  
  # Pagel’s lambda
  lambda ~ dunif(0, 1)
  
  # Phylogenetic scale
  tau_phy   ~ dgamma(0.001, 0.001)
  sigma_phy <- 1/sqrt(tau_phy)
  
  # Latent standard normals
  for(i in 1:nAnimal) {
    z[i] ~ dnorm(0, sd = 1)
    # build the adjusted eigenvalues
    d_lambda[i] <- lambda * d[i] + (1 - lambda)
    sd_lambda[i] <- sqrt(d_lambda[i])
  }
  
  # Construct the phylogenetic effect via spectral trick:
  #animal = sigma_phy * U %*% (sd_lambda * z)
  for(i in 1:nAnimal) {
    # dot‐product of row i of U with elementwise sd_lambda*z
    animal[i] <- sigma_phy * inprod(U[i, 1:nAnimal], sd_lambda[1:nAnimal] * z[1:nAnimal])
  }
  
  # likelihood
  for(i in 1:nObs) {
    logit(prob[i]) <- inprod(X[i, 1:p], beta[1:p]) + animal[animalID[i]]
    y[i] ~ dbern(prob[i])
  }
})

# compile & configure
Rmodel <- nimbleModel(phyloCode, 
                      data = nim_data,
                      constants = constants,
                      inits = inits)


Cmodel <- compileNimble(Rmodel)
conf   <- configureMCMC(Rmodel,
                        monitors = c("beta","lambda","sigma_phy"),
                        useConjugacy = TRUE,
                        enableAutomatedBlocking = TRUE)

conf$removeSampler(c("beta", "lambda", "tau_phy", "z"))
conf$addSampler(target = "lambda", type = "slice")
conf$addSampler(target = "tau_phy", type = "slice")
conf$addSampler(target = "beta", type = "AF_slice", control = list(adapt = TRUE))
conf$addSampler(target = "z", type = "AF_slice", control = list(adapt = TRUE))
conf$printSamplers()      

Rmcmc  <- buildMCMC(conf)
Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(Cmcmc, 
                   niter = 50000, 
                   thin = 100,
                   nburnin = 9000,
                   nchains = 2)

saveRDS(samples, "Output/LambdaModel_n500.rds")


