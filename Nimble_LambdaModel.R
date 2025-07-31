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

# Pick 100 species at random
sample_spp <- sample(unique(data$animal), size = 500)

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

# 1. Read & prepare your pruned tree (only tips you use):
tree <- read.tree("Data/ult_5k_tree.trees")
tree$edge.length[tree$edge.length <= 0] <- 1e-8
keep_sp <- intersect(data$animal, tree$tip.label)
tree    <- drop.tip(tree, setdiff(tree$tip.label, keep_sp))

# 2. Get the Brownian‐motion correlation matrix among tips:
R <- vcv.phylo(tree, corr = TRUE)

# 3. Eigen‐decompose it once:
E  <- eigen(R, symmetric = TRUE)
U  <- E$vectors      # orthonormal matrix (n×n)
d  <- E$values       # eigenvalues (length n)

# 4. Drop tiny negative eigenvalues (numerical noise)
d[d < 0] <- 0

data$animalID <- as.numeric(factor(data$animal, levels = tree$tip.label))
nAnimal <- length(unique(data$animal))

# Design matrix
X_mat  <- model.matrix(~ migratory * Hemisphere * flyway2 + EOO_cent, data)

# pack data & constants
constants <- list(
  nObs     = nrow(data),
  p        = ncol(X_mat),
  nAnimal  = nAnimal,
  animalID = data$animalID,
  U        = U,
  d        = d,
  zeroMean = rep(0, nAnimal)
)

nim_data <- list(X = X_mat, y = data$y)

inits <- list(beta = rep(0, ncol(X_mat)),
              lambda = 0.5,
              tau_phy = 1,
              z = rnorm(nAnimal))

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
Rmodel <- nimbleModel(phyloCode, data = nim_data,
                      constants = constants, inits = inits)
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
                   niter = 10000, 
                   thin = 10,
                   nburnin = 1000,
                   nchains = 2)

saveRDS(samples, "Output/LambdaModel_n500.rds")
