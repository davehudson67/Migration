# Load required libraries
library(ape)
library(nimble)
library(dplyr)
library(MCMCglmm)
setwd("~/FraserCode")

# Read and clean the data
data <- read.csv("Data/IUCN_data_B.csv", header = TRUE, na.strings = "", stringsAsFactors = FALSE)

# Recode flyway variables
data$flyway  <- data$American + 2 * data$AfroPal + 4 * data$Asian
data$flyway2 <- data$flyway
data$flyway2[data$flyway2 %in% c(3, 5, 6, 7)] <- 5
data$flyway2[data$flyway2 == 4]           <- 3
data$flyway2[data$flyway2 == 5]           <- 4

# Turn selected columns into factors
fixedEffects <- c("decline","migratory","Hemisphere","Artificial","Forest","Grassland",
                  "Savanna","Shrubland","Wetlands","American","AfroPal","Asian","flyway","flyway2")
data[fixedEffects] <- lapply(data[fixedEffects], factor)

# Transform and center EOO
data$EOO            <- as.numeric(data$EOO)
data$EOO_log        <- log(data$EOO + 1)
data$EOO_log_centered <- as.numeric(scale(data$EOO_log, center = TRUE, scale = FALSE))

# Re‐level factors
data$migratory  <- relevel(data$migratory, ref = "0")
data$Forest     <- relevel(data$Forest,    ref = "1")
data$flyway     <- relevel(data$flyway,    ref = "1")
data$flyway2    <- relevel(data$flyway2,   ref = "1")
data$Hemisphere <- factor(data$Hemisphere, levels = c("1","3","2"))

# Convert response to integer 1/C for categorical
data$decline_int <- as.integer(data$decline)

# Load and process phylogeny
tree <- read.tree("Data/ult_5k_tree.trees")
if(!is.ultrametric(tree)) {
  tree$edge.length[tree$edge.length <= 0] <- 0
  tree <- di2multi(tree)
}

eps <- 1e-8

# Add epsilon to all branch lengths,
# which guarantees no zero lengths remain:
tree$edge.length <- tree$edge.length + eps

# Now inverseA should work:
Ainv_list <- inverseA(tree)
Ainv      <- Ainv_list$Ainv

# Prepare matrices and indices for NIMBLE
X_mat    <- model.matrix(~ migratory * Hemisphere * flyway2 + EOO_log_centered, data = data)
nObs     <- nrow(data)
p        <- ncol(X_mat)
animalID <- factor(data$animal)
nAnimal  <- nlevels(animalID)
animalIdx<- as.integer(animalID)

constants <- list(
  nObs      = nObs,
  p         = p,
  nAnimal   = nAnimal,
  animalID  = animalIdx,       
  Ainv      = as.matrix(Ainv), 
  zeroMean  = rep(0, nAnimal)           
)

nim_data <- list(
  y = data$decline_int,
  X = X_mat
)

inits <- list(
  beta   = rep(0, p),
  tau    = 1,
  animal = rep(0, nAnimal)
)

# Define the NIMBLE model
phyloCode <- nimbleCode({
  for (j in 1:p) {
    beta[j] ~ dnorm(0, sd = 10)
  }
  tau ~ dgamma(0.001, 0.001)
  
  # Build the phylogenetic precision matrix
  prec_mat[1:nAnimal, 1:nAnimal] <- tau * Ainv[1:nAnimal, 1:nAnimal]
  
  animal[1:nAnimal] ~ dmnorm(
    mean = zeroMean[1:nAnimal],
    prec = prec_mat[1:nAnimal, 1:nAnimal]
  )
  
  for (i in 1:nObs) {
    logit(prob[i]) <- inprod(X[i, 1:p], beta[1:p]) + animal[animalID[i]]
    y[i]           ~ dbern(prob[i])
  }
})

# 5. Build, compile, and run
Rmodel <- nimbleModel(code      = phyloCode,
                      data      = nim_data,
                      constants = constants,
                      inits     = inits)

Cmodel <- compileNimble(Rmodel)

conf  <- configureMCMC(Rmodel, monitors = c("beta","tau","animal"))

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc,
                   niter    = 50000,
                   thin     = 100,
                   nburnin  = 1000,
                   nchains = 2
                   )

