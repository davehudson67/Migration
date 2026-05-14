library(phytools)
library(nimble)
library(MCMCglmm)

# Load dataset
data <- read.csv("IUCN_data_B.csv", header=T, na.strings="", stringsAsFactors = FALSE)

# Convert relevant columns to factors
fixedEffects <- c("decline", "migratory", "Hemisphere","Artificial", "Forest","Grassland", "Savanna",
                  "Shrubland", "Wetlands", "American", "AfroPal", "Asian")
data[fixedEffects] <- lapply(data[fixedEffects], factor)
data$EOO <- as.numeric(data$EOO)

# Log-transform and center EOO
data$EOO_log <- log(data$EOO + 1)
data$EOO_log_centered <- scale(data$EOO_log, center = TRUE, scale = FALSE)

# Re-leveling factors
data$migratory <- relevel(data$migratory, ref='0')
data$Forest <- relevel(data$Forest, ref='1')
data$American <- relevel(data$American, ref='1')
data$Hemisphere <- relevel(data$Hemisphere, ref='2')

# Load and process the phylogenetic tree
tree <- read.tree("ult_5k_tree.trees")
is.ultrametric(tree)
tree$edge.length[tree$edge.length <= 0] <- 0
tree <- di2multi(tree)  # Convert to a tree with zero-length branches collapsed

# Use the inverseA function to generate the A-inverse matrix
Ainv <- inverseA(tree)$Ainv 


code <- nimbleCode({
  for (i in 1:N) {
    # Likelihood for categorical response
    decline[i] ~ dcat(p[i, 1:K])  
    logit(p[i, 1:K]) <- beta0 + 
      beta_migratory * migratory[i] + 
      beta_American * American[i] + 
      beta_AfroPal * AfroPal[i] +
      beta_Asian * Asian[i] +
      beta_Hemisphere * Hemisphere[i] +
      beta_EOO * EOO_log_centered[i] +
      #beta_habitat * habitat[i] +  #
      phylo_effect[species[i]]  # Phylogenetic random effect
  }
  
  # Priors for fixed effects
  beta0 ~ dnorm(0, 0.001)
  beta_migratory ~ dnorm(0, 0.001)
  beta_American ~ dnorm(0, 0.001)
  beta_AfroPal ~ dnorm(0, 0.001)
  beta_Asian ~ dnorm(0, 0.001)
  beta_Hemisphere ~ dnorm(0, 0.001)
  beta_EOO ~ dnorm(0, 0.001)
  #beta_habitat ~ dnorm(0, 0.001)
  
  # Phylogenetic random effects
  phylo_effect[1:N_species] ~ dmnorm(rep(0, N_species), precision = Ainv)
  
  # Prior for phylogenetic variance
  sigma_phylo ~ dunif(0, 10)
  Ainv <- solve(sigma_phylo * A)
})


Ainv_dense <- as.matrix(Ainv)
dimnames(Ainv_dense) <- NULL

# Prepare data list for NIMBLE
modeldata <- list(decline = as.numeric(data$decline),
             migratory = as.numeric(data$migratory),
             American = as.numeric(data$American),
             AfroPal = as.numeric(data$AfroPal),
             Asian = as.numeric(data$Asian),
             Hemisphere = as.numeric(data$Hemisphere),
             EOO_log_centered = as.numeric(data$EOO_log_centered),
             habitat = as.numeric(data$decline),
             Ainv = as.numeric(Ainv_dense))

# Constants
constants <- list(N = 9544,
                  N_species = 8651,
                  K = 2)

# Build and compile the model
model <- nimbleModel(code, data = data, constants = constants)
mcmc_conf <- configureMCMC(model)

# You may add more monitors if necessary
# mcmc_conf$addMonitors(c('phylo_effect'))

mcmc <- buildMCMC(mcmc_conf)
compiled_model <- compileNimble(model)
compiled_mcmc <- compileNimble(mcmc, project = model)

# Run MCMC
samples <- runMCMC(compiled_mcmc, niter = 50000, nburnin = 10000, thin = 10)
