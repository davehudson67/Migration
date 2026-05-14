

#####################################################################################################
#####################################################################################################


#############################################################################
#############################################################################
#############################################################################
###                                                                       ###
###                   IUCN data Analysis 2024: Model 1                    ###
###                                                                       ###
#############################################################################
#############################################################################
#############################################################################


#####################################################################################################
#####################################################################################################


##########################
##########################


# remove any old variables from the workspace
#rm(list=ls()) 
# Check by listing if any variables still in the workspace
#ls()          


#########################
# Required packages

library(phytools)
library(MCMCglmm)
library(phangorn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

options(max.print=1000000)


# setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024")


# First dataset
#data<- read.csv("IUCN_data_A.csv", header=T, na.strings="", stringsAsFactors = FALSE)

# Other datasets:
 data<- read.csv("IUCN_data_B.csv", header=T, na.strings="", stringsAsFactors = FALSE)
# data<- read.csv("IUCN_data_C.csv", header=T, na.strings="", stringsAsFactors = FALSE)


####
#### Prepare data
####

# Sort data structure - turn fixed effects into factors (from integers)
fixedEffects <- c("decline", "migratory", "Hemisphere","Artificial", "Forest","Grassland", "Savanna",
                  "Shrubland", "Wetlands", "American", "AfroPal", "Asian")
data[fixedEffects] <- lapply(data[fixedEffects], factor)  
data$EOO <- as.numeric(data$EOO)

# check
str(data) 

# Apply log transformation (add 1 to avoid log(0) if necessary)
data$EOO_log <- log(data$EOO + 1)

# Center the log-transformed EOO variable
data$EOO_log_centered <- scale(data$EOO_log, center = TRUE, scale = FALSE)


####
#### Re-level factors to sort intercept issue
####

# Re-leveling to use 'South American forest specialists'
data$migratory<-relevel(data$migratory, ref='0')
data$Forest<-relevel(data$Forest, ref='1')
data$American<-relevel(data$American, ref='1')
data$Hemisphere<-relevel(data$Hemisphere, ref='2')

####
#### Phylogeny 
####

### Load and sort phylogenetic tree
tree <- read.tree("ult_5k_tree.trees")
# Check if tree is ultrametric
is.ultrametric(tree)
tree$edge.length[tree$edge.length <=0]<-0
tree<-di2multi(tree)



####
#### MCMCglmm 
####

###
### Model 1: decline ~ migratory * flyway * hemisphere * habitat * EOO
###

MCMC_1 <- MCMCglmm(decline ~ migratory*(American + AfroPal + Asian) + 
                           migratory*Hemisphere + 
                           migratory*EOO_log_centered +
                           migratory*(Artificial + Forest + Grassland + Savanna + Shrubland + Wetlands) + 
                           (American + AfroPal + Asian)*Hemisphere + 
                           (American + AfroPal + Asian)*(Artificial + Forest + Grassland + Savanna + Shrubland + Wetlands) + 
                           (American + AfroPal + Asian)*EOO_log_centered + 
                           Hemisphere*(Artificial + Forest + Grassland + Savanna + Shrubland + Wetlands) +
                           Hemisphere*EOO_log_centered +
                           (Artificial + Forest + Grassland + Savanna + Shrubland + Wetlands)*EOO_log_centered,
                         data=data, family = "categorical",
                         random=~animal, ginverse=list(animal=inverseA(tree)$Ainv),
                         prior=list(R=list(V=1, fix=1),
                                    G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1))),
                         pl=TRUE, nitt=5000000,thin=100,burnin=1000, verbose =TRUE)


summary.MCMCglmm(MCMC_1)


##
## Save model outputs
##

# Export R workspace
save(MCMC_1, file = "MCMC_1_B.Rdata")
# Export model as object
saveRDS(MCMC_1, "MCMC_1_B.rds")

# Export plots
dev.off()
par(mfrow=c(8,2), mar=c(2,2,1,0))
pdf(file = "MCMC_1_B.pdf")
plot(MCMC_1$Sol, auto.layout=F)
dev.off()


##
## Export 'sol' as csv (less essential)
##

# https://github.com/tmalsburg/MCMCglmm-intro

clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # add the model name

oneModel <- clean.MCMC(MCMC_1)  # get all the info from summary(modelName)
oneModel$modelName <- getName.MCMC(MCMC_1)  # add the model's name in a new column
oneModel  # check out the created dataframe

dataList <- list(MCMC_1)

dataListNames <- list("MCMC_1_B")

readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

mcmcOutputs <- as.data.frame(do.call(rbind, readyList), stringsAsFactors = FALSE)

write.csv (mcmcOutputs, "MCMC_1_B.csv")

run <- readRDS("MCMC_1.rds")

library(Mcmc)
