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
# First dataset
data<- read.csv("Data/IUCN_data_B.csv", header=T, na.strings="", stringsAsFactors = FALSE)

# Other datasets:
# data<- read.csv("IUCN_data_C.csv", header=T, na.strings="", stringsAsFactors = FALSE)
data$rowID <- paste0("sp_", seq_len(nrow(data)))

data$flyway <- data$American + 2 * data$AfroPal + 4 * data$Asian
data$flyway2 <- data$flyway
data <- data[!(data$flyway2 %in% c(6, 7) & data$Hemisphere==2), ]
data$flyway2[data$flyway2==3]<-2
data$flyway2[data$flyway2==4]<-3
data$flyway2[data$flyway2==5]<-4
data$flyway2[data$flyway2==6]<-5
data$flyway2[data$flyway2==7]<-6
data$HemFly <- paste("Hem", data$Hemisphere, "Fly", data$flyway2, sep="")

####
#### Prepare data
####

# Sort data structure - turn fixed effects into factors (from integers)
fixedEffects <- c("decline", "migratory", "Hemisphere","Artificial", "Forest","Grassland", "Savanna",
                  "Shrubland", "Wetlands", "American", "AfroPal", "Asian","flyway","flyway2")
data[fixedEffects] <- lapply(data[fixedEffects], factor)  
data$EOO <- as.numeric(data$EOO)

# check
str(data) 

# Apply log transformation (add 1 to avoid log(0) if necessary)
data$EOO_log <- log(data$EOO + 1)

# Center the log-transformed EOO variable
data$EOO_log_cent <- scale(data$EOO_log, center = TRUE, scale = FALSE)


####
#### Re-level factors to sort intercept issue
####

# Re-leveling to use 'South American forest specialists'
data$migratory <- relevel(data$migratory, ref='1')
data$Forest <- relevel(data$Forest, ref='1')
data$flyway <- relevel(data$flyway, ref='1')
data$flyway2 <- relevel(data$flyway2, ref='1')
data$Hemisphere <- factor(data$Hemisphere, levels=c('1','3','2'))

####
#### Phylogeny 
####

### Load and sort phylogenetic tree
#tree <- read.tree("ult_5k_tree.trees")
tree <-readRDS("FixedTree.rds")
# Check if tree is ultrametric
#is.ultrametric(tree)
#tree$edge.length[tree$edge.length <= 0] <- 0
#tree <- di2multi(tree)

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
rownames(data) <- data$rowID

# allow up to 2 GiB of globals
options(future.globals.maxSize = 2 * 1024^3)

# Fit phylo logistic regression

###This version is GLM-style but it makes it harder to judge the PR(decline) differences (I think) 
m1 <- phyloglm(
  decline ~ migratory * Hemisphere * flyway2 + EOO_cent,
  data   = data,
  phy    = tree,
  method = "logistic_MPLE",
  btol   = 20,
  log.alpha.bound = 4,
  boot = 100
)
summary(m1)


####This model useful for plotting, and gives clarity on each category-combination' Pr(decline)
m2 <- phyloglm(
  decline ~ migratory:HemFly -1 + EOO_log_centered,
  data   = data,
  phy    = tree,
  method = "logistic_MPLE",
  btol   = 20,
  log.alpha.bound = 4,
  boot = 0
)
summary(m2)

