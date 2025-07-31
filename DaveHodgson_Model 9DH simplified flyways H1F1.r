

#####################################################################################################
#####################################################################################################


#############################################################################
#############################################################################
#############################################################################
###                                                                       ###
###                   IUCN data Analysis 2024: Model 2                    ###
###                                                                       ###
#############################################################################
#############################################################################
#############################################################################


#####################################################################################################
#####################################################################################################


##########################
##########################


# remove any old variables from the workspace
rm(list=ls()) 
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
data <- read.csv("Data/IUCN_data_B.csv", header=T, na.strings="", stringsAsFactors = FALSE)

data$flyway <- data$American + 2 * data$AfroPal + 4 * data$Asian
data$flyway2 <- data$flyway
data$flyway2[data$flyway2 %in% c(3, 5, 6, 7)] <- 5
data$flyway2[data$flyway2 == 4] <- 3
data$flyway2[data$flyway2 == 5] <- 4

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
data$EOO_log_centered <- scale(data$EOO_log, center = TRUE, scale = FALSE)


####
#### Re-level factors to sort intercept issue
####

# Re-leveling to use 'South American forest specialists'
data$migratory <- relevel(data$migratory, ref = '0')
data$Forest <- relevel(data$Forest, ref = '1')
data$flyway <- relevel(data$flyway, ref = '1')
data$flyway2 <- relevel(data$flyway2, ref = '1')
data$Hemisphere <- factor(data$Hemisphere, levels = c('1','3','2'))

####
#### Phylogeny 
####

### Load and sort phylogenetic tree
tree <- read.tree("Data/ult_5k_tree.trees")
# Check if tree is ultrametric
is.ultrametric(tree)
tree$edge.length[tree$edge.length <=0]<-0
tree<-di2multi(tree)



####
#### MCMCglmm 
####

###
### Model 1: decline ~ migratory * flyway * hemisphere + EOO
###

MCMC_9DH <- MCMCglmm(decline ~ migratory * Hemisphere * flyway2 + EOO_log_centered,
                         data=data, family = "categorical",
                         random=~animal, ginverse=list(animal=inverseA(tree)$Ainv),
                         prior=list(R=list(V=1, fix=1),
                                    G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1))),
                         pl=TRUE, nitt=5000,thin=100,burnin=1000, verbose =TRUE)


summary.MCMCglmm(MCMC_9DH)
levels(data$decline)

##
## Save model outputs
##

# Export R workspace
save(MCMC_9DH, file = "MCMC_9DH.Rdata")
# Export model as object
saveRDS(MCMC_9DH, "MCMC_9DH.rds")

# Export plots
#dev.off()
par(mfrow=c(8,2), mar=c(2,2,1,0))
pdf(file = "MCMC_9DH.pdf")
plot(MCMC_9DH$Sol, auto.layout=F)
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

oneModel <- clean.MCMC(MCMC_9DH)  # get all the info from summary(modelName)
oneModel$modelName <- getName.MCMC(MCMC_9DH)  # add the model's name in a new column
oneModel  # check out the created dataframe

dataList <- list(MCMC_9DH)

dataListNames <- list("MCMC_9DH")

readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

mcmcOutputs <- as.data.frame(do.call(rbind, readyList), stringsAsFactors = FALSE)

write.csv (mcmcOutputs, "MCMC_9DH.csv")

###get the geographic estimates of Pr(decline)
H1F1M0<-MCMC_9DH$Sol[,1]
H3F1M0<-apply(MCMC_9DH$Sol[,c(1,3)],1,sum)
H2F1M0<-apply(MCMC_9DH$Sol[,c(1,4)],1,sum)
H1F2M0<-apply(MCMC_9DH$Sol[,c(1,5)],1,sum)
H3F2M0<-apply(MCMC_9DH$Sol[,c(1,3,5,14)],1,sum)
H2F2M0<-apply(MCMC_9DH$Sol[,c(1,4,5,15)],1,sum)
H1F3M0<-apply(MCMC_9DH$Sol[,c(1,6)],1,sum)
H3F3M0<-apply(MCMC_9DH$Sol[,c(1,3,6,16)],1,sum)
H2F3M0<-apply(MCMC_9DH$Sol[,c(1,4,6,17)],1,sum)
H1F4M0<-apply(MCMC_9DH$Sol[,c(1,7)],1,sum)
H3F4M0<-apply(MCMC_9DH$Sol[,c(1,3,7,18)],1,sum)
H2F4M0<-apply(MCMC_9DH$Sol[,c(1,4,7,19)],1,sum)
H1F1M1<-apply(MCMC_9DH$Sol[,c(1,2)],1,sum)
H3F1M1<-apply(MCMC_9DH$Sol[,c(1,3,2,9)],1,sum)
H2F1M1<-apply(MCMC_9DH$Sol[,c(1,4,2,10)],1,sum)
H1F2M1<-apply(MCMC_9DH$Sol[,c(1,5,2,11)],1,sum)
H3F2M1<-apply(MCMC_9DH$Sol[,c(1,3,5,14,2,9,11,20)],1,sum)
H2F2M1<-apply(MCMC_9DH$Sol[,c(1,4,5,15,2,10,11,21)],1,sum)
H1F3M1<-apply(MCMC_9DH$Sol[,c(1,6,2,12)],1,sum)
H3F3M1<-apply(MCMC_9DH$Sol[,c(1,3,6,16,2,9,12,22)],1,sum)
H2F3M1<-apply(MCMC_9DH$Sol[,c(1,4,6,17,2,10,12,23)],1,sum)
H1F4M1<-apply(MCMC_9DH$Sol[,c(1,7,2,13)],1,sum)
H3F4M1<-apply(MCMC_9DH$Sol[,c(1,3,7,18,2,9,13,24)],1,sum)
H2F4M1<-apply(MCMC_9DH$Sol[,c(1,4,7,19,2,10,13)],1,sum)

pdf("model9 H1F1 intercept.pdf")
opar<-par(mfrow=c(3,4))
plot(density(H1F1M0),main="",xlab="H1F1")
lines(density(H1F1M1),col="red")
plot(density(H1F2M0),main="",xlab="H1F2")
lines(density(H1F2M1),col="red")
plot(density(H1F3M0),main="",xlab="H1F3")
lines(density(H1F3M1),col="red")
plot(density(H1F4M0),main="",xlab="H1F4")
lines(density(H1F4M1),col="red")
plot(density(H3F1M0),main="",xlab="H3F1")
lines(density(H3F1M1),col="red")
plot(density(H3F2M0),main="",xlab="H3F2")
lines(density(H3F2M1),col="red")
plot(density(H3F3M0),main="",xlab="H3F3")
lines(density(H3F3M1),col="red")
plot(density(H3F4M0),main="",xlab="H3F4")
#lines(density(H3F4M1),col="red")
plot(density(H2F1M0),main="",xlab="H2F1")
lines(density(H2F1M1),col="red")
plot(density(H2F2M0),main="",xlab="H2F2")
lines(density(H2F2M1),col="red")
plot(density(H2F3M0),main="",xlab="H2F3")
lines(density(H2F3M1),col="red")
plot(density(H2F4M0),main="",xlab="H2F4")
lines(density(H2F4M1),col="red")
par(opar)
dev.off()

###get the estimates of delta-Pr(decline)
dH1F1<-H1F1M1-H1F1M0
dH1F2<-H1F2M1-H1F2M0
dH1F3<-H1F3M1-H1F3M0
dH1F4<-H1F4M1-H1F4M0
dH2F1<-H2F1M1-H2F1M0
dH2F2<-H2F2M1-H2F2M0
dH2F3<-H2F3M1-H2F3M0
dH2F4<-H2F4M1-H2F4M0
dH3F1<-H3F1M1-H3F1M0
dH3F2<-H3F2M1-H3F2M0
dH3F3<-H3F3M1-H3F3M0
dH3F4<-H3F4M1-H3F4M0

pdf("model9 migration diff H1F1 intercept.pdf")
opar<-par(mfrow=c(3,4))
plot(density(dH1F1),main=paste("pMCMC = ",round(length(dH1F1[dH1F1<0])/length(dH1F1),3)),xlab="H1F1")
lines(c(0,0),c(0,2),col="red")
plot(density(dH1F2),main=paste("pMCMC = ",round(length(dH1F2[dH1F2<0])/length(dH1F2),3)),xlab="H1F2")
lines(c(0,0),c(0,2),col="red")
plot(density(dH1F3),main=paste("pMCMC = ",round(length(dH1F3[dH1F3<0])/length(dH1F3),3)),xlab="H1F3")
lines(c(0,0),c(0,2),col="red")
plot(density(dH1F4),main=paste("pMCMC = ",round(length(dH1F4[dH1F4<0])/length(dH1F4),3)),xlab="H1F4")
lines(c(0,0),c(0,2),col="red")
plot(density(dH3F1),main=paste("pMCMC = ",round(length(dH3F1[dH3F1<0])/length(dH3F1),3)),xlab="H3F1")
lines(c(0,0),c(0,2),col="red")
plot(density(dH3F2),main=paste("pMCMC = ",round(length(dH3F2[dH3F2<0])/length(dH3F2),3)),xlab="H3F2")
lines(c(0,0),c(0,2),col="red")
plot(density(dH3F3),main=paste("pMCMC = ",round(length(dH3F3[dH3F3<0])/length(dH3F3),3)),xlab="H3F3")
lines(c(0,0),c(0,2),col="red")
#plot(density(dH3F4),paste("pMCMC = ",round(length(dH3F4[dH3F4<0])/length(dH3F4),3)),xlab="H3F4")
plot(density(dH3F4),col="white",paste(""),xlab="H3F4")
lines(c(0,0),c(0,2),col="red")
plot(density(dH2F1),main=paste("pMCMC = ",round(length(dH2F1[dH2F1<0])/length(dH2F1),3)),xlab="H2F1")
lines(c(0,0),c(0,2),col="red")
plot(density(dH2F2),main=paste("pMCMC = ",round(length(dH2F2[dH2F2<0])/length(dH2F2),3)),xlab="H2F2")
lines(c(0,0),c(0,2),col="red")
plot(density(dH2F3),main=paste("pMCMC = ",round(length(dH2F3[dH2F3<0])/length(dH2F3),3)),xlab="H2F3")
lines(c(0,0),c(0,2),col="red")
plot(density(dH2F4),main=paste("pMCMC = ",round(length(dH2F4[dH2F4<0])/length(dH2F4),3)),xlab="H2F4")
lines(c(0,0),c(0,2),col="red")
par(opar)
dev.off()
