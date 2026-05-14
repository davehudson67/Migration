
# Start with empty workspace
rm(list=ls(all=TRUE)) 


# Load packages
library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(hrbrthemes)
library(tibble)
library(ggtext)
library(bayesplot)
library(MCMCvis)
library(patchwork)
library(ggmcmc)
library(effectsize)
library(coda)

### Set Working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024")

## Check that unknown pop trend species are unbiased
dataold<- read.csv("Previous_IUCN_dataset.csv", header=T, na.strings="", stringsAsFactors = FALSE)
datanew<- read.csv("Final_IUCN_dataset.csv", header=T, na.strings="", stringsAsFactors = FALSE)
shared_animals <- intersect(dataold$animal, datanew$animal)
data <- dataold[!(dataold$animal %in% shared_animals), ]

# write.csv(data, "UnknownPopTrends.csv")

head(data)

# Using matrix() for a concise approach:
table <- matrix(c(359, 461, 195, 625, 360, 460), nrow = 3, ncol = 2, byrow = TRUE)
rownames(table) <- c("American", "AfroPal", "Asian")
colnames(table) <- c("Present", "Not Present")

# Print the table to verify its structure:
print(table)

# Apply the chisq.test() function:
chisq.test(table)



#


###################################################
###################################################





### Set Working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1")

# load ("MCMC_0_Oct.RData") wrong data, ignore
# load ("MCMC_0.RData")
# load ("MCMC_1.RData")

# Function to turn into probability on logit scale
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


###################################################
###################################################

# New S1

data<- read.csv("Summary_S1.csv", header=T, na.strings="", stringsAsFactors = FALSE)

data$Migratory     <- as.factor(data$Migratory)
data$Decline <- as.factor(data$Decline)

#contingency_table <- table(data$decline, data$migratory)

# Reorder the levels of the "Migratory" variable
data <- data %>%
  mutate(Migratory = reorder(Migratory, -as.numeric(Migratory)))

# Define custom colors for each combination of "Migratory" and "Decline"
custom_colors <- c(
    "#f2e3c0", "#a1ddd8",           
    "#E1BE6A", "#40B0A6" 
)




ggplot(data, aes(x = Migratory, y = Prop, fill = interaction(Migratory, Decline))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  scale_x_discrete(labels = c("Non-migratory \n(7,984)", "Migratory \n(1,560)")) +  # Specify custom x-axis labels
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  )


###################################################
###################################################


#### Plot of phylogenetic signal

# calculate heritability (i.e. phylogenetic signal) for probit models (see: http://devillemereuil.legtux.org/wp-content/uploads/2012/12/tuto_en.pdf)
herit_probit<-function(MCMC_1){ # x = an MCMCglmm model object
  MCMC_1$VCV[,1]/(MCMC_1$VCV[,1]+MCMC_1$VCV[,2]+1) # phylogenetic variance divided by total variance, +1 due to additional variance of the probit link
}

# calculate mean heritability + 95% credible intervals
mean(herit_probit(MCMC_0))
quantile(herit_probit(MCMC_1), c(0.025, 0.975))
heritMCMC<- herit_probit(MCMC_1) #can be then plotted as a density plot if we want.

PM <- posterior.mode(heritMCMC)

ggplot(as.data.frame(heritMCMC), aes(x = heritMCMC))+
  geom_density(fill="#F4BAB8", alpha=0.4, size = 0.8) +
  
  geom_segment(x = posterior.mode(heritMCMC), xend = posterior.mode(heritMCMC), y = -1.5, yend = 45, colour = "#001219", size = 1.2) +
  geom_segment(x = HPDinterval(heritMCMC)[[1]], xend = HPDinterval(heritMCMC)[[1]], y = -1.5, yend = 45, colour = "#001219", size = 1.2, linetype = "dotted") +
  geom_segment(x = HPDinterval(heritMCMC)[[2]], xend = HPDinterval(heritMCMC)[[2]], y = -1.5, yend = 45, colour = "#001219", size = 1.2, linetype = "dotted") +

  xlab("\nPhylogenetic signal ")+
  ylab(" Density\n")+
  ylim(0,42) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



###################################################
###################################################


#### Figure 1


# Get data from R object
data<-as.data.frame(MCMC_0[["Sol"]])

# Rename columns
names(data)[1] <- "nonmigratory"
names(data)[2] <- "mig"
# Create migratory values (by adding to intercept)
data$migratory <- data$nonmigratory + data$mig
# Turn values into probability
data$nonmigratory <- logit2prob(data$nonmigratory)
data$migratory <- logit2prob(data$migratory)

# Turn df from wide to long
dataM <- as.data.frame(data$migratory)
dataM$term <- "migratory"
names(dataM)[1] <- "value"
medianM <- median(dataM$value)

dataNM <- as.data.frame(data$nonmigratory)
dataNM$term <- "nonmigratory"
names(dataNM)[1] <- "value"
medianNM <- median(dataNM$value)

data<- rbind(dataM,dataNM)

# Then subset by term
dataM <- subset(data, data$term== "migratory")
dataM <- dataM[order(dataM$value),]

dataNM <- subset(data, data$term== "nonmigratory")
dataNM <- dataNM[order(dataNM$value),]

# Plot
F1b <- ggplot(data) +
  xlim(0,1) +
  ylim(0.5,4) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.025)], xmax=dataM$value[round(nrow(dataM)*0.975)],
           ymin=1, ymax=2, fill="#40B0A6", alpha=0.25) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.25)], xmax=dataM$value[round(nrow(dataM)*0.75)],
           ymin=1, ymax=2, fill="#40B0A6", alpha=0.5) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.45)], xmax=dataM$value[round(nrow(dataM)*0.55)],
           ymin=1, ymax=2, fill="#40B0A6", alpha=1) +
  annotate("rect", xmin=dataNM$value[round(nrow(dataNM)*0.025)], xmax=dataNM$value[round(nrow(dataNM)*0.975)],
           ymin=2.5, ymax=3.5, fill="#E1BE6A", alpha=0.25) +
  annotate("rect", xmin=dataNM$value[round(nrow(dataNM)*0.25)], xmax=dataNM$value[round(nrow(dataNM)*0.75)],
           ymin=2.5, ymax=3.5, fill="#E1BE6A", alpha=0.5) +
  annotate("rect", xmin=dataNM$value[round(nrow(dataNM)*0.45)], xmax=dataNM$value[round(nrow(dataNM)*0.55)],
           ymin=2.5, ymax=3.5, fill="#E1BE6A", alpha=1) +
  geom_segment(x=medianM, xend=medianM, y=1, yend=2, size = 2) +
  geom_segment(x=medianNM, xend=medianNM, y=2.5, yend=3.5, size = 2) +
  xlab("\n Probability of decline") +
  ggtitle("b)") +
  scale_x_continuous(labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3), labels = c("Migratory", "Non-migratory"))



# Calculate quantiles for both datasets upfront
quantiles_M <- quantile(dataM$value, probs = c(0.025, 0.25, 0.45, 0.55, 0.75, 0.975))
quantiles_NM <- quantile(dataNM$value, probs = c(0.025, 0.25, 0.45, 0.55, 0.75, 0.975))

# Create the plot with optimized annotations
F1b <- ggplot(data) +
  xlim(0, 1) +
  ylim(0.5, 4) +
  annotate("rect", xmin = quantiles_M[1], xmax = quantiles_M[6], ymin = 1, ymax = 2, fill = "#40B0A6", alpha = 0.25) +
  annotate("rect", xmin = quantiles_M[2], xmax = quantiles_M[5], ymin = 1, ymax = 2, fill = "#40B0A6", alpha = 0.5) +
  annotate("rect", xmin = quantiles_M[3], xmax = quantiles_M[4], ymin = 1, ymax = 2, fill = "#40B0A6", alpha = 1) +
  annotate("rect", xmin = quantiles_NM[1], xmax = quantiles_NM[6], ymin = 2.5, ymax = 3.5, fill = "#E1BE6A", alpha = 0.25) +
  annotate("rect", xmin = quantiles_NM[2], xmax = quantiles_NM[5], ymin = 2.5, ymax = 3.5, fill = "#E1BE6A", alpha = 0.5) +
  annotate("rect", xmin = quantiles_NM[3], xmax = quantiles_NM[4], ymin = 2.5, ymax = 3.5, fill = "#E1BE6A", alpha = 1) +
  geom_segment(x = medianM, xend = medianM, y = 1, yend = 2, size = 2) +
  geom_segment(x = medianNM, xend = medianNM, y = 2.5, yend = 3.5, size = 2) +
  xlab("\n Probability of decline") +
  ggtitle("b)") +
  scale_x_continuous(labels = c("0.00" = "0%", "0.25" = "25%", "0.50" = "50%", "0.75" = "75%", "1.00" = "100%")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(breaks = c(1.5, 3), labels = c("Migratory", "Non-migratory"))


##
### Plot parameter estimates
##

# Use GGMCMC package to plot using GG graphics
mc <- MCMC_0[["Sol"]]
S <- ggs(mc)

S1 <- subset(S, Parameter == "migratory1")


#F1a <- ggs_caterpillar(S1) +
#  geom_segment(x=0, xend=0, y=-Inf, yend=Inf, size = 1, colour = "black", linetype = 2, alpha = 0.7) +
#  aes(color = Parameter) +
#  scale_colour_manual(values = c("#40B0A6")) +
#  scale_y_discrete(labels=c("Migratory")) +
#  xlim(-0.9,0.2) +
#  ggtitle("a)") +
#  ylab("") +
#  xlab("\n Parameter estimate") +
#  theme_bw() +
#  theme(panel.grid.major.x = element_blank(), 
#        panel.grid.minor = element_blank(),
#        legend.position="none",
#        plot.title = element_text(size=18,face="bold"),
#        axis.text = element_text(size=14,face="bold"),
#        plot.margin = margin(1,1,1,1, "cm"),
#        axis.title = element_text(size=16,face="bold"))




# 2024 Fig 1a



# Calculate quantiles for both datasets upfront
quantiles <- quantile(S1$value, probs = c(0.025, 0.25, 0.45, 0.55, 0.75, 0.975))
median <- median(S1$value)

# Create the plot with optimized annotations
F1a <- ggplot(data) +
  scale_x_continuous(limits = c(-0.9, 0.3), breaks = c(-0.75, -0.50, -0.25, 0.00, 0.25)) +
  scale_y_continuous(limits = c(0.5, 4), breaks = c(2.25), labels = c("Difference")) +
  annotate("rect", xmin = quantiles[1], xmax = quantiles[6], ymin = 1.5, ymax = 3, fill = "#40B0A6", alpha = 0.25) +
  annotate("rect", xmin = quantiles[2], xmax = quantiles[5], ymin = 1.5, ymax = 3, fill = "#40B0A6", alpha = 0.5) +
  annotate("rect", xmin = quantiles[3], xmax = quantiles[4], ymin = 1.5, ymax = 3, fill = "#40B0A6", alpha = 1) +
  geom_segment(x = median, xend = median, y = 1.5, yend = 3, size = 2) +
  geom_segment(x=0, xend=0, y=-2, yend=5, size = 1, colour = "black", linetype = 2, alpha = 0.7) +
  ggtitle("a)") +
  xlab("\nLog-odds of parameter estimate") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))


# F1a | F1b  


# Set the working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024/graphics")

# Combine the plots using patchwork
combined_plot <- F1a | F1b

# Save the plot as a PNG file with a resolution of 300 dpi
#png("Fig1.png", width = 12, height = 6, units = "in", res = 300)
#print(combined_plot)
#dev.off()



#





# Violin plot

F1c <- ggplot(S1, aes(x = "", y = value)) +
  geom_violin(fill = "#40B0A6", alpha = 0.6, colour = "#40B0A6", size = 1.1)+
#  geom_segment(x=0, xend=0, y=-Inf, yend=Inf, size = 1, colour = "black", linetype = 2, alpha = 0.7) +
  xlim(-0.9,0.2) +
  coord_flip() +
  scale_x_discrete(labels=c("Migratory")) +
  ggtitle("a)") +
  xlab("") +
  ylab("\n Parameter estimate") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none",
        plot.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=14,face="bold"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title = element_text(size=16,face="bold"))


F1a | F1c





###################################################
###################################################



## Check how in decline the intercept category is
logit2prob(2.996309379) #95% as in text



###################################################
###################################################



#### Figure 2

plot.new()

# load ("MCMC_1.RData")

# For SD, need coda then:
mcmcTrace <- mcmc(MCMC_1[["Sol"]])
# Convert the summary output to a data frame
mcmcSummary <- summary(mcmcTrace)

# Extract summary statistics
summary_data <- as.data.frame(mcmcSummary$statistics)


# Write to a CSV file
# write.csv(summary_data, file = "mcmc_summary.csv")

# Get data from R object
# data<-as.data.frame(MCMC_1[["Sol"]])
# data <- as.mcmc(data)


##
### 2a
##

# Use GGMCMC package to plot using GG graphics
mc <- MCMC_1[["Sol"]]
S <- ggs(mc)

S1 <- subset(S, Parameter == "migratory1" 
             | Parameter == "Artificial1"
             | Parameter == "Savanna1"
             | Parameter == "Forest0"
             | Parameter == "Shrubland1"
             | Parameter == "Hemisphere1"
             | Parameter == "Hemisphere5"
             | Parameter == "Asian1"
             #| Parameter == "migratory1:Asian1"
             | Parameter == "migratory1:Savanna1"
             | Parameter == "migratory1:Hemisphere1"
             | Parameter == "Asian1:Artificial1"
             | Parameter == "Asian1:Savanna1")



F2b <- ggs_caterpillar(S1) +
  geom_segment(x=0, xend=0, y=-Inf, yend=Inf, size = 1, colour = "black", linetype = 2, alpha = 0.7) +
  aes(color = Parameter) +

  scale_y_discrete(limits = rev(c("migratory1", "Hemisphere1",
                                  "Hemisphere5","Forest0", "Shrubland1",
                                  "Savanna1",  "Artificial1", 
                                  "Asian1", 
                                  #"migratory1:Asian1", 
                                  "migratory1:Savanna1", "migratory1:Hemisphere1",
                                  "Asian1:Artificial1", "Asian1:Savanna1")),
                   
                   labels = c("migratory1" = "Migratory species",
                              "Artificial1" = "Artificial habitat",
                              "Savanna1" = "Savanna habitat", 
                              "Forest0" = "Non-forest habitats", 
                              "Shrubland1" = "Shrubland habitat", 
                              "Hemisphere1" = "Northern-only species",
                              "Hemisphere5" = "Widespread species", 
                              "Asian1" = "Asian system", 
                              #"migratory1:Asian1" = "Migratory : Asian", 
                              "migratory1:Savanna1" = "Migratory : Savanna", 
                              "migratory1:Hemisphere1" = "Migratory : Northern-only",
                              "Asian1:Artificial1" = "Asian : Artificial", 
                              "Asian1:Savanna1" = "Asian: Savanna")) +
  
  scale_colour_manual(values = c("migratory1" = "#40B0A6",
                                 "Artificial1" = "#C4900F",
                                 "Savanna1" = "#FECC4A", 
                                 "Forest0" = "lightgoldenrod2", 
                                 "Shrubland1" = "lightgoldenrod2", 
                                 "Hemisphere1" = "#30A8DF",
                                 "Hemisphere5" = "#a2d6f9", 
                                 "Asian1" = "#941b0c", 
                                 #"migratory1:Asian1" = "#264653", 
                                 "migratory1:Savanna1" = "#264653", 
                                 "migratory1:Hemisphere1" = "#264653",
                                 "Asian1:Artificial1" = "#264653", 
                                 "Asian1:Savanna1" = "#264653")) +
  
  geom_text(x= 3.1, y = 12, label = "**", colour = "#40B0A6", size = 7) +
  geom_text(x= 3.1, y = 10, label = "*", colour = "#a2d6f9", size = 7) +
  geom_text(x= 3.1, y = 9, label = "*", colour = "lightgoldenrod2", size = 7) +
  geom_text(x= 3.1, y = 8, label = "***", colour = "lightgoldenrod2", size = 7) +
  geom_text(x= 3.1, y = 7, label = "*", colour = "#FECC4A", size = 7) +
  geom_text(x= 3.1, y = 6, label = "***", colour = "#C4900F", size = 7) +
  geom_text(x= 3.1, y = 4, label = "*", colour = "#264653", size = 7) +
  geom_text(x= 3.1, y = 3, label = "*", colour = "#264653", size = 7) +
  geom_text(x= 3.1, y = 2, label = "*", colour = "#264653", size = 7) +
  geom_text(x= 3.1, y = 1, label = "*", colour = "#264653", size = 7) +
  
  xlim(-6,3.2) +
  ggtitle("b)") +
  ylab("") +
  xlab("\n Parameter estimate") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none",
        plot.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=14,face="bold"),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.title = element_text(size=16,face="bold"))


##
### 2c
##

# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_2c = subset(data, select = c("(Intercept)", "migratory1", "Hemisphere1", "migratory1:Hemisphere1"))

# Rename columns
names(data_2c)[1] <- "intercept"
names(data_2c)[2] <- "mig"
names(data_2c)[3] <- "north"
names(data_2c)[4] <- "mig_north"

# Create correct values (by adding to intercept)
data_2c$migratory <- data_2c$intercept + data_2c$mig
data_2c$northern <- data_2c$intercept + data_2c$north
data_2c$mn_interaction <- data_2c$intercept + data_2c$mig + data_2c$north + data_2c$mig_north

# Subset to keep new data
data_2c = subset(data_2c, select = -c(1,2,3,4))

# Turn values into probability
data_2c$migratory <- logit2prob(data_2c$migratory)
data_2c$northern <- logit2prob(data_2c$northern)
data_2c$mn_interaction <- logit2prob(data_2c$mn_interaction)

# Turn df from wide to long
dataM <- as.data.frame(data_2c$migratory)
dataM$term <- "migratory"
names(dataM)[1] <- "value"
median1 <- median(dataM$value)

dataN <- as.data.frame(data_2c$northern)
dataN$term <- "northern"
names(dataN)[1] <- "value"
median2 <- median(dataN$value)

dataMN <- as.data.frame(data_2c$mn_interaction)
dataMN$term <- "interaction"
names(dataMN)[1] <- "value"
median3 <- median(dataMN$value)

data<- rbind(dataM,dataN, dataMN)


# Then subset by term
dataM <- subset(data, data$term== "migratory")
dataM <- dataM[order(dataM$value),]

dataN <- subset(data, data$term== "northern")
dataN <- dataN[order(dataN$value),]

dataMN <- subset(data, data$term== "interaction")
dataMN <- dataMN[order(dataMN$value),]


# Plot
F2c <-  ggplot(data) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.025)], xmax=dataM$value[round(nrow(dataM)*0.975)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=0.4) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.25)], xmax=dataM$value[round(nrow(dataM)*0.75)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=0.6) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.45)], xmax=dataM$value[round(nrow(dataM)*0.55)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=1) +
  
  annotate("rect", xmin=dataMN$value[round(nrow(dataMN)*0.025)], xmax=dataMN$value[round(nrow(dataMN)*0.975)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.4) +
  annotate("rect", xmin=dataMN$value[round(nrow(dataMN)*0.25)], xmax=dataMN$value[round(nrow(dataMN)*0.75)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.6) +
  annotate("rect", xmin=dataMN$value[round(nrow(dataMN)*0.45)], xmax=dataMN$value[round(nrow(dataMN)*0.55)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=1) +
  
  annotate("rect", xmin=dataN$value[round(nrow(dataN)*0.025)], xmax=dataN$value[round(nrow(dataN)*0.975)],
           ymin=1, ymax=2, fill="#30A8DF", alpha=0.4) +
  annotate("rect", xmin=dataN$value[round(nrow(dataN)*0.25)], xmax=dataN$value[round(nrow(dataN)*0.75)],
           ymin=1, ymax=2, fill="#30A8DF", alpha=0.6) +
  annotate("rect", xmin=dataN$value[round(nrow(dataN)*0.45)], xmax=dataN$value[round(nrow(dataN)*0.55)],
           ymin=1, ymax=2, fill="#30A8DF", alpha=1) +
  
  geom_segment(x=median1, xend=median1, y=4, yend=5, size = 2) +
  geom_segment(x=median3, xend=median3, y=2.5, yend=3.5, size = 2) +
  geom_segment(x=median2, xend=median2, y=1, yend=2, size = 2) +

  xlab("\n Probability of decline") +
  ggtitle("c)") +
  scale_x_continuous(limits = c(0, 1),
    labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3, 4.5), labels = c("Northern-only\n species", "Northern-only\n migratory\n species", "Migratory\n species"))


F2c


##
### 2d
##

# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_2d = subset(data, select = c("(Intercept)", "migratory1", "Savanna1", "migratory1:Savanna1"))

# Rename columns
names(data_2d)[1] <- "intercept"
names(data_2d)[2] <- "mig"
names(data_2d)[3] <- "Savanna"
names(data_2d)[4] <- "mig_savanna"

# Create correct values (by adding to intercept)
data_2d$migratory <- data_2d$intercept + data_2d$mig
data_2d$savanna <- data_2d$intercept + data_2d$Savanna
data_2d$interaction <- data_2d$intercept + data_2d$mig + data_2d$Savanna + data_2d$mig_savanna

# Subset to keep new data
data_2d = subset(data_2d, select = -c(1,2,3,4))

# Turn values into probability
data_2d$migratory <- logit2prob(data_2d$migratory)
data_2d$savanna <- logit2prob(data_2d$savanna)
data_2d$interaction <- logit2prob(data_2d$interaction)

# Turn df from wide to long
dataM <- as.data.frame(data_2d$migratory)
dataM$term <- "migratory"
names(dataM)[1] <- "value"
median1 <- median(dataM$value)

dataS <- as.data.frame(data_2d$savanna)
dataS$term <- "savanna"
names(dataS)[1] <- "value"
median2 <- median(dataS$value)

dataI <- as.data.frame(data_2d$interaction)
dataI$term <- "interaction"
names(dataI)[1] <- "value"
median3 <- median(dataI$value)

data<- rbind(dataM,dataS,dataI)


# Then subset by term
dataM <- subset(data, data$term== "migratory")
dataM <- dataM[order(dataM$value),]

dataS <- subset(data, data$term== "savanna")
dataS <- dataS[order(dataS$value),]

dataI <- subset(data, data$term== "interaction")
dataI <- dataI[order(dataI$value),]


# Plot
F2d <- ggplot(data) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.025)], xmax=dataM$value[round(nrow(dataM)*0.975)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=0.4) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.25)], xmax=dataM$value[round(nrow(dataM)*0.75)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=0.6) +
  annotate("rect", xmin=dataM$value[round(nrow(dataM)*0.45)], xmax=dataM$value[round(nrow(dataM)*0.55)],
           ymin=4, ymax=5, fill="#40B0A6", alpha=1) +
  
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.025)], xmax=dataS$value[round(nrow(dataS)*0.975)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=0.4) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.25)], xmax=dataS$value[round(nrow(dataS)*0.75)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=0.6) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.45)], xmax=dataS$value[round(nrow(dataS)*0.55)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=1) +
  
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.025)], xmax=dataI$value[round(nrow(dataI)*0.975)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.4) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.25)], xmax=dataI$value[round(nrow(dataI)*0.75)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.6) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.45)], xmax=dataI$value[round(nrow(dataI)*0.55)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=1) +
  
  geom_segment(x=median1, xend=median1, y=4, yend=5, size = 2) +
  geom_segment(x=median2, xend=median2, y=1, yend=2, size = 2) +
  geom_segment(x=median3, xend=median3, y=2.5, yend=3.5, size = 2) +
  
  xlab("\n Probability of decline") +
  ggtitle("d)") +
  scale_x_continuous(limits = c(0, 1),
                     labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3, 4.5), labels = c("Savanna\n habitat", "Migratory\n species in\n savanna\n habitat", "Migratory\n species"))


# F2b | F2c | F2d



# Set the working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024/graphics")

# Combine the plots using patchwork
combined_plot <- F2b | F2c | F2d

# Save the plot as a PNG file with a resolution of 300 dpi
#png("Fig2.png", width = 18, height = 6, units = "in", res = 300)
#print(combined_plot)
#dev.off()









##
### Alternative figures c - humped line plot... looks shite
##

# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_2c = subset(data, select = c("(Intercept)", "migratory1", "Hemisphere1", "migratory1:Hemisphere1"))

# Rename columns
names(data_2c)[1] <- "intercept"
names(data_2c)[2] <- "mig"
names(data_2c)[3] <- "north"
names(data_2c)[4] <- "mig_north"

# Create correct values (by adding to intercept)
data_2c$migratory <- data_2c$intercept + data_2c$mig
data_2c$northern <- data_2c$intercept + data_2c$north
data_2c$mn_interaction <- data_2c$intercept + data_2c$mig + data_2c$north + data_2c$mig_north

# Subset to keep new data
data_2c = subset(data_2c, select = -c(1,2,3,4))

# Turn values into probability
#data_2c$migratory <- logit2prob(data_2c$migratory)
#data_2c$northern <- logit2prob(data_2c$northern)
#data_2c$mn_interaction <- logit2prob(data_2c$mn_interaction)


# Turn df from wide to long
dataM <- as.data.frame(data_2c$migratory)
dataM$term <- "migratory"
names(dataM)[1] <- "value"
median1 <- median(dataM$value)

dataN <- as.data.frame(data_2c$northern)
dataN$term <- "northern"
names(dataN)[1] <- "value"
median2 <- median(dataN$value)

dataMN <- as.data.frame(data_2c$mn_interaction)
dataMN$term <- "interaction"
names(dataMN)[1] <- "value"
median3 <- median(dataMN$value)

data<- rbind(dataM,dataN, dataMN)

#data_2c$Value <- logit2prob(data_2c$Value)


ggplot(data)  +
  geom_density(aes(x=value, colour=term), size = 1) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  ggtitle ("c)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))








###################################################
###################################################



#### Figure 3


# load ("MCMC_1.RData")


##
### 3a
##

# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_3a = subset(data, select = c("(Intercept)", "Asian1", "Artificial1", "Asian1:Artificial1"))

# Rename columns
names(data_3a)[1] <- "intercept"
names(data_3a)[2] <- "asia"
names(data_3a)[3] <- "artificial"
names(data_3a)[4] <- "asia_artificial"

# Create correct values (by adding to intercept)
data_3a$Asia <- data_3a$intercept + data_3a$asia
data_3a$Arti <- data_3a$intercept + data_3a$artificial
data_3a$interaction <- data_3a$intercept + data_3a$asia + data_3a$artificial + data_3a$asia_artificial

# Subset to keep new data
data_3a = subset(data_3a, select = -c(1,2,3,4))

# Turn values into probability
data_3a$Asia <- logit2prob(data_3a$Asia)
data_3a$Arti <- logit2prob(data_3a$Arti)
data_3a$interaction <- logit2prob(data_3a$interaction)

# Turn df from wide to long
dataA <- as.data.frame(data_3a$Asia)
dataA$term <- "Asia"
names(dataA)[1] <- "value"
median1 <- median(dataA$value)

dataArt <- as.data.frame(data_3a$Arti)
dataArt$term <- "Arti"
names(dataArt)[1] <- "value"
median2 <- median(dataArt$value)

dataI <- as.data.frame(data_3a$interaction)
dataI$term <- "interaction"
names(dataI)[1] <- "value"
median3 <- median(dataI$value)

data<- rbind(dataA,dataArt, dataI)


# Then subset by term
dataA <- subset(data, data$term== "Asia")
dataA <- dataA[order(dataA$value),]

dataArt <- subset(data, data$term== "Arti")
dataArt <- dataArt[order(dataArt$value),]

dataI <- subset(data, data$term== "interaction")
dataI <- dataI[order(dataI$value),]


# Plot

F3a <-  ggplot(data) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.025)], xmax=dataA$value[round(nrow(dataA)*0.975)],
           ymin=4, ymax=5, fill="#941b0c", alpha=0.4) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.25)], xmax=dataA$value[round(nrow(dataA)*0.75)],
           ymin=4, ymax=5, fill="#941b0c", alpha=0.6) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.45)], xmax=dataA$value[round(nrow(dataA)*0.55)],
           ymin=4, ymax=5, fill="#941b0c", alpha=1) +
  
  annotate("rect", xmin=dataArt$value[round(nrow(dataArt)*0.025)], xmax=dataArt$value[round(nrow(dataArt)*0.975)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.4) +
  annotate("rect", xmin=dataArt$value[round(nrow(dataArt)*0.25)], xmax=dataArt$value[round(nrow(dataArt)*0.75)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.6) +
  annotate("rect", xmin=dataArt$value[round(nrow(dataArt)*0.45)], xmax=dataArt$value[round(nrow(dataArt)*0.55)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=1) +
  
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.025)], xmax=dataI$value[round(nrow(dataI)*0.975)],
           ymin=1, ymax=2, fill="#C4900F", alpha=0.4) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.25)], xmax=dataI$value[round(nrow(dataI)*0.75)],
           ymin=1, ymax=2, fill="#C4900F", alpha=0.6) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.45)], xmax=dataI$value[round(nrow(dataI)*0.55)],
           ymin=1, ymax=2, fill="#C4900F", alpha=1) +
  
  geom_segment(x=median1, xend=median1, y=4, yend=5, size = 2) +
  geom_segment(x=median2, xend=median2, y=2.5, yend=3.5, size = 2) +
  geom_segment(x=median3, xend=median3, y=1, yend=2, size = 2) +
  
  xlab("\n Probability of decline") +
  ggtitle("a)") +
  scale_x_continuous(limits = c(0, 1),
                     labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3, 4.5), labels = c("Artificial\n habitat", "Asian species in \n Artificial habitat", "Asian\n species"))


# F3a




##
### 3b
##


# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_3b = subset(data, select = c("(Intercept)", "Asian1", "Savanna1", "Asian1:Savanna1"))

# Rename columns
names(data_3b)[1] <- "intercept"
names(data_3b)[2] <- "asia"
names(data_3b)[3] <- "savanna"
names(data_3b)[4] <- "asia_savanna"

# Create correct values (by adding to intercept)
data_3b$Asia <- data_3b$intercept + data_3b$asia
data_3b$Sava <- data_3b$intercept + data_3b$savanna
data_3b$interaction <- data_3b$intercept + data_3b$asia + data_3b$savanna + data_3b$asia_savanna

# Subset to keep new data
data_3b = subset(data_3b, select = -c(1,2,3,4))

# Turn values into probability
data_3b$Asia <- logit2prob(data_3b$Asia)
data_3b$Sava <- logit2prob(data_3b$Sava)
data_3b$interaction <- logit2prob(data_3b$interaction)

# Turn df from wide to long
dataA <- as.data.frame(data_3b$Asia)
dataA$term <- "Asia"
names(dataA)[1] <- "value"
median1 <- median(dataA$value)

dataS <- as.data.frame(data_3b$Sava)
dataS$term <- "Sava"
names(dataS)[1] <- "value"
median2 <- median(dataS$value)

dataI <- as.data.frame(data_3b$interaction)
dataI$term <- "interaction"
names(dataI)[1] <- "value"
median3 <- median(dataI$value)

data<- rbind(dataA,dataS, dataI)


# Then subset by term
dataA <- subset(data, data$term== "Asia")
dataA <- dataA[order(dataA$value),]

dataS <- subset(data, data$term== "Sava")
dataS <- dataS[order(dataS$value),]

dataI <- subset(data, data$term== "interaction")
dataI <- dataI[order(dataI$value),]


# Plot

F3b <-  ggplot(data) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.025)], xmax=dataA$value[round(nrow(dataA)*0.975)],
           ymin=4, ymax=5, fill="#941b0c", alpha=0.4) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.25)], xmax=dataA$value[round(nrow(dataA)*0.75)],
           ymin=4, ymax=5, fill="#941b0c", alpha=0.6) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.45)], xmax=dataA$value[round(nrow(dataA)*0.55)],
           ymin=4, ymax=5, fill="#941b0c", alpha=1) +
  
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.025)], xmax=dataS$value[round(nrow(dataS)*0.975)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.4) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.25)], xmax=dataS$value[round(nrow(dataS)*0.75)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=0.6) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.45)], xmax=dataS$value[round(nrow(dataS)*0.55)],
           ymin=2.5, ymax=3.5, fill="#264653", alpha=1) +
  
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.025)], xmax=dataI$value[round(nrow(dataI)*0.975)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=0.4) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.25)], xmax=dataI$value[round(nrow(dataI)*0.75)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=0.6) +
  annotate("rect", xmin=dataI$value[round(nrow(dataI)*0.45)], xmax=dataI$value[round(nrow(dataI)*0.55)],
           ymin=1, ymax=2, fill="#FECC4A", alpha=1) +
  
  geom_segment(x=median1, xend=median1, y=4, yend=5, size = 2) +
  geom_segment(x=median2, xend=median2, y=2.5, yend=3.5, size = 2) +
  geom_segment(x=median3, xend=median3, y=1, yend=2, size = 2) +
  
  xlab("\n Probability of decline") +
  ggtitle("b)") +
  scale_x_continuous(limits = c(0, 1),
                     labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3, 4.5), labels = c("Savanna\n habitat", "Asian species in \n Savanna habitat", "Asian\n species"))




# F3a | F3b



# Set the working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024/graphics")

# Combine the plots using patchwork
combined_plot <- F3a | F3b

# Save the plot as a PNG file with a resolution of 300 dpi
# png("Fig3.png", width = 12, height = 6, units = "in", res = 300)
# print(combined_plot)
# dev.off()








####### Figure S2

##
### Habitats and widespread in one figure
##

# Get data from R object
data<-as.data.frame(MCMC_1[["Sol"]])

# Subset necessary data
data_S2 = subset(data, select = c("(Intercept)", "Hemisphere5", "Forest0", "Shrubland1", "Savanna1", "Artificial1"))

# Rename columns
names(data_S2)[1] <- "intercept"
names(data_S2)[2] <- "Widespread"
names(data_S2)[3] <- "Not_Forest"
names(data_S2)[4] <- "Shrubland"
names(data_S2)[5] <- "Savanna"
names(data_S2)[6] <- "Artificial"

# Create correct values (by adding to intercept)
data_S2$Widespread <- data_S2$intercept + data_S2$Widespread
data_S2$Not_Forest <- data_S2$intercept + data_S2$Not_Forest
data_S2$Shrubland <- data_S2$intercept + data_S2$Shrubland
data_S2$Savanna <- data_S2$intercept + data_S2$Savanna
data_S2$Artificial <- data_S2$intercept + data_S2$Artificial

# Turn values into probability
data_S2$Widespread <- logit2prob(data_S2$Widespread)
data_S2$Not_Forest <- logit2prob(data_S2$Not_Forest)
data_S2$Shrubland <- logit2prob(data_S2$Shrubland)
data_S2$Savanna <- logit2prob(data_S2$Savanna)
data_S2$Artificial <- logit2prob(data_S2$Artificial)

# Turn df from wide to long
data1 <- as.data.frame(data_S2$Widespread)
data1$term <- "Widespread"
names(data1)[1] <- "value"
median1 <- median(data1$value)

data2 <- as.data.frame(data_S2$Not_Forest)
data2$term <- "Not_Forest"
names(data2)[1] <- "value"
median2 <- median(data2$value)

data3 <- as.data.frame(data_S2$Shrubland)
data3$term <- "Shrubland"
names(data3)[1] <- "value"
median3 <- median(data3$value)

data4 <- as.data.frame(data_S2$Savanna)
data4$term <- "Savanna"
names(data4)[1] <- "value"
median4 <- median(data4$value)

data5 <- as.data.frame(data_S2$Artificial)
data5$term <- "Artificial"
names(data5)[1] <- "value"
median5 <- median(data5$value)

data<- rbind(data1, data2, data3, data4, data5)


# Then subset by term
dataW <- subset(data, data$term== "Widespread")
dataW <- dataW[order(dataW$value),]

dataF <- subset(data, data$term== "Not_Forest")
dataF <- dataF[order(dataF$value),]

dataB <- subset(data, data$term== "Shrubland")
dataB <- dataB[order(dataB$value),]

dataS <- subset(data, data$term== "Savanna")
dataS <- dataS[order(dataS$value),]

dataA <- subset(data, data$term== "Artificial")
dataA <- dataA[order(dataA$value),]


# Plot
FS2 <-  ggplot(data) +
  annotate("rect", xmin=dataW$value[round(nrow(dataW)*0.025)], xmax=dataW$value[round(nrow(dataW)*0.975)],
           ymin=7, ymax=8, fill="#30A8DF", alpha=0.4) +
  annotate("rect", xmin=dataW$value[round(nrow(dataW)*0.25)], xmax=dataW$value[round(nrow(dataW)*0.75)],
           ymin=7, ymax=8, fill="#30A8DF", alpha=0.6) +
  annotate("rect", xmin=dataW$value[round(nrow(dataW)*0.45)], xmax=dataW$value[round(nrow(dataW)*0.55)],
           ymin=7, ymax=8, fill="#30A8DF", alpha=1) +
  
  annotate("rect", xmin=dataF$value[round(nrow(dataF)*0.025)], xmax=dataF$value[round(nrow(dataF)*0.975)],
           ymin=5.5, ymax=6.5, fill="lightgoldenrod2", alpha=0.4) +
  annotate("rect", xmin=dataF$value[round(nrow(dataF)*0.25)], xmax=dataF$value[round(nrow(dataF)*0.75)],
           ymin=5.5, ymax=6.5, fill="lightgoldenrod2", alpha=0.6) +
  annotate("rect", xmin=dataF$value[round(nrow(dataF)*0.45)], xmax=dataF$value[round(nrow(dataF)*0.55)],
           ymin=5.5, ymax=6.5, fill="lightgoldenrod2", alpha=1) +
  
  annotate("rect", xmin=dataB$value[round(nrow(dataB)*0.025)], xmax=dataB$value[round(nrow(dataB)*0.975)],
           ymin=4, ymax=5, fill="lightgoldenrod2", alpha=0.4) +
  annotate("rect", xmin=dataB$value[round(nrow(dataB)*0.25)], xmax=dataB$value[round(nrow(dataB)*0.75)],
           ymin=4, ymax=5, fill="lightgoldenrod2", alpha=0.6) +
  annotate("rect", xmin=dataB$value[round(nrow(dataB)*0.45)], xmax=dataB$value[round(nrow(dataB)*0.55)],
           ymin=4, ymax=5, fill="lightgoldenrod2", alpha=1) +
  
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.025)], xmax=dataS$value[round(nrow(dataS)*0.975)],
           ymin=2.5, ymax=3.5, fill="#FECC4A", alpha=0.4) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.25)], xmax=dataS$value[round(nrow(dataS)*0.75)],
           ymin=2.5, ymax=3.5, fill="#FECC4A", alpha=0.6) +
  annotate("rect", xmin=dataS$value[round(nrow(dataS)*0.45)], xmax=dataS$value[round(nrow(dataS)*0.55)],
           ymin=2.5, ymax=3.5, fill="#FECC4A", alpha=1) +
  
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.025)], xmax=dataA$value[round(nrow(dataA)*0.975)],
           ymin=1, ymax=2, fill="#C4900F", alpha=0.4) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.25)], xmax=dataA$value[round(nrow(dataA)*0.75)],
           ymin=1, ymax=2, fill="#C4900F", alpha=0.6) +
  annotate("rect", xmin=dataA$value[round(nrow(dataA)*0.45)], xmax=dataA$value[round(nrow(dataS)*0.55)],
           ymin=1, ymax=2, fill="#C4900F", alpha=1) +
  
  geom_segment(x=median1, xend=median1, y=7, yend=8, size = 2) +
  geom_segment(x=median2, xend=median2, y=5.5, yend=6.5, size = 2) +
  geom_segment(x=median3, xend=median3, y=4, yend=5, size = 2) +
  geom_segment(x=median4, xend=median4, y=2.5, yend=3.5, size = 2) +
  geom_segment(x=median5, xend=median5, y=1, yend=2, size = 2) +
  
  xlab("\n Probability of decline") +
  scale_x_continuous(limits = c(0, 1),
                     labels = c("0.00" = "0","0.25" = "0.25", "0.50" = "0.5","0.75" = "0.75","1.00" = "1")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))+
  scale_y_continuous(breaks=c(1.5, 3, 4.5, 6, 7.5), labels = c("Artificial habitat", "Savanna habitat", "Shrubland habitat", "Not in forest habitats", "Widespread species"))


# FS2



# Set the working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024/graphics")

# Save the plot as a PNG file with a resolution of 300 dpi
# png("FigS3.png", width = 6, height = 6, units = "in", res = 300)
# print(FS2)
# dev.off()









####### Figure S1

##
### Histograms!
##


data<-as.data.frame(MCMC_1[["Sol"]])

# Migratory
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$migratory1
data1$Type <- c("Migratory")
data1 <- data1[,c(83:84)]

H1 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#40B0A6")) +
  ggtitle ("Migratory species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



# Northern-only
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Hemisphere1
data1$Type <- c("Hemisphere1")
data1 <- data1[,c(83:84)]

H2 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#30A8DF")) +
  ggtitle ("Northern-only species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



# Widespread
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Hemisphere5
data1$Type <- c("Hemisphere5")
data1 <- data1[,c(83:84)]

H3 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#a2d6f9")) +
  ggtitle ("Widespread species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



# Not in Forest
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Forest0
data1$Type <- c("Forest0")
data1 <- data1[,c(83:84)]

H4 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("lightgoldenrod2")) +
  ggtitle ("Not in forest habitats") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Shrubland
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Shrubland1
data1$Type <- c("Shrubland1")
data1 <- data1[,c(83:84)]

H5 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("lightgoldenrod2")) +
  ggtitle ("Shrubland habitats") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Savanna
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Savanna1
data1$Type <- c("Savanna1")
data1 <- data1[,c(83:84)]

H6 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#FECC4A")) +
  ggtitle ("Savanna habitats") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



# Artificial
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Artificial1
data1$Type <- c("Artificial1")
data1 <- data1[,c(83:84)]

H7 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#C4900F")) +
  ggtitle ("Artificial habitats") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Asian
data1<- data
names(data1)[1] <- "Intercept"

data1$Value <- data1$Intercept + data1$Asian1
data1$Type <- c("Asian1")
data1 <- data1[,c(83:84)]

H8 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#941b0c")) +
  ggtitle ("Asian migratory system") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Migratory Savanna
data1<- data
names(data1)[1] <- "Intercept"
names(data1)[26] <- "Interest"

data1$Value <- data1$Intercept + data1$migratory1 + data1$Savanna1 + data1$Interest
data1$Type <- c("Migratory1: Savanna1")
data1 <- data1[,c(83:84)]

H9 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#264653")) +
  ggtitle ("Migratory : Savanna") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))



# Migratory northern only
data1<- data
names(data1)[1] <- "Intercept"
names(data1)[19] <- "Interest"

data1$Value <- data1$Intercept + data1$migratory1 + data1$Hemisphere1 + data1$Interest
data1$Type <- c("Migratory1: Hemisphere1")
data1 <- data1[,c(83:84)]

H10 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#264653")) +
  ggtitle ("Migratory : Northern-only") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Asian Artificial
data1<- data
names(data1)[1] <- "Intercept"
names(data1)[53] <- "Interest"

data1$Value <- data1$Intercept + data1$Asian1 + data1$Artificial1 + data1$Interest
data1$Type <- c("Asian1: Artificial1")
data1 <- data1[,c(83:84)]

H11 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#264653")) +
  ggtitle ("Asian : Artificial") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Asian Savanna
data1<- data
names(data1)[1] <- "Intercept"
names(data1)[56] <- "Interest"

data1$Value <- data1$Intercept + data1$Asian1 + data1$Savanna1 + data1$Interest
data1$Type <- "Asian1: Savanna1"
data1 <- data1[,c(83:84)]

H12 <- data1 %>%
  ggplot( aes(x=Value, fill=Type)) +
  geom_histogram(alpha=0.85) +
  geom_vline(xintercept=c(0),linetype="solid",size = 1.2) +
  scale_fill_manual(values=c("#264653")) +
  ggtitle ("Asian : Savanna") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size=14,face="bold"),
        plot.title = element_text(size=18,face="bold"),
        axis.title = element_text(size=16,face="bold"))


# Set the working directory
setwd("C:/Users/fraserbell/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/PhD Manuscripts/Chap 1/2024/graphics")

# Combine the plots using patchwork
combined_plot <- (H1|H2|H3) / (H4|H5|H6) / (H7|H8|H9) / (H10|H11|H12)

# Save the plot as a PNG file with a resolution of 300 dpi
png("FigS2.png", width = 18, height = 18, units = "in", res = 300)
print(combined_plot)
dev.off()


