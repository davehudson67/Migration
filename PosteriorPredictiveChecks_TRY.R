# Load the necessary libraries
library(coda)
library(ggplot2)

# Load the posterior samples from the MCMCglmm model
posterior_samples <- readRDS("MCMC_3.rds")  # Assuming you saved the model output
posterior_samples <- as.matrix(posterior_samples$Sol)

data <- read.csv("IUCN_data_B.csv", header=T, na.strings="", stringsAsFactors = FALSE)
# Apply log transformation (add 1 to avoid log(0) if necessary)
data$EOO_log <- log(data$EOO + 1)

# Center the log-transformed EOO variable
data$EOO_log_centered <- scale(data$EOO_log, center = TRUE, scale = FALSE)

# Check if there are any NAs in the relevant columns
any(is.na(data$migratory))
any(is.na(data$Hemisphere))
any(is.na(data$EOO_log_centered))

# Number of posterior samples
n_samples <- nrow(posterior_samples)

# Create storage for the simulated data (for posterior predictive check)
y_sim <- matrix(NA, nrow = n_samples, ncol = nrow(data))

# Loop through each posterior sample to generate predictions
for (i in 1:n_samples) {
  # For each sample, extract the parameters
  intercept <- posterior_samples[i, "(Intercept)"]
  migratory_effect <- posterior_samples[i, "migratory1"]
  hemisphere_effect <- posterior_samples[i, "Hemisphere1"]
  eoo_effect <- posterior_samples[i, "EOO_log_centered"]
  
  # Check if interactions are present
  migratory_hemisphere_interaction <- posterior_samples[i, "migratory1:Hemisphere1"]
  migratory_eoo_interaction <- posterior_samples[i, "migratory1:EOO_log_centered"]
  
  # Predict the log-odds of species decline for each row in the data
  log_odds <- intercept + 
    migratory_effect * data$migratory + 
    hemisphere_effect * (data$Hemisphere == 1) +  # Adjust for Hemisphere1
    eoo_effect * data$EOO_log_centered +
    migratory_hemisphere_interaction * data$migratory * (data$Hemisphere == 1) + 
    migratory_eoo_interaction * data$migratory * data$EOO_log_centered
  
  # Convert log-odds to probabilities
  prob <- exp(log_odds) / (1 + exp(log_odds))
  
  # Simulate binary outcomes (decline or not) based on these probabilities
  y_sim[i, ] <- rbinom(n = nrow(data), size = 1, prob = prob)
}

# Proportion of species in decline in the observed data
observed_decline <- mean(data$decline)

# Proportion of species in decline for each simulated dataset
simulated_decline <- apply(y_sim, 1, mean)

# Plot comparison of observed vs simulated data
hist(simulated_decline, breaks = 30, main = "Posterior Predictive Distribution of Decline",
     xlab = "Proportion of species in decline")
abline(v = observed_decline, col = "red", lwd = 2)  # Observed value

# Simulated mean decline probabilities
simulated_means <- apply(y_sim, 2, mean)

# Plot observed vs simulated probabilities of decline
plot(simulated_means, data$decline, 
     xlab = "Simulated Decline Probability", ylab = "Observed Decline",
     main = "Observed vs Simulated Decline")
abline(0, 1, col = "red")  # Add a reference line


