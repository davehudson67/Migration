# Load the necessary libraries
library(coda)
library(ggplot2)

# Load the posterior samples from the MCMCglmm model
posterior_samples <- readRDS("MCMC_2.rds")  # Assuming you saved the model output
posterior_samples <- as.matrix(posterior_samples$Sol)

colnames(posterior_samples)

# Function to calculate the individual effect for a given factor, accounting for interactions
calculate_individual_effect <- function(samples, factor, interactions = NULL) {
  # Start with the intercept for each sample
  effect <- samples[, "(Intercept)"]
  
  # Add the main effect for the factor, if it exists in the samples
  if (factor %in% colnames(samples)) {
    effect <- effect + samples[, factor]
  }
  
  # Add interaction terms if provided and present in the samples
  if (!is.null(interactions) && length(interactions) > 0) {
    for (interaction in interactions) {
      if (interaction %in% colnames(samples)) {
        effect <- effect + samples[, interaction]
      }
    }
  }
  
  return(effect)
}

# Reference Category is...
# Non migratory
# Forest habitat (not included in model 2)
# American flyway (American, AfroPal, Asian)
# Hemisphere 2 (1, 2, 3)

# Define the factors and their relevant interaction terms for model 2
factors <- list(
  
  EOO_Hem1_nonMigratory = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered")),
  
  EOO_Hem1_Migratory = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered", "migratory1",
                                                  "migratory1:EOO_log_centered", "migratory1:Hemisphere1")),
  EOO_Hem2_nonMigratory = list("EOO_log_centered", NULL),
  
  EOO_Hem2_Migratory = list("EOO_log_centered", c("migratory1")),
  
  EOO_Hem3_nonMigratory = list("EOO_log_centered", c("Hemisphere3", "Hemisphere3:EOO_log_centered")),
  
  EOO_Hem3_Migratory = list("EOO_log_centered", c("Hemisphere3", "Hemisphere3:EOO_log_centered", "migratory1",
                                                  "migratory3:EOO_log_centered", "migratory3:Hemisphere1"))
)

factors <- list(
  
  EOO_Hem1_nonMigratory_AfroPal = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered", "American0", "AfroPal1",
                                                             "Hemisphere1:American0", "Hemisphere1:AfroPal1", "EOO_log_centered:American0",
                                                             "EOO_log_centered:AfroPal1")),
  
  EOO_Hem1_Migratory_AfroPal = list("EOO_log_centered", c("Hemisphere1", "migratory1", "American0", "AfroPal1", 
                                                          "Hemisphere1:EOO_log_centered", "Hemisphere1:American0", "Hemisphere1:AfroPal1", 
                                                          "migratory1:Hemisphere1", "migratory1:America0", "migratory1:AfroPal1", "migratory1:EOO_log_centered",
                                                          "EOO_log_centered:American0", "EOO_log_centered:AfroPal1")),
  
  EOO_Hem1_nonMigratory_Asian = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered", "American0", "Asian1",
                                                             "Hemisphere1:American0", "Hemisphere1:Asian1", "EOO_log_centered:American0",
                                                             "EOO_log_centered:Asian1")),
  
  EOO_Hem1_Migratory_Asian = list("EOO_log_centered", c("Hemisphere1", "migratory1", "American0", "Asian1", 
                                                          "Hemisphere1:EOO_log_centered", "Hemisphere1:American0", "Hemisphere1:Asian1", 
                                                          "migratory1:Hemisphere1", "migratory1:America0", "migratory1:Asian1", "migratory1:EOO_log_centered",
                                                          "EOO_log_centered:American0", "EOO_log_centered:Asian1")),
  
  EOO_Hem1_nonMigratory_American = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered")),
  
  EOO_Hem1_Migratory_American = list("EOO_log_centered", c("Hemisphere1", "migratory1", 
                                                        "Hemisphere1:EOO_log_centered", 
                                                        "migratory1:Hemisphere1", "migratory1:EOO_log_centered")),
  
  EOO_Hem2_nonMigratory_AfroPal = list("EOO_log_centered", c("American0", "AfroPal1", 
                                                             "EOO_log_centered:American0", "EOO_log_centered:AfroPal1")),
  
  EOO_Hem2_Migratory_AfroPal = list("EOO_log_centered", c("migratory1", "American0", "AfroPal1", 
                                                          "migratory1:American0", "migratory1:AfroPal1", "migratory1:EOO_log_centered",
                                                          "EOO_log_centered:American0", "EOO_log_centered:AfroPal1")),
  
  EOO_Hem2_nonMigratory_Asian = list("EOO_log_centered", c("American0", "Asian1", 
                                                           "EOO_log_centered:American0", "EOO_log_centered:Asian1")),
  
  EOO_Hem2_Migratory_Asian = list("EOO_log_centered", c("migratory1", "American0", "Asian1", 
                                                        "migratory1:American0", "migratory1:Asian1", "migratory1:EOO_log_centered",
                                                        "EOO_log_centered:American0", "EOO_log_centered:Asian1")),
  
  EOO_Hem2_nonMigratory_American = list("EOO_log_centered", NULL),
  
  EOO_Hem2_Migratory_American = list("EOO_log_centered", c("migratory1", "migratory1:EOO_log_centered")),
  
  EOO_Hem3_nonMigratory_AfroPal = list("EOO_log_centered", c("Hemisphere3", "Hemisphere3:EOO_log_centered", "American0", "AfroPal1",
                                                             "Hemisphere3:American0", "Hemisphere3:AfroPal1", "EOO_log_centered:American0",
                                                             "EOO_log_centered:AfroPal1")),
  
  EOO_Hem3_Migratory_AfroPal = list("EOO_log_centered", c("Hemisphere3", "migratory1", "American0", "AfroPal1", 
                                                          "Hemisphere3:EOO_log_centered", "Hemisphere3:American0", "Hemisphere3:AfroPal1", 
                                                          "migratory1:Hemisphere3", "migratory1:America0", "migratory1:AfroPal1", "migratory1:EOO_log_centered",
                                                          "EOO_log_centered:American0", "EOO_log_centered:AfroPal1")),
  
  EOO_Hem3_nonMigratory_Asian = list("EOO_log_centered", c("Hemisphere3", "Hemisphere3:EOO_log_centered", "American0", "Asian1",
                                                           "Hemisphere3:American0", "Hemisphere3:Asian1", "EOO_log_centered:American0",
                                                           "EOO_log_centered:Asian1")),
  
  EOO_Hem3_Migratory_Asian = list("EOO_log_centered", c("Hemisphere3", "migratory1", "American0", "Asian1", 
                                                        "Hemisphere3:EOO_log_centered", "Hemisphere3:American0", "Hemisphere3:Asian1", 
                                                        "migratory1:Hemisphere3", "migratory1:America0", "migratory1:Asian1", "migratory1:EOO_log_centered",
                                                        "EOO_log_centered:American0", "EOO_log_centered:Asian1")),
  
  EOO_Hem3_nonMigratory_American = list("EOO_log_centered", c("Hemisphere3", "Hemisphere3:EOO_log_centered")),
  
  EOO_Hem3_Migratory_American = list("EOO_log_centered", c("Hemisphere3", "migratory1", 
                                                           "Hemisphere3:EOO_log_centered", 
                                                           "migratory1:Hemisphere3", "migratory1:EOO_log_centered"))
)

# Initialize a data frame to store the results
effects_summary <- data.frame()

# Loop through each factor and calculate the individual effect sizes
for (factor in names(factors)) {
  main_effect <- factors[[factor]][[1]]
  interaction_terms <- factors[[factor]][[2]]
  
  # Calculate the individual effect for this factor
  individual_effect <- calculate_individual_effect(posterior_samples, main_effect, interaction_terms)
  
  # Summarize the effect (mean and credible intervals)
  mean_effect <- mean(individual_effect)
  lower_ci <- quantile(individual_effect, 0.025)
  upper_ci <- quantile(individual_effect, 0.975)
  
  # Store the results in a data frame
  effects_summary <- rbind(effects_summary, data.frame(
    Factor = factor,
    Mean = mean_effect,
    `Lower 95% CI` = lower_ci,
    `Upper 95% CI` = upper_ci
  ))
}

rownames(effects_summary) <- NULL
print(effects_summary)

# Function to plot the histogram/density plot with CIs for a given effect
plot_effect_distribution <- function(effect_values, effect_name) {
  
  effect_df <- data.frame(Effect = effect_values)
  
  # Calculate the 95% credible intervals
  lower_ci <- quantile(effect_values, 0.025)
  upper_ci <- quantile(effect_values, 0.975)
  
  # Create the plot: Histogram with density curve and CIs
  ggplot(effect_df, aes(x = Effect)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.05, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "blue", size = 1) +
    geom_vline(aes(xintercept = lower_ci), linetype = "dashed", color = "red", size = 1) +
    geom_vline(aes(xintercept = upper_ci), linetype = "dashed", color = "red", size = 1) +
    theme_bw() +
    labs(title = paste("Distribution of the Effect of", effect_name),
         x = "Effect Size",
         y = "Density") +
    annotate("text", x = lower_ci, y = 0.05, label = paste("2.5% CI:", round(lower_ci, 3)), hjust = 1.1, color = "red") +
    annotate("text", x = upper_ci, y = 0.05, label = paste("97.5% CI:", round(upper_ci, 3)), hjust = -0.1, color = "red")
}

dev.off()
# Save the individual effect plots to a PDF
pdf("individualEffect_MCMC2.pdf")

# Iterate through each factor
for (factor in names(factors)) {
  main_effect <- factors[[factor]][[1]]
  interaction_terms <- factors[[factor]][[2]]
  
  individual_effect <- calculate_individual_effect(posterior_samples, main_effect, interaction_terms)
  
  print(plot_effect_distribution(individual_effect, factor))
}

dev.off()
