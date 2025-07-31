library(coda)
library(ggplot2)

# Load the posterior samples from the MCMCglmm model
posterior_samples <- readRDS("MCMC_1.rds")
posterior_samples <- as.matrix(posterior_samples$Sol)

colnames(posterior_samples)

# Function to calculate the individual effect for a given factor, accounting for interactions
calculate_individual_effect <- function(samples, factor, interactions = NULL) {
  #browser()
  # Intercept for each sample
  effect <- samples[, "(Intercept)"]
  
  # Add the main effect for the factor
  if (factor %in% colnames(samples)) {
    effect <- effect + samples[, factor]
  }
  
  # Add any interaction terms involving this factor
  if (!is.null(interactions)) {
    for (interaction in interactions) {
      if (interaction %in% colnames(samples)) {
        effect <- effect + samples[, interaction]
      }
    }
  }
  
  return(effect)
}

# Define the factors and their relevant interaction terms
factors <- list(
  migratory1 = list("migratory1", c("migratory1:American1", "migratory1:AfroPal1", "migratory1:Asian1",
                                    "migratory1:Hemisphere1", "migratory1:Hemisphere3", "migratory1:EOO_log_centered", 
                                    "migratory1:Artificial1", "migratory1:Forest1", "migratory1:Grassland1", 
                                    "migratory1:Savanna1", "migratory1:Shrubland1", "migratory1:Wetlands1")),
  
  Hemisphere1 = list("Hemisphere1", c("Hemisphere1:American1", "Hemisphere1:AfroPal1", "Hemisphere1:Asian1", 
                                      "Hemisphere1:Artificial1", "Hemisphere1:Forest1", "Hemisphere1:Grassland1", 
                                      "Hemisphere1:Savanna1", "Hemisphere1:Shrubland1", "Hemisphere1:Wetlands1", 
                                      "Hemisphere1:EOO_log_centered")),
  
  Hemisphere3 = list("Hemisphere3", c("Hemisphere3:American1", "Hemisphere3:AfroPal1", "Hemisphere3:Asian1", 
                                      "Hemisphere3:Artificial1", "Hemisphere3:Forest1", "Hemisphere3:Grassland1", 
                                      "Hemisphere3:Savanna1", "Hemisphere3:Shrubland1", "Hemisphere3:Wetlands1", 
                                      "Hemisphere3:EOO_log_centered")),
  
  American1 = list("American1", c("American1:Hemisphere1", "American1:Hemisphere3", "American1:EOO_log_centered", 
                                  "American1:Artificial1", "American1:Forest1", "American1:Grassland1", 
                                  "American1:Savanna1", "American1:Shrubland1", "American1:Wetlands1")),
  
  AfroPal1 = list("AfroPal1", c("AfroPal1:Hemisphere1", "AfroPal1:Hemisphere3", "AfroPal1:EOO_log_centered", 
                                "AfroPal1:Artificial1", "AfroPal1:Forest1", "AfroPal1:Grassland1", 
                                "AfroPal1:Savanna1", "AfroPal1:Shrubland1", "AfroPal1:Wetlands1")),
  
  Asian1 = list("Asian1", c("Asian1:Hemisphere1", "Asian1:Hemisphere3", "Asian1:EOO_log_centered", 
                            "Asian1:Artificial1", "Asian1:Forest1", "Asian1:Grassland1", 
                            "Asian1:Savanna1", "Asian1:Shrubland1", "Asian1:Wetlands1")),
  
  EOO_log_centered = list("EOO_log_centered", c("migratory1:EOO_log_centered", "American1:EOO_log_centered", 
                                                "AfroPal1:EOO_log_centered", "Asian1:EOO_log_centered", 
                                                "Hemisphere1:EOO_log_centered", "Hemisphere3:EOO_log_centered", 
                                                "Artificial1:EOO_log_centered", "Forest1:EOO_log_centered", 
                                                "Grassland1:EOO_log_centered", "Savanna1:EOO_log_centered", 
                                                "Shrubland1:EOO_log_centered", "Wetlands1:EOO_log_centered")),
  
  Artificial1 = list("Artificial1", c("Artificial1:Hemisphere1", "Artificial1:Hemisphere3", "Artificial1:American1", 
                                      "Artificial1:AfroPal1", "Artificial1:Asian1", "Artificial1:EOO_log_centered")),
  
  Forest1 = list("Forest1", c("Forest1:Hemisphere1", "Forest1:Hemisphere3", "Forest1:American1", 
                              "Forest1:AfroPal1", "Forest1:Asian1", "Forest1:EOO_log_centered")),
  
  Grassland1 = list("Grassland1", c("Grassland1:Hemisphere1", "Grassland1:Hemisphere3", "Grassland1:American1", 
                                    "Grassland1:AfroPal1", "Grassland1:Asian1", "Grassland1:EOO_log_centered")),
  
  Savanna1 = list("Savanna1", c("Savanna1:Hemisphere1", "Savanna1:Hemisphere3", "Savanna1:American1", 
                                "Savanna1:AfroPal1", "Savanna1:Asian1", "Savanna1:EOO_log_centered")),
  
  Shrubland1 = list("Shrubland1", c("Shrubland1:Hemisphere1", "Shrubland1:Hemisphere3", "Shrubland1:American1", 
                                    "Shrubland1:AfroPal1", "Shrubland1:Asian1", "Shrubland1:EOO_log_centered")),
  
  Wetlands1 = list("Wetlands1", c("Wetlands1:Hemisphere1", "Wetlands1:Hemisphere3", "Wetlands1:American1", 
                                  "Wetlands1:AfroPal1", "Wetlands1:Asian1", "Wetlands1:EOO_log_centered"))
)



factors <- list(
  
  EOO_Hem1_nonMigratory = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered")),
  
  EOO_Hem1_Migratory = list("EOO_log_centered", c("Hemisphere1", "Hemisphere1:EOO_log_centered", "migratory1",
                                                  "migratory1:EOO_log_centered", "migratory1:Hemisphere1"))
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

pdf("individualEffect_subset.pdf")

# Iterate through each factor
for (factor in names(factors)) {
  main_effect <- factors[[factor]][[1]]
  interaction_terms <- factors[[factor]][[2]]
  
  individual_effect <- calculate_individual_effect(posterior_samples, main_effect, interaction_terms)
  
  print(plot_effect_distribution(individual_effect, factor))
}

dev.off()

