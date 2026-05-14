# Load the necessary libraries
library(coda)
library(ggplot2)

# Load the posterior samples from the MCMCglmm model
posterior_samples <- readRDS("MCMC_2.rds")  # Assuming you saved the model output
posterior_samples <- as.matrix(posterior_samples$Sol)

# Inspect column names to identify terms in the model
colnames(posterior_samples)

# Adjusted function to calculate individual effects
calculate_individual_effect <- function(samples, factor, interactions = NULL, include_intercept = TRUE) {
  # Start with or without the intercept for each sample
  effect <- if (include_intercept) samples[, "(Intercept)"] else 0
  
  # Safely handle the factor if it is not NULL
  if (!is.null(factor)) {
    if (factor %in% colnames(samples)) {
      effect <- effect + samples[, factor]
    }
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

###############################################################################

# Define the factors and their relevant interaction terms for model 2
factors <- list(
  
  EOO1 = list(
    factor = "EOO_log_centered", 
    interactions = NULL,
    include_intercept = FALSE
  ),
  
  EOO2_migratory = list(
    factor = "EOO_log_centered", 
    interactions = c("migratory1"),
    include_intercept = FALSE
  ),
  
  EOO3_nonMigratory = list(
    factor = "Eoo_log_centered",
    interactions = NULL,
    include_intercept = TRUE
  ),
  
  EOO4_Hem1_nonMigratory = list(
    factor = "EOO_log_centered",
    interactions = c("Hemisphere1", "Hemisphere1:EOO_log_centered"),
    include_intercept = FALSE
  ),
  EOO5_Hem1_Migratory = list(
    factor = "EOO_log_centered",
    interactions = c("Hemisphere1", "Hemisphere1:EOO_log_centered", "migratory1", "migratory1:Hemisphere1",
                     "migratory1:EOO_log_centered"),
    include_intercept = FALSE
  ),
  EOO6_Hem2_nonMigratory = list(
    factor = "EOO_log_centered",
    interactions = NULL,
    include_intercept = TRUE
  ),
  EOO7_Hem2_Migratory = list(
    factor = "EOO_log_centered",
    interactions = c("migratory1"),
    include_intercept = TRUE
  ),
  EOO8_Hem3_nonMigratory = list(
    factor = "EOO_log_centered",
    interactions = c("Hemisphere3", "Hemisphere3:EOO_log_centered"),
    include_intercept = FALSE
  ),
  EOO9_Hem3_Migratory = list(
    factor = "EOO_log_centered",
    interactions = c("Hemisphere3", "Hemisphere3:EOO_log_centered", "migratory1", "migratory1:Hemisphere3",
                     "migratory1:EOO_log_centered"),
    include_intercept = FALSE
  )
)

# Loop through factors and calculate effects
effects <- lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  effect_samples <- calculate_individual_effect(
    samples = posterior_samples,
    factor = f$factor,
    interactions = f$interactions,
    include_intercept = f$include_intercept
  )
  
  # Summarize posterior with mean and 95% CI
  data.frame(
    Mean = mean(effect_samples),
    Lower = quantile(effect_samples, 0.025),
    Upper = quantile(effect_samples, 0.975),
    Category = names(factors)[i]
  )
})

# Combine into a single data frame for summary statistics
summary_data <- do.call(rbind, effects)

# Prepare data for density plotting
density_data <- do.call(rbind, lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  data.frame(
    Effect = calculate_individual_effect(
      samples = posterior_samples,
      factor = f$factor,
      interactions = f$interactions,
      include_intercept = f$include_intercept
    ),
    Category = names(factors)[i]
  )
}))

# Plot posterior distributions with 95% CIs as vertical lines
ggplot() +
  # Add density plots
  geom_density(data = density_data, aes(x = Effect, fill = Category), alpha = 0.6) +
  # Add 95% CIs as vertical lines
  geom_vline(data = summary_data, aes(xintercept = Lower, color = Category), linetype = "dotted", size = 0.8) +
  geom_vline(data = summary_data, aes(xintercept = Upper, color = Category), linetype = "dotted", size = 0.8) +
  theme_minimal() +
  labs(
    title = "Posterior Distributions and 95% Credible Intervals of Effects of EOO",
    x = "Effect Size",
    y = "Density",
    fill = "Category",
    color = "Category"
  ) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(~ Category, scales = "free")

################################################################################

## Effect of henisphere
factors <- list(
  a.Hemisphere1 = list(
    factor = "Hemisphere1",
    interactions = NULL,
    include_intercept = FALSE
  ),
  
  b.Hemisphere2 = list(
    factor = NULL,
    interactions = c("migratory1"),
    include_intercept = TRUE
  ),
  
  c.Hemisphere3 = list(
    factor = "Hemisphere3",
    interactions = NULL,
    include_intercept = FALSE
  ),
  
  Hemisphere1_NonMigratory = list(
    factor = "Hemisphere1",
    interactions = NULL,  # No interaction terms for non-migratory species
    include_intercept = FALSE  # Hemisphere 1 is not the reference category
  ),
  Hemisphere1_Migratory = list(
    factor = "Hemisphere1",
    interactions = c("migratory1", "migratory1:Hemisphere1"),
    include_intercept = FALSE  # Hemisphere 1 is not the reference category
  ),
  Hemisphere2_NonMigratory = list(
    factor = NULL,  # No hemisphere-specific factor since Hemisphere 2 is the reference
    interactions = NULL,  # No interaction terms for the reference hemisphere
    include_intercept = TRUE  # Include intercept for Hemisphere 2 (reference category)
  ),
  Hemisphere2_Migratory = list(
    factor = "migratory1",
    interactions = NULL,  # No hemisphere interaction since Hemisphere 2 is the reference
    include_intercept = TRUE  # Include intercept for Hemisphere 2 (reference category)
  ),
  Hemisphere3_NonMigratory = list(
    factor = "Hemisphere3",
    interactions = NULL,  # No interaction terms for non-migratory species
    include_intercept = FALSE  # Hemisphere 3 is not the reference category
  ),
  Hemisphere3_Migratory = list(
    factor = "Hemisphere3",
    interactions = c("migratory1", "migratory1:Hemisphere3"),
    include_intercept = FALSE  # Hemisphere 3 is not the reference category
  )
)

# Loop through factors and calculate effects
effects <- lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  effect_samples <- calculate_individual_effect(
    samples = posterior_samples,
    factor = f$factor,
    interactions = f$interactions,
    include_intercept = f$include_intercept
  )
  
  # Summarize posterior with mean and 95% CI
  data.frame(
    Mean = mean(effect_samples),
    Lower = quantile(effect_samples, 0.025),
    Upper = quantile(effect_samples, 0.975),
    Category = names(factors)[i]
  )
})

# Combine into a single data frame for summary statistics
summary_data <- do.call(rbind, effects)

# Prepare data for density plotting
density_data <- do.call(rbind, lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  data.frame(
    Effect = calculate_individual_effect(
      samples = posterior_samples,
      factor = f$factor,
      interactions = f$interactions,
      include_intercept = f$include_intercept
    ),
    Category = names(factors)[i]
  )
}))

# Plot posterior distributions with 95% CIs as vertical lines
ggplot() +
  # Add density plots
  geom_density(data = density_data, aes(x = Effect, fill = Category), alpha = 0.6) +
  # Add 95% CIs as vertical lines
  geom_vline(data = summary_data, aes(xintercept = Lower, color = Category), linetype = "dotted", size = 0.8) +
  geom_vline(data = summary_data, aes(xintercept = Upper, color = Category), linetype = "dotted", size = 0.8) +
  theme_minimal() +
  labs(
    title = "Posterior Distributions and 95% Credible Intervals of Effects of Hemisphere",
    x = "Effect Size",
    y = "Density",
    fill = "Category",
    color = "Category"
  ) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(~ Category, scales = "free")

################################################################################

## Effect of migration on decline in different hemispheres
factors <- list(
  
  Migration = list(
    factor = "migratory1",
    interactions = NULL,
    include_intercept = FALSE
  ),
  Migration_Hem1 = list(
    factor = "migratory1",
    interactions = c("Hemisphere1", "migratory1:Hemisphere1"),
    include_intercept = FALSE  # Hemisphere 1 is not the reference category
  ),
  Migration_Hem2 = list(
    factor = "migratory1",
    interactions = NULL,  # No Hemisphere interaction since Hem2 is the reference category
    include_intercept = TRUE  # Include intercept for Hemisphere 2 (reference category)
  ),
  Migration_Hem3 = list(
    factor = "migratory1",
    interactions = c("Hemisphere3", "migratory1:Hemisphere3"),
    include_intercept = FALSE  # Hemisphere 3 is not the reference category
  )
)

# Loop through factors and calculate effects
effects <- lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  effect_samples <- calculate_individual_effect(
    samples = posterior_samples,
    factor = f$factor,
    interactions = f$interactions,
    include_intercept = f$include_intercept
  )
  
  # Summarize posterior with mean and 95% CI
  data.frame(
    Mean = mean(effect_samples),
    Lower = quantile(effect_samples, 0.025),
    Upper = quantile(effect_samples, 0.975),
    Category = names(factors)[i]
  )
})

# Combine into a single data frame for summary statistics
summary_data <- do.call(rbind, effects)

# Prepare data for density plotting
density_data <- do.call(rbind, lapply(seq_along(factors), function(i) {
  f <- factors[[i]]
  data.frame(
    Effect = calculate_individual_effect(
      samples = posterior_samples,
      factor = f$factor,
      interactions = f$interactions,
      include_intercept = f$include_intercept
    ),
    Category = names(factors)[i]
  )
}))

# Plot posterior distributions with 95% CIs as vertical lines
ggplot() +
  # Add density plots
  geom_density(data = density_data, aes(x = Effect, fill = Category), alpha = 0.6) +
  # Add 95% CIs as vertical lines
  geom_vline(data = summary_data, aes(xintercept = Lower, color = Category), linetype = "dotted", size = 0.8) +
  geom_vline(data = summary_data, aes(xintercept = Upper, color = Category), linetype = "dotted", size = 0.8) +
  theme_minimal() +
  labs(
    title = "Posterior Distributions and 95% Credible Intervals of Effects of Migration",
    x = "Effect Size",
    y = "Density",
    fill = "Category",
    color = "Category"
  ) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(~ Category, scales = "free")


################################################################################