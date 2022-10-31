# sensitivity analyses MIMIC and eICU individually

library(tidyr)
library(dplyr)
library(arrow)
library(ggpubr)
library(cmdstanr)
# CJY working drive:
#setwd("C:/Git/ARIVe")
#setwd("~/GitHub/ARIVe")
options(dplyr.summarise.inform = FALSE)

source("ARIVe_stan_piecewise_functions_joint.R")

# mimic with insurance and WOB covariate

mimic_12 <- read_data(readRDS("data/mimic_12.rds"))

standat <- list(P = ncol(mimic_12) - 2,
                N_cen = sum(1-mimic_12$status),
                N_obs = sum(mimic_12$status),
                X_cen = select(filter(mimic_12, status == 0), -status, -duration),
                X_obs = select(filter(mimic_12, status == 1), -status, -duration),
                y_cen = mimic_12$duration[mimic_12$status == 0],
                y_obs = mimic_12$duration[mimic_12$status == 1]
)

library(cmdstanr)

mod <- cmdstan_model("models/ARIVe_O2IMV_msm.stan")

samples1 <- mod$sample(data = standat,
                      iter_warmup = 250, 
                      iter_sampling = 250,
                      chains = 4,
                      parallel_chains = 4)

samples1$save_object(file = "mimic_sens_samples.rds")

samples1 <- readRDS("mimic_sens_samples.rds")

draws <- samples1$draws(variable = "beta")

hrdf <- data.frame(variable = names(mimic_12)[c(-1,-2)],
                   mean = apply(draws, 3, mean),
                   lb95 = apply(draws, 3, quantile, 0.025),
                   ub95 = apply(draws, 3, quantile, 0.975)) %>%
  mutate(mean = exp(mean),
         lb95 = exp(lb95),
         ub95 = exp(ub95))

hrdf %>%
  filter(variable != "intercept") %>%
  forest_mimic() +
  labs(title = "MIMIC-IV: Hazard ratios for invasive ventilation")
ggsave("figures/mimic_forest.svg", height = 10, width = 8)


# eICU with hospital ICU proportion non-white patients covariate

standat <- readRDS("data/eicu_sens_dat.rds")

library(cmdstanr)

samples <- mod$sample(data = standat,
                      iter_warmup = 250, 
                      iter_sampling = 250,
                      chains = 4,
                      parallel_chains = 4)

samples$save_object(file = "eicu_sens_samples.rds")

samples <- readRDS("eicu_sens_samples.rds")

draws <- samples$draws(variable = "beta")

hrdf <- data.frame(variable = names(esens)[c(-1,-2)],
                 mean = apply(draws, 3, mean),
                 lb95 = apply(draws, 3, quantile, 0.025),
                 ub95 = apply(draws, 3, quantile, 0.975)) %>%
  mutate(mean = exp(mean),
         lb95 = exp(lb95),
         ub95 = exp(ub95))

hrdf %>%
  filter(variable != "intercept") %>%
  forest_eicu() +
  labs(title = "eICU: Hazard ratios for invasive ventilation")
ggsave("figures/eicu_forest.svg", height = 10, width = 8)

# eicu simple (same covariates as primary analysis, minus year of admission and database)

eicu_s <- read_data(readRDS("data/arive_joint_12.rds")) %>%
    filter(eicu == 1) %>%
    select(-eicu, -G1, -G2, -G3)

standat <- list(P = ncol(eicu_s) - 2,
                N_cen = sum(1-eicu_s$status),
                N_obs = sum(eicu_s$status),
                X_cen = select(filter(eicu_s, status == 0), -status, -duration),
                X_obs = select(filter(eicu_s, status == 1), -status, -duration),
                y_cen = eicu_s$duration[eicu_s$status == 0],
                y_obs = eicu_s$duration[eicu_s$status == 1]
)

mod <- cmdstan_model("models/ARIVe_O2IMV_msm.stan")

samples_es <- mod$sample(data = standat,
                      iter_warmup = 250, 
                      iter_sampling = 250,
                      chains = 4,
                      parallel_chains = 4)

samples_es$save_object(file = "eicu_simple_sens_samples.rds")

# mimic simple (no insurance status, no work of breathing, no region, no teaching status)


mimic_s <- read_data(readRDS("data/arive_joint_12.rds")) %>%
    filter(eicu == 0) %>%
    select(-(teachingstatus:eicu))

standat <- list(P = ncol(mimic_s) - 2,
                N_cen = sum(1-mimic_s$status),
                N_obs = sum(mimic_s$status),
                X_cen = select(filter(mimic_s, status == 0), -status, -duration),
                X_obs = select(filter(mimic_s, status == 1), -status, -duration),
                y_cen = mimic_s$duration[mimic_s$status == 0],
                y_obs = mimic_s$duration[mimic_s$status == 1]
)

mod <- cmdstan_model("models/ARIVe_O2IMV_msm.stan")

samples_ms <- mod$sample(data = standat,
                      iter_warmup = 250, 
                      iter_sampling = 250,
                      chains = 4,
                      parallel_chains = 4)

samples_ms$save_object(file = "mimic_simple_sens_samples.rds")
