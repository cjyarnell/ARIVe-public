# this file performs the analysis on the ARIVe cohort
# including the Bayesian multistate model analysis
# and all of the postprocessing to generate mediation analysis outputs

library(tidyverse)
library(cmdstanr)
library(bayesplot)

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}
ggtext_size <- function(base_size, ratio = 0.8) {
  ratio * base_size / ggplot2::.pt
}

# colours
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

o_dark <- "#DD571C"
o_dark_highlight <- "#FF7417"
o_mid <- "#EF820D"
o_mid_highlight <- "#FDA50F"
o_light <- "#FFCC80"  #"#FFBF00"


c_blue = "cornflower blue"

# CJY working drive:
setwd("C:/Git/ARIVe")
#setwd("~/GitHub/ARIVe")
options(dplyr.summarise.inform = FALSE)

# functions used to complete the analysis, figures, and post-processing
source("ARIVe_stan_piecewise_functions_joint.R")

# load the stan file with / without within-chain threading
fullmod <- cmdstan_model("ARIVe_4state_msm.stan")#,cpp_options = list(stan_threads = TRUE)
                         #)

###################################################
# prepare data for Stan

stanlist <- list()

transitions = list(c(1,2),
                   c(1,3),
                   c(1,4),
                   c(2,3),
                   c(2,4),
                   c(3,4))

for(i in 1:6){
  filename = paste0("data/arive_joint_",
                    transitions[[i]][1],
                    transitions[[i]][2],
                    ".rds")
  stanlist[[i]] <- read_data(readRDS(filename))
  if(transitions[[i]][1] > 1){
    stanlist[[i]] <- stanlist[[i]] %>%
      mutate(duration = duration/28)
  }
}

standat <- make_full_stanlist(stanlist)

#################################################################

# Stan analysis
# better to do this on cluster (takes 45min with 20 threads per chain)
post <- fullmod$sample(data = standat, parallel_chains = 4, iter_warmup = 5,
                       iter_sampling = 5, adapt_delta = 0.8,
#                       threads_per_chain = 20, 
                       refresh = 5)

# save output
post$save_object(file = "post_joint_full.rds")

# reload if doing post-processing from completed analysis
post <- as_cmdstan_fit(files = c("fits/postfull1.csv",
                                 "fits/postfull2.csv",
                                 "fits/postfull3.csv",
                                 "fits/postfull4.csv"))

# this makes the mediation analysis faster
postcsv <- draws_to_csv(post$draws(), dir = getwd(),
                        basename = "post")


#########################################################

# generate hazard ratios and causal mediation estimates

# first, for a reference patient with FiO2 0.80 
# this obtains both hazard ratios (independent of reference patient)
# and mediation outputs (dependent on reference patient characteristics)

gqdat <- reference_covariates() %>% generate_X()
gqdat$N_pred <- nrow(gqdat$X_pred12)
gqdat$N_rep <- 10000

gqmod <- cmdstan_model("models/ARIVe_4state_msm_gq.stan")

post_hrcm <- gqmod$generate_quantities(fitted_params = post,
                                       data = gqdat)

post_hrcm$save_object(file = "post_hrcm.rds")

post_hrcm <- readRDS("fits/post_hrcm.rds")

# probability of less IMV based on hazard ratios (which are independent of reference patient)

apply(post_hrcm$draws("hazard_ratio12")[,,1:3]<1, 3, mean) 

############################################################

# Figure 1

textsize = 12

hr12 %>% 
  mutate(label = paste0(round(mean, 2), " (",
                        round(lb95, 2)," to ", 
                        round(ub95, 2), ")")) %>%
  bind_rows(data.frame(variable = "White",
                       mean = 1, lb95 = 1, ub95 = 1,
                       label = "Reference")) %>%
  mutate(variable = factor(variable, 
                           levels = c("hispanic","black","asian", "White"),
                           labels = c("Hispanic", "Black", "Asian", "White"))) %>%
  ggplot(aes(x = variable, y = mean,
             ymin = lb95, ymax = ub95)) +
  geom_hline(yintercept = 1, color = "grey", alpha = 0.5,
             size = 1) +
  geom_pointrange(size = 0.8) +
  coord_flip() +
  theme_minimal(base_size = textsize) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.6, 0.8, 1),
                     limits = c(0.6, 1.25)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_text(aes(label = label,
                y = 1.03), hjust = 0,
            size = ggtext_size(textsize)) +
  labs(x = "", 
       y = "Posterior hazard ratio (95% credible interval)")#,
       #title = "Hazard ratio for rate of invasive ventilation")

ggsave("figures/figure1.svg", width = 6, height = 4, dpi = 300)

# Hazard ratios - big forest plots for other transitions

output_hazard_ratio(post_hrcm,
                    "hazard_ratio12",
                    names(stanlist[[1]])[-c(1:3)]) %>%
  forest1()+ 
  labs(title = "Hazard ratios for transition from oxygen therapy to invasive ventilation")

ggsave("forest12.svg", height = 10, width = 8)

output_hazard_ratio(post_hrcm,
                    "hazard_ratio13",
                    names(stanlist[[2]])[-c(1:3)]) %>%
  forest1() + 
  labs(title = "Hazard ratios for transition from oxygen therapy to ICU discharge")
ggsave("forest13.svg", height = 10, width = 8)

output_hazard_ratio(post_hrcm,
                    "hazard_ratio14",
                    names(stanlist[[3]])[-c(1:3)]) %>%
  forest1() + 
  labs(title = "Hazard ratios for transition from oxygen therapy to death")
ggsave("forest14.svg", height = 10, width = 8)

output_hazard_ratio(post_hrcm,
                    "hazard_ratio23",
                    names(stanlist[[4]])[-c(1:3)]) %>%
  forest2() + 
  labs(title = "Hazard ratios for transition from invasive ventilation to ICU discharge")
ggsave("forest23.svg", height = 8, width = 8)

output_hazard_ratio(post_hrcm,
                    "hazard_ratio24",
                    names(stanlist[[5]])[-c(1:3)]) %>%
  forest2() + 
  labs(title = "Hazard ratios for transition from invasive ventilation to death")
ggsave("forest24.svg", height = 8, width = 8)

output_hazard_ratio(post_hrcm,
                    "hazard_ratio34",
                    names(stanlist[[6]])[-c(1:3)]) %>%
  forest3() + 
  labs(title = "Hazard ratios for transition from ICU discharge to death")
ggsave("forest34.svg", height = 8, width = 8)


# weibull shapes

alpha <- post_hrcm$draws(variables = c("alpha"))
data.frame(variable = dimnames(alpha)$variable, 
           mean = apply(alpha, 3, mean),
           lb95 = apply(alpha, 3, quantile, 0.025),
           ub95 = apply(alpha, 3, quantile, 0.975))

######################################################################

# average patient causal mediation

# timeweighted mean by patient
ptlist <- 
  timeweighted_mean("data/arive_joint_12.rds")

# time-weighted mean of time-varying variables table
tw_means <- ptlist %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  pivot_longer(cols = c(white, black, hispanic, asian), names_to = "race") %>%
  filter(value == 1) %>%
  group_by(race) %>%
  summarise(across(heart_rate:pressor, mean))
tw_75 <- ptlist %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  pivot_longer(cols = c(white, black, hispanic, asian), names_to = "race") %>%
  filter(value == 1) %>%
  group_by(race) %>%
  summarise(across(heart_rate:pressor, quantile, 0.975))
tw_25 <- ptlist %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  pivot_longer(cols = c(white, black, hispanic, asian), names_to = "race") %>%
  filter(value == 1) %>%
  group_by(race) %>%
  summarise(across(heart_rate:pressor, quantile, 0.025))

bind_rows(tw_means, tw_75, tw_25, .id = "measure") %>%
  mutate(measure = factor(measure,
                          levels = 1:3,
                          labels = c("mean","ub95","lb95"))) %>%
  pivot_longer(heart_rate:pressor, names_to = "variable") %>%
  left_join(readRDS("CreateCohort/joint_mean_sd.rds"),
            by = "variable") %>%
  mutate(value = ifelse(variable != "pressor",
                        value*sd + mean, value)) %>%
  mutate(value = case_when(
    variable %in% c("heart_rate", "resp_rate") ~ exp(value),
    variable == "spo2" ~ round(inv_logit(value)*101),
    variable == "fio2" ~ round(inv_logit(value)*79.3 + 20.9),
    variable == "gcs"  ~ round(inv_logit(value)*12.2 + 2.9),
    variable == "pressor" ~ value
  )) %>%
  select(measure:value) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(label = paste0(round(mean,2), " (",
                        round(lb95,2), " to ",
                        round(ub95,2), ")")) %>%
  select(race, variable, label) %>%
  pivot_wider(names_from = race, values_from = label) %>%
  relocate(variable, white, black, hispanic, asian) %>%
  write_csv(file = "time-weighted-means.csv")

# time-weighted mean of time-varying variables table
tw_means <- ptlist %>%
  summarise(across(heart_rate:pressor, mean))
tw_75 <- ptlist %>%
  summarise(across(heart_rate:pressor, quantile, 0.975))
tw_25 <- ptlist %>%
  summarise(across(heart_rate:pressor, quantile, 0.025))

bind_rows(tw_means, tw_75, tw_25, .id = "measure") %>%
  mutate(measure = factor(measure,
                          levels = 1:3,
                          labels = c("mean","ub95","lb95"))) %>%
  pivot_longer(heart_rate:pressor, names_to = "variable") %>%
  left_join(readRDS("CreateCohort/joint_mean_sd.rds"),
            by = "variable") %>%
  mutate(value = ifelse(variable != "pressor",
                        value*sd + mean, value)) %>%
  mutate(value = case_when(
    variable %in% c("heart_rate", "resp_rate") ~ exp(value),
    variable == "spo2" ~ round(inv_logit(value)*101),
    variable == "fio2" ~ round(inv_logit(value)*79.3 + 20.9),
    variable == "gcs"  ~ round(inv_logit(value)*12.2 + 2.9),
    variable == "pressor" ~ value
  )) %>%
  select(measure:value) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(label = paste0(round(mean,2), " (",
                        round(lb95,2), " to ",
                        round(ub95,2), ")")) %>%
  select(variable, label) %>%
  write_csv(file = "time-weighted-means_wc.csv")

readRDS("data/arive_joint_12.rds") %>%
  summarise(imv = sum(status == 1)) %>%
  ungroup() %>%
  summarise(sum(imv),
            mean(imv))

readRDS("data/arive_joint_12.rds") %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  select(status, white, black, hispanic, asian) %>%
  pivot_longer(white:asian, names_to = "race") %>%
  filter(value == 1) %>%
  group_by(id, race) %>%
  summarise(imv = sum(status == 1)) %>%
  ungroup() %>%
  group_by(race) %>%
  summarise(sum(imv),
            mean(imv))

##############################################

# Causal mediation for the average patients

# average Black patient

bptlist <- ptlist %>% filter(black == 1) %>%
  select(-id, -duration) %>%
  ungroup() %>%
  summarise_all(mean)

sb <- generate_X(bptlist)
sb$N_pred <- nrow(sb$X_pred12)
sb$N_rep <- 10000

postb <- gqmod$generate_quantities(fitted_params = postcsv,
                                   data = sb,
                                   parallel_chains = 4)



# average Asian patient

aptlist <- ptlist %>% filter(asian == 1) %>%
  select(-id, -duration) %>%
  ungroup() %>%
  summarise_all(mean)

sa <- generate_X(aptlist)
sa$N_pred <- nrow(sa$X_pred12)
sa$N_rep <- 10000

posta <- gqmod$generate_quantities(fitted_params = postcsv,
                                   data = sa,
                                   parallel_chains = 4)

# Average Hispanic patient

hptlist <- ptlist %>% filter(hispanic == 1) %>%
  select(-id, -duration) %>%
  ungroup() %>%
  summarise_all(mean)

sh <- generate_X(hptlist)
sh$N_pred <- nrow(sh$X_pred12)
sh$N_rep <- 10000

posth <- gqmod$generate_quantities(fitted_params = postcsv,
                                   data = sh,
                                   parallel_chains = 4)

package_gq(posth)

# Average white patient

wptlist <- ptlist %>% filter(black == 0, asian == 0, hispanic == 0) %>%
  select(-id, -duration) %>%
  ungroup() %>%
  summarise_all(mean) 

sw <- generate_X(wptlist)
sw$N_pred <- nrow(sw$X_pred12)
sw$N_rep <- 10000

postw <- gqmod$generate_quantities(fitted_params = postcsv,
                                   data = sw,
                                   parallel_chains = 4)

package_gq(postw)


# average patient characteristics table

mean_sd <- readRDS("CreateCohort/mean_sd.rds")

mean_age = 66.6
sd_age = 15.46

bind_rows(white = wptlist,
          black = bptlist,
          hispanic = hptlist,
          asian = aptlist,
          .id = "race") %>%
  select(-(intercept:hispanic)) %>%
  mutate(age = age*sd_age+mean_age,
         heart_rate = exp(heart_rate*filter(mean_sd, variable == "heart_rate")$sd +
                            filter(mean_sd, variable == "heart_rate")$mean),
         resp_rate = exp(resp_rate*filter(mean_sd, variable == "resp_rate")$sd +
                           filter(mean_sd, variable == "resp_rate")$mean),
         spo2 = inv_logit(spo2*filter(mean_sd, variable == "spo2")$sd +
                                     filter(mean_sd, variable == "spo2")$mean)*101,
         fio2 = inv_logit(fio2*filter(mean_sd, variable == "fio2")$sd +
                            filter(mean_sd, variable == "fio2")$mean)*79.3+20.9,
         gcs = inv_logit(gcs*filter(mean_sd, variable == "GCS")$sd +
                           filter(mean_sd, variable == "GCS")$mean)*12.2 + 2.9) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  write_csv(file = "average_pts.csv")

# Subgroups by sex

gqdat_f <- 
  ptlist %>% filter(gender == 0) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_f <- outcomes_by_race(gqdat_f)

saveRDS(outcomes_f, "fits/outcomes_f.rds")

gqdat_m <- ptlist %>% filter(gender == 1) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_m <- outcomes_by_race(gqdat_m)

saveRDS(outcomes_m, "fits/outcomes_m.rds")

# Subgroups by age

agebreaks <- (c(50, 65, 80)-66.6)/15.46

gqdat_a1 <- ptlist %>% filter(age < agebreaks[1]) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_a1 <- outcomes_by_race(gqdat_a1)
saveRDS(outcomes_a1, "fits/outcomes_a1.rds")

gqdat_a2 <- ptlist %>% filter(age < agebreaks[2], 
                              age >= agebreaks[1]) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_a2 <- outcomes_by_race(gqdat_a2)
saveRDS(outcomes_a2, "fits/outcomes_a2.rds")

gqdat_a3 <- ptlist %>% filter(age < agebreaks[3],
                              age >= agebreaks[2]) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_a3 <- outcomes_by_race(gqdat_a3)
saveRDS(outcomes_a3, "fits/outcomes_a3.rds")

gqdat_a4 <- ptlist %>% filter(age >= agebreaks[3]) %>%
  select(-id, -duration) %>%
  mutate(race = case_when(black == 1 ~ "black",
                          asian == 1 ~ "asian",
                          hispanic == 1~  "hispanic",
                          black == 0 ~ "white")) %>%
  group_by(race) %>%
  summarise_all(mean)

outcomes_a4 <- outcomes_by_race(gqdat_a4)
saveRDS(outcomes_a4, "fits/outcomes_a4.rds")


bind_rows(tableready(readRDS("fits/outcomes_f.rds")),
          tableready(readRDS("fits/outcomes_m.rds")),
          tableready(readRDS("fits/outcomes_a1.rds")),
          tableready(readRDS("fits/outcomes_a2.rds")),
          tableready(readRDS("fits/outcomes_a3.rds")),
          tableready(readRDS("fits/outcomes_a4.rds"))) %>%
  write_csv("subgroups_tbl2.csv")

# average patient outcomes

w <- package_gq(readRDS("fits/postw.rds"))
b <- package_gq(readRDS("fits/postb.rds"))
a <- package_gq(readRDS("fits/posta.rds"))
h <- package_gq(readRDS("fits/posth.rds"))


bind_rows(bind_rows(w,.id = "outcome"),
          bind_rows(b,.id = "outcome"),
          bind_rows(h,.id = "outcome"),
          bind_rows(a,.id = "outcome"),
          .id = "race") %>%
  mutate(race = factor(race,
                       levels = 1:4,
                       labels = c("White",
                                  "Black",
                                  "Hispanic",
                                  "Asian")),
         label = paste0(
           round(100*mean,1), " (",
           round(100*lb95, 1), " to ",
           round(100*ub95, 1), ")")
         ) %>%
  filter(race == subgroup) %>%
  select(outcome, race, label) %>%
  pivot_wider(names_from = race,
              values_from = label) %>% 
  filter(outcome != "ttimv") %>%
  write.csv("table2_average_patient.csv")

# reference patient part of table 2

bind_rows(package_gq(post_hrcm), .id = "outcome") %>% 
  mutate(label = paste0(
    round(100*mean,1), " (",
    round(100*lb95, 1), " to ",
    round(100*ub95, 1), ")")
  ) %>%
  select(outcome, subgroup, label) %>%
  pivot_wider(names_from = subgroup,
              values_from = label) %>%
  write.csv("tbl2_reference_patient80.csv")

bind_rows(package_gq(readRDS("fits/post_hrcm100.rds")), .id = "outcome") %>% 
  mutate(label = paste0(
    round(100*mean,1), " (",
    round(100*lb95, 1), " to ",
    round(100*ub95, 1), ")")
  ) %>%
  select(outcome, subgroup, label) %>%
  pivot_wider(names_from = subgroup,
              values_from = label) %>%
  write.csv("tbl2_reference_patient100.csv")

bind_rows(package_gq(readRDS("fits/post_hrcm50.rds")), .id = "outcome") %>% 
  mutate(label = paste0(
    round(100*mean,1), " (",
    round(100*lb95, 1), " to ",
    round(100*ub95, 1), ")")
  ) %>%
  select(outcome, subgroup, label) %>%
  pivot_wider(names_from = subgroup,
              values_from = label) %>%
  write.csv("tbl2_reference_patient50.csv")

package_gq(readRDS("fits/post_hrcm50.rds"))

# reference patient by fio2

gq50 <- reference_covariates() %>% mutate(fio2 = 0.68) %>% 
  single_patient_gq("post_hrcm50.rds")

gq60 <- reference_covariates() %>% mutate(fio2 = 0.99) %>% 
  single_patient_gq("post_hrcm60.rds")

gq70 <- reference_covariates() %>% mutate(fio2 = 1.305) %>% 
  single_patient_gq("post_hrcm70.rds")

# 0.80 already done

gq90 <- reference_covariates() %>% mutate(fio2 = 2.17) %>% 
  single_patient_gq("post_hrcm90.rds")

gq100 <- reference_covariates() %>% mutate(fio2 = 4.65) %>% 
  single_patient_gq("post_hrcm100.rds")

bind_gq_fio2 <- function(name,x){
  bind_cols(data.frame(fio2 = x),
            bind_rows(package_gq(readRDS(name))[1:4],
            .id = "variable"))
}

refptdf <- bind_rows(
  bind_gq_fio2("fits/post_hrcm50.rds", 50),
  bind_gq_fio2("fits/post_hrcm60.rds", 60),
  bind_gq_fio2("fits/post_hrcm70.rds", 70),
  bind_gq_fio2("fits/post_hrcm.rds", 80),
  bind_gq_fio2("fits/post_hrcm90.rds", 90),
  bind_gq_fio2("fits/post_hrcm100.rds", 100),
) %>%
  bind_rows(data.frame(fio2 = seq(from = 50, to = 100, by = 10),
                       subgroup = "White",
                       variable = "ie",
                       mean = 0,
                       lb95 = 0,
                       ub95 = 0)) %>%
  mutate(fio2 = fio2/100,
         subgroup = factor(subgroup,
                           levels = c("White","Asian",
                                      "Black","Hispanic"))) %>%
  mutate(mean = ifelse(variable == "ttimv", mean, 100*mean),
         ub95 = ifelse(variable == "ttimv", ub95, 100*ub95),
         lb95 = ifelse(variable == "ttimv", lb95, 100*lb95)) %>%
  mutate(mean = ifelse(variable == "ie", -mean, mean),
         ub95 = ifelse(variable == "ie", -ub95, ub95),
         lb95 = ifelse(variable == "ie", -lb95, lb95)) %>%
  mutate(variable = factor(variable,
                           levels = c("imv","ttimv","te","ie"),
                           labels = c(
                             "Invasive ventilation by day 28 (%)",
                             "Average time to invasive ventilation (days)",
                             "Survival to day 28 (%)",
                             paste(sep = "\n",
                             "Absolute change in 28-day survival mediated by differences in invasive ventilation rate"
                             )
                             )))

# figure 2

refptdf %>%
  filter(variable != "Average time to invasive ventilation (days)") %>%
  ggplot(aes(y = mean, 
             x = fio2,
             ymin = lb95, 
             ymax = ub95, 
             color = subgroup)) +
  geom_pointrange(position = position_dodge(width = 0.05),
                  alpha = 0.6,
                  size = 0.4) +
  facet_wrap(.~variable,
             scales = "free",
             ncol = 1) +
  scale_color_manual(values = c("grey50",c_blue, c_dark, "black"),
                     name = "Race/Ethnicity") +
  theme_minimal() +
  scale_x_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9,1)) +
  labs(y = "", x = "Inspired oxygen fraction",
       title = "Outcomes for a reference patient on high-flow nasal cannula") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(2, "line"))
  ggsave(filename = "figures/reference_patient_outcomes.svg", dpi = 300,
         height = 10, width = 7)

  
  refptdf %>%
    filter(!(variable %in% c(
      "Invasive ventilation by day 28 (%)",
      "Average time to invasive ventilation (days)",
      "Survival to day 28 (%)"))) %>%
    ggplot(aes(y = mean, 
               x = fio2,
               ymin = lb95, 
               ymax = ub95, 
               color = subgroup)) +
    geom_pointrange(position = position_dodge(width = 5),
                    alpha = 0.6,
                    size = 1) +
    scale_color_manual(values = c("grey50",c_blue, c_dark, "black"),
                       name = "Race/Ethnicity") +
    theme_minimal() +
    scale_x_continuous(breaks = c(50,60,70,80,90,100)) +
    labs(y = "Survival change (absolute, %)", x = "Inspired oxygen fraction (%)") +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(2, "line"))
  ggsave(filename = "figures/reference_patient_outcomes_tghriposter.svg", dpi = 300,
         height = 3, width = 8)

  # sensitivity analysis plots
  
  postmimic <- readRDS("fits/mimic_O2toIMV.rds")
  mimic_12 <- readRDS("data/arive_12.rds") %>%
    mutate(intercept = 1,
           duration = stop-start) %>%
    ungroup() %>%
    select(-stay_id, -start, -stop) %>%
    relocate(duration, status, intercept)
  
  mimic_draws <- postmimic$draws(variables = "beta")
  
  df <- data.frame(variable = names(mimic_12)[c(-1,-2)],
                   mean = apply(mimic_draws, 3, mean),
                   lb95 = apply(mimic_draws, 3, quantile, 0.025),
                   ub95 = apply(mimic_draws, 3, quantile, 0.975)) %>%
    mutate(mean = exp(mean),
           lb95 = exp(lb95),
           ub95 = exp(ub95))

  
  df %>%
    filter(variable != "intercept") %>%
    forest_mimic()
  
  posteicu <- readRDS("fits/eicu_sens_samples.rds")
  eicu_12 <- readRDS("data/eicu_sens_dat.rds")
  
  eicu_draws <- posteicu$draws(variables = "beta")[,,-1]
  
  df2 <- data.frame(variable = names(eicu_12$X_cen)[-1],
                   mean = apply(eicu_draws, 3, mean),
                   lb95 = apply(eicu_draws, 3, quantile, 0.025),
                   ub95 = apply(eicu_draws, 3, quantile, 0.975)) %>%
    mutate(mean = exp(mean),
           lb95 = exp(lb95),
           ub95 = exp(ub95))

  df2 %>%
    forest_eicu() +
    labs(title = "eICU: Hazard ratios for invasive ventilation")
  ggsave("figures/eicu_forest.svg", height = 10, width = 8)

# simple eicu
  
  posteicu_s <- readRDS("fits/eicu_simple_sens_samples.rds")
  names <- read_data(readRDS("data/arive_joint_12.rds")) %>%
    filter(eicu == 1) %>%
    select(-eicu,-G1, -G2, -G3)
  
  eicu_draws_s <- posteicu_s$draws(variables = "beta")[,,-1]
  
  dfs <- data.frame(variable = names(names)[-c(1,2,3)],
                    mean = apply(eicu_draws_s, 3, mean),
                    lb95 = apply(eicu_draws_s, 3, quantile, 0.025),
                    ub95 = apply(eicu_draws_s, 3, quantile, 0.975)) %>%
    mutate(mean = exp(mean),
           lb95 = exp(lb95),
           ub95 = exp(ub95))
  
  dfs %>% write_csv("eicu_simple.csv")
  
  
  dfs %>%
    forest_eicu2() +
    labs(title = "eICU: Hazard ratios for invasive ventilation")
  ggsave("figures/eicu_forest_simple.svg", height = 10, width = 8)
  
  
##   Pairs plots
  
post12 <- post$draws(variables = "beta12", format = "draws_matrix")
dimnames(post12)$variable <- names(stanlist[[1]])[-c(1:2)]   

cormat <- cor(post12)
library(corrplot)
corrplot(cormat, 
         tl.col = "grey50")
ggsave(filename = "figures/parameters_correlationplot.svg",
       width = 8, height = 8)


# all fio2 rows with fio2 in particular brackets


fio2_breaks <- c(45, 55, 65, 75, 85, 95) 
fio2_mean <- -1.54
fio2_sd <- 1.53
fio2_breaks_transformed <- (logit((fio2_breaks-20.9)/79.3)-fio2_mean)/fio2_sd


# 45-55
subset50 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[1],
                   fio2 <  fio2_breaks_transformed[2]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()

o50 <- outcomes_by_race(subset50)

  
# 55-65
subset60 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[2],
                   fio2 <  fio2_breaks_transformed[3]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()


# 65-75
subset70 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[3],
                   fio2 <  fio2_breaks_transformed[4]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()

# 75-85
subset80 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[3],
                   fio2 <  fio2_breaks_transformed[4]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()
  

# 85-95
subset90 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[4],
                   fio2 <  fio2_breaks_transformed[5]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()

# 95-100
subset100 <- filter(stanlist[[1]],
                   fio2 >= fio2_breaks_transformed[5]) %>%
  timeweighted_mean_fromdf() %>%
  select(-duration) %>%
  ready_for_gq()

o100 <- outcomes_by_race(subset100)
