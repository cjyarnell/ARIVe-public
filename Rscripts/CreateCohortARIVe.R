###################################################
#
# ARIVe Cohort
# for the stricter NIV | HFNC & FiO2 > 40 & sat <= 97
# builds cohort database from 2 parquet files:
# PAHRC_Timevarying
# PAHRC_Baseline
# builds a 96h cohort


# First, MIMIC

library(tidyverse)
library(arrow)

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

# CJY working drive:
setwd("C:/Git/ARIVe")
# Fred's working drive:
#setwd("C:/Users/freda/Documents/GitHub/ARIVe")
ariv <- read_parquet(
  file = "CreateCohort/ARIVe_MIMIC_Timevarying"
) 

ariv_b <- read_parquet(
  file = "CreateCohort/ARIVe_MIMIC_Baseline"
) %>%
  filter(ethnicity %in% c("ASIAN",
                          "BLACK/AFRICAN AMERICAN",
                          "WHITE",
                          "HISPANIC/LATINO"))

# filter to patients with known ethnicity
ariv <- filter(ariv,
               stay_id %in% ariv_b$stay_id)

# combine GCS components into single score
gcs <- filter(ariv,
              grepl("GCS", variable)) %>%
  group_by(stay_id, time) %>%
  summarise(variable = "GCS",
            value = sum(value)) %>%
  distinct()

# add back to dataframe

ariv <- ariv %>%
  bind_rows(gcs) %>%
  filter(!(variable %in% c("GCS_eyes", "GCS_motor", "GCS_verbal"))) %>%
  arrange(stay_id, time, variable, value)


# suggested transformation function 
# use log for all continuous variables 
# except for spo2, fio2, GCS
# for those three, use the below

transform <- function(row, log_group){
  variable <- row[3]
  value <- as.numeric(row[4])
  if(variable %in% log_group){log(value)}
  # use logit transformation functions for GCS, spo2, and fio2
  # because they have a restricted range 
  else if (variable == "GCS"){logit((value-2.9)/12.2)}
  else if (variable == "spo2"){logit(value/101)} # take the floor when transforming back
  else if (variable == "fio2"){logit((value-20.9)/79.3)} # take the floor when transforming back
  else value
}

log_group <- c("heart_rate", "resp_rate")

ariv_transformed <- ariv

ariv_transformed$value <- apply(ariv, 1, transform, log_group)

# how many NaN introduced by log?
sum(is.na(ariv_transformed$value)) # 140, mostly GCS < 3 which is impossible anyway
length(ariv_transformed$value) # > 4.2 million, so filter them out
ariv_transformed <- filter(ariv_transformed, 
                           !is.na(value))


# then record mean and standard deviation of each continuous variable
# then rescale each continuous variable by its mean and standard deviation
# code for standardization is now here, adjusted from your perfectly reasonable version


# rescale time to hours
ariv_transformed <- ariv_transformed %>%
  mutate(time = time / (60*24))

####################################################

# Prep for primary analysis (4 states, O2 / IMV / ICU DC / Death)
# Steps:
# 1) convert to wide
# 2) filter data consistent with absorbing states
# 3) add final transitions to censorship or death after IMV and ICU discharge
###### save this as raw data so we can look at how many obs per patient etc ####
# 4) add in relevant baseline covariates
# 5) impute missing data

# 1) convert to wide 
ariv_wide <- 
  ariv_transformed %>%
  group_by(stay_id) %>%
  pivot_wider(names_from = variable,
              values_from = value) %>%
  # forward fill o2 device and fio2 and pressor because they 
  # only change when the clinician changes them
  fill(o2_device, fio2, pressor) %>%
  mutate(pressor = ifelse(is.na(pressor), 0, pressor),
         icu_dc = ifelse(is.na(icu_dc), 0, icu_dc),
         death = ifelse(is.na(death), 0, death),
         state = case_when(
           death == 1 ~ 4,
           icu_dc == 1 ~ 3,
           o2_device == 6 ~ 2,
           o2_device < 6 ~ 1)) %>%
  filter(time >= 0) %>%
  arrange(stay_id, time)

# 2) Filter data to make consistent with absorbing states

# death before IMV

death_times <-
  ariv_wide %>%
  filter(death == 1) %>%
  summarise(deathtime = first(time))


ariv_wide <-
  ariv_wide %>%
  left_join(death_times, 
            by = "stay_id") %>%
  # then filter to the times <= deathtime
  # (and don't drop all the patients that never died)
  filter(time <= deathtime | is.na(deathtime)) %>%
  select(-deathtime)

# ICU discharge before death 

dc_times <-
  ariv_wide %>%
  filter(icu_dc == 1) %>%
  summarise(dctime = first(time))

ariv_wide <-
  ariv_wide %>%
  left_join(dc_times, 
            by = "stay_id") %>%
  # then filter to the times <= deathtime
  # (and don't drop all the patients that never died)
  filter(time <= dctime | is.na(dctime)) %>%
  select(-dctime)

# No time-varying covariates after IMV 

imv_times <-
  ariv_wide %>%
  filter(state == 2) %>%
  summarise(imvtime = first(time))

ariv_wide <-
  ariv_wide %>%
  left_join(imv_times, 
            by = "stay_id") %>%
  # then filter to the times <= deathtime
  # (and don't drop all the patients that never died)
  filter(time <= imvtime | is.na(imvtime)) %>%
  select(-imvtime)

# 3) add in final transitions after IMV or ICU discharge or censorship 96h

temp <- ariv_wide
ariv_wide <- temp

end <- 28

last_state <- ariv_wide %>%
  group_by(stay_id) %>%
  summarise(state = last(state),
            lasttime = last(time))

# have a look - but remember this is last state as of 
# the earliest of 96h, IMV, ICU discharge, and death
table(last_state$state)
round(table(last_state$state)/nrow(last_state),2)

# add in time of icu discharge and death
# icu_discharge_time and time_to_death are relative 
# to admission so subtract time to eligibility 
last_state <- 
  last_state %>%
  left_join(select(ariv_b, stay_id, 
                   icu_discharge_time, time_to_death, 
                   eligibility_time), 
            by = "stay_id") %>%
  mutate(time_to_death = case_when(
    is.na(time_to_death) ~ time_to_death,
    !is.na(time_to_death) ~ time_to_death - eligibility_time),
    dc_time = case_when(
      is.na(icu_discharge_time) ~ icu_discharge_time,
      !is.na(icu_discharge_time) ~ icu_discharge_time - eligibility_time)
  ) %>%
  mutate(dc_time = ifelse(dc_time > end*60*24, NA, dc_time),
         time_to_death = ifelse(time_to_death > end*60*24, NA, time_to_death)) %>%
  select(-eligibility_time, -icu_discharge_time) %>%
  # patients still in ICU at 28 days
  mutate(imv_till_28 = ifelse(is.na(time_to_death) &
                                is.na(dc_time) &
                                state == 2,
                              end*60*24,
                              NA),
         o2_till_28 = ifelse(is.na(time_to_death) &
                               is.na(dc_time) &
                               state == 1,
                             end*60*24,
                             NA)) %>%
  #filter(is.na(time_to_death), is.na(dc_time)) %>%
  # put events in order
  pivot_longer(cols = time_to_death:o2_till_28) %>%
  mutate(nexttime = value/(24*60)) %>%
  select(-value) %>%
  rename(nextstate = name) %>%
  mutate(nextstate = case_when(
    grepl("death", nextstate) ~ 4,
    grepl("dc", nextstate)    ~ 3,
    grepl("imv", nextstate)   ~ 2,
    grepl("o2",  nextstate)   ~ 1)) %>%
  arrange(stay_id, nexttime, desc(state)) %>%
  filter(!is.na(nexttime)) %>%
  group_by(stay_id) %>%
  mutate(priorstate = lag(nextstate)) %>%
  filter((priorstate <= nextstate) | is.na(priorstate)) %>%
  select(-priorstate)

# if no death occurs then patients remain in the same state
# probably will have to increase window beyond 96h
icudc28day <-
  last_state %>%
  arrange(stay_id, nexttime) %>%
  group_by(stay_id) %>%
  summarise(nextstate = last(nextstate),
            nexttime = end) %>%
  filter(nextstate != 4)

last_state <- 
  last_state %>% 
  bind_rows(icudc28day) %>%
  filter(!is.na(nexttime)) %>%
  arrange(stay_id, nexttime) %>%
  group_by(stay_id) %>%
  mutate(priorstate = lag(nextstate)) %>%
  filter(nextstate >= priorstate)

# remove all states after death
deathtimes2 <- 
  last_state %>%
  filter(nextstate == 4) %>%
  mutate(deathtime = nexttime) %>%
  select(-state, -nexttime, -nextstate)

last_state <- 
  last_state %>%
  left_join(deathtimes2, by = "stay_id") %>%
  filter(is.na(deathtime) | nexttime <= deathtime) %>%
  select(-deathtime, -state) %>%
  rename(time = nexttime,
         state = nextstate) 

ariv_wide <-
  bind_rows(ariv_wide, last_state) %>%
  arrange(stay_id, time) %>%
  # add covariates for transitions from discharge to death
  # and for imv to death
  group_by(stay_id) %>%
  mutate(imv_time = ifelse(state == 2,
                           time, NA),
         dc_time = ifelse(state == 3,
                          time, NA)
  ) %>%
  mutate(imv_time = ifelse(time < imv_time, NA, imv_time),
         dc_time = ifelse(time < dc_time, NA, dc_time)) %>%
  select(stay_id, time, state, heart_rate:pressor,imv_time, dc_time)

# remove impossible transitions that snuck in somehow
ariv_wide <-
  ariv_wide %>%
  group_by(stay_id) %>%
  arrange(stay_id, time) %>%
  mutate(lagstate = lag(state)) %>%
  mutate(lagstate = ifelse(is.na(lagstate), state, lagstate)) %>%
  filter(lagstate <= state & 
           lagstate != 4) %>%
  select(-lagstate)

# remove patients with initial absorbing state
initialabsorbing <- 
  ariv_wide %>%
  arrange(stay_id, time) %>%
  group_by(stay_id) %>%
  summarise(drop = first(state) != 1) # 6 patients total

ariv_wide <- 
  ariv_wide %>%
  left_join(initialabsorbing, by = "stay_id") %>%
  filter(drop == FALSE) %>%
  select(-drop)

# how many subjects with only one observation after processing?
#onlyone <- 
ariv_wide %>%
  group_by(stay_id) %>%
  summarise(count = n()) %>%
  arrange(count) %>%
  ungroup() %>%
  summarise(sum(count == 1)) # none

ariv_wide <- distinct(ariv_wide) # no duplicated rows

# 4) add in relevant baseline covariates

ariv_bl <- select(ariv_b,
                  stay_id, careunit, 
                  anchor_age,
                  anchor_year_group,
                  gender, insurance,
                  ethnicity,
                  chf, copd, cancer, dementia)

# missingness in baseline covariates? 

apply(is.na(ariv_bl), 2, mean)
# only in chf / copd / cancer / dementia
# and very minimal 
# we will assume those are 0

ariv_bl[is.na(ariv_bl)] <- 0

# combine with ariv_wide 

ariv_wide <- left_join(ariv_wide, ariv_bl, by = "stay_id")

# proportion of missing observations by variable
apply(is.na(ariv_wide), 2, mean)

# proportion missing after forward fill 
ariv_wide_ff <- ariv_wide %>%
  fill(heart_rate:pressor,
       .direction = "down")

apply(is.na(ariv_wide_ff), 2, mean)

# and fill in the mean for the tiny number  left unfilled after that
ariv_wide_ff <- ariv_wide_ff %>%  
  mutate(heart_rate = ifelse(is.na(heart_rate), 0, heart_rate),
         resp_rate = ifelse(is.na(resp_rate), 0, resp_rate),
         spo2 = ifelse(is.na(spo2), 0, spo2),
         GCS = ifelse(is.na(GCS), 0, GCS),
         wob = ifelse(is.na(wob), 0, wob),
         imv_time = ifelse(is.na(imv_time), -1, imv_time),
         dc_time = ifelse(is.na(dc_time), 0, dc_time))

apply(is.na(ariv_wide_ff), 2, mean)

# prep for Stan

# categorical covariates need to be broken down to binary

ariv_wide_stan <- 
  ariv_wide_ff %>%
  mutate(ra = ifelse(o2_device == 0,1,0),
         np = ifelse(o2_device == 1,1,0),
         fm = ifelse(o2_device == 2,1,0),
         nrb = ifelse(o2_device ==3,1,0),
         #hfnc=ifelse(o2_device == 4,1,0), #hfnc reference category
         niv =ifelse(o2_device == 5,1,0),
         G1 = ifelse(anchor_year_group == "2008 - 2010",1,0),
         G2 = ifelse(anchor_year_group == "2011 - 2013",1,0),
         G3 = ifelse(anchor_year_group == "2014 - 2016",1,0),
         #G4 = ifelse(anchor_year_group == "2017-2019",1,0) # reference category
         black = ifelse(ethnicity == "BLACK/AFRICAN AMERICAN", 1, 0),
         asian = ifelse(ethnicity == "ASIAN", 1, 0),
         hispanic = ifelse(ethnicity == "HISPANIC/LATINO",1, 0),
         # white reference category
         cardiac_icu = ifelse(careunit %in% c(
           "Coronary Care Unit (CCU)",
           "Cardiac Vascular Intensive Care Unit (CVICU)"),
           1,0),
         neuro_trauma_icu = ifelse(careunit %in% c(
           "Trauma SICU (TSICU)",
           "Neuro Surgical Intensive Care Unit (Neuro SICU)",
           "Neuro Intermediate",
           "Neuro Stepdown"), 1,0)
         # med-surg reference
  ) %>%
  select(-o2_device, -anchor_year_group, -careunit) %>%
  select(stay_id:pressor, ra:niv, anchor_age:neuro_trauma_icu, imv_time, dc_time) %>%
  mutate(gender = ifelse(gender == "M", 1, 0), # female reference category
         medicaid = ifelse(insurance == "Medicaid",1,0),
         otherins = ifelse(insurance == "Other", 1, 0) # MEdicare reference category
  ) %>%
  select(-insurance, -ethnicity)

saveRDS(ariv_wide_stan, "ariv_mimic_stan.rds")

################################################################################
#
#
#
#
# eICU segment
#
#

ariv_b <- read_parquet(
  file = "CreateCohort/ARIVe_eICU_Baseline"
) %>%
  # filter out to known race/ethnicities
  filter(ethnicity %in% c("African American",
                          "Asian",
                          "Caucasian",
                          "Hispanic")) 



ariv <- 
  read_parquet(
    file = "CreateCohort/ARIVe_eICU_Timevarying"
  ) %>% arrange(patientunitstayid, charttime, variable) %>% 
  filter(patientunitstayid %in% ariv_b$patientunitstayid) %>%
  # not using sbp
  filter(variable != "sbp") %>%
  filter(
    (variable == "fio2_recorded" & (
      (value >= 0.21 & value <= 1) |
        (value >= 21 & value <= 100)
    )
    ) | 
      (variable == "fio2_calculated" & (
        (value >= 0.21 & value <= 1) |
          (value >= 21 & value <= 100)
      )
      ) | 
      (variable == "heart_rate" & (
        value > 0 & value < 300
      )
      ) |
      (variable == "resp_rate" & (
        value > 0 & value < 100
      )
      ) |
      (variable == "gcs" & (
        value >= 3 & value <= 15
      )) |
      (variable == "spo2" & (value > 0 & value <= 100))
    | (
      variable %in% c("pressor", "o2_device")
    )) %>%
  rename(time = charttime) %>%
  # remove a patient who died at moment of eligibility (according to the data...)
  filter(patientunitstayid != 3066737)

# combine recorded and calculated fio2

fio2 <- 
  ariv %>%
  filter(grepl("fio2", variable)) %>%
  # take worst value where there are simultaneous duplicate entries in each category
  # ie if there are two simultaneous "fio2_recorded" then take the worst
  # there are 4572 (out of 491930) such rows ie < 1%
  group_by(patientunitstayid, time, variable) %>%
  summarise(value = max(value)) %>%
  ungroup() %>%
  group_by(patientunitstayid) %>%
  # now combine by using calculated only where no recorded value exists
  pivot_wider(values_from = value, 
              names_from = variable) %>%
  mutate(fio2 = case_when(
    !is.na(fio2_recorded) ~ fio2_recorded,
    is.na(fio2_recorded)  ~ fio2_calculated
  )) %>%
  select(-fio2_calculated, -fio2_recorded) %>%
  rename(value = fio2) %>%
  mutate(variable = "fio2") %>%
  select(patientunitstayid, time, variable, value)

ariv <- ariv %>%
  filter(!(grepl("fio2", variable))) %>%
  bind_rows(fio2) %>%
  arrange(patientunitstayid, time, variable)


log_group <- c("heart_rate", "resp_rate")

ariv_transformed <- ariv
ariv_transformed$value <- apply(ariv, 1, transform, log_group)

# how many NaN introduced by log?
sum(is.na(ariv_transformed$value)) # 33
length(ariv_transformed$value) # > 4.8 million, so filter them out
ariv_transformed <- filter(ariv_transformed, 
                           !is.na(value))

# rescale time to days
ariv_transformed <- ariv_transformed %>%
  mutate(time = time / (60*24))


# make into wide-r format 
# with a row for every patient and every observation time
# a column for every variable


####################################################

# Prep for primary analysis (4 states, O2 / IMV / ICU DC / Death)
# Steps:
# 1) convert to wide
# 2) filter data consistent with absorbing states
# 3) add final transitions to censorship or death after IMV and ICU discharge
###### save this as raw data so we can look at how many obs per patient etc ####
# 4) add in relevant baseline covariates
# 5) impute missing data

doubles <- ariv_transformed %>%
  group_by(patientunitstayid, time) %>%
  pivot_wider(names_from = variable,
              values_from = value,
              values_fn = length) 

apply(doubles[,-c(1:2)] > 1, 2, sum, na.rm = T)

# oxygen device simultaneous measurements easily resolved by taking the max
# (it's ordinal)
# take the max value for the others (because so few it doesn't matter)

# 1) convert to wide 
ariv_wide <- 
  ariv_transformed %>%
  group_by(patientunitstayid) %>%
  pivot_wider(names_from = variable,
              values_from = value,
              values_fn = max) %>% 
  # forward fill o2 device and fio2 and pressor because they 
  # only change when the clinician changes them
  fill(o2_device, fio2, pressor) %>%
  mutate(pressor = ifelse(is.na(pressor), 0, pressor),
         state = case_when(
           o2_device == 6 ~ 2,
           o2_device < 6 ~ 1)) %>%
  arrange(patientunitstayid, time)

# 2) Filter data to make consistent with absorbing states

# death and discharge

death_dc <-
  ariv_b %>%
  select(patientunitstayid, 
         eligibletime,
         unitdischargeoffset,
         unitdischargelocation,
         hospitaldischargeoffset,
         hospitaldischargestatus) %>%
  mutate(icu_dc_alive = case_when(
    unitdischargelocation == "Death" ~ 0,
    unitdischargelocation != "Death" ~ 1),
  ) %>%
  mutate(icu_dc = ifelse(icu_dc_alive == 1,
                         unitdischargeoffset - eligibletime,
                         NA),
         death = ifelse(hospitaldischargestatus == "Expired",
                        hospitaldischargeoffset - eligibletime,
                        NA)) %>%
  select(patientunitstayid,
         icu_dc,
         death) %>%
  pivot_longer(icu_dc:death,
               names_to = "state",
               values_to = "time") %>%
  filter(!is.na(time)) %>%
  mutate(state = case_when(
    state == "death" ~ 4,
    state == "icu_dc"~ 3
  ),
  time = time/(60*24))

# 3) add add dc and death transitions
# drop timevarying observations after discharge
# drop all observations after death

ariv_wide <-  ariv_wide %>%
  bind_rows(death_dc) %>%
  arrange(patientunitstayid, time)

death_times <- 
  death_dc %>%
  filter(state == 4) %>%
  select(patientunitstayid, time) %>%
  rename(deathtime = time)

dc_times <- 
  death_dc %>%
  filter(state == 3) %>%
  select(patientunitstayid, time) %>%
  rename(dctime = time)

ariv_wide <- 
  ariv_wide %>%
  left_join(dc_times,
            by = "patientunitstayid") %>%
  filter(is.na(dctime) | state >= 3 | time <= dctime) %>%
  select(-dctime) %>%
  left_join(death_times,
            by = "patientunitstayid") %>%
  filter(is.na(deathtime) | 
           time <= deathtime)

# add in censored observations at day 28

day28state = ariv_wide %>% 
  filter(time <= 28) %>%
  group_by(patientunitstayid) %>%
  summarise(laststate = last(state),
            lasttime = last(time)) %>%
  mutate(state = ifelse(laststate == 4, NA, laststate),
         time = 28) %>%
  select(-laststate, -lasttime) %>%
  filter(!is.na(state))

# drop all the data beyond day 28

ariv_wide <- 
  ariv_wide %>%
  bind_rows(day28state) %>%
  filter(time <= 28)

# 4) add in relevant baseline covariates

ariv_bl <- select(ariv_b,
                  patientunitstayid, unittype, 
                  age,
                  region, teachingstatus,
                  gender, 
                  ethnicity,
                  chf, copd, cancer, dementia) %>%
  mutate(age = ifelse(age == "> 89", 90,
                      as.numeric(age))) %>%
  filter(age >= 18) # 14 patients without age or with age = 0

apply(is.na(ariv_bl), 2, sum) # only region is missing so just add missing category

ariv_bl <- mutate(ariv_bl,
                  region = ifelse(is.na(region), "Unknown", region))

# combine with ariv_wide 
# and filter out the patients who were dropped for no age

ariv_wide <- filter(ariv_wide, 
                    patientunitstayid %in% ariv_bl$patientunitstayid) 
ariv_wide <- left_join(ariv_wide, ariv_bl, by = "patientunitstayid")

# proportion of missing observations by variable 
apply(is.na(filter(ariv_wide, time >= 0)), 2, mean)

# proportion missing after forward fill 
ariv_wide_ff <- ariv_wide %>%
  fill(heart_rate:pressor,
       .direction = "down") %>%
  filter(time >= 0)

apply(is.na(ariv_wide_ff), 2, mean)

# and fill in the mean for the tiny number  left unfilled after that
# GCS is 0.21 missing after that but everything else < 5%
ariv_wide_ff <- ariv_wide_ff %>%  
  mutate(heart_rate = ifelse(is.na(heart_rate), 0, heart_rate),
         resp_rate = ifelse(is.na(resp_rate), 0, resp_rate),
         spo2 = ifelse(is.na(spo2), 0, spo2),
         gcs = ifelse(is.na(gcs), 0, gcs)) %>%
  filter(time >= 0)

apply(is.na(ariv_wide_ff), 2, mean)

# drop tiny number of rows with no state, fio2 or o2 device (0.6% and 0.1%)

ariv_wide_ff <- filter(ariv_wide_ff,
                       !is.na(state), !is.na(o2_device),
                       !is.na(fio2))

# prep for Stan

# categorical covariates need to be broken down to binary

ariv_wide_stan <- 
  ariv_wide_ff %>%
  mutate(ra = ifelse(o2_device == 0,1,0),
         np = ifelse(o2_device == 1,1,0),
         fm = ifelse(o2_device == 2,1,0),
         nrb = ifelse(o2_device ==3,1,0),
         #hfnc=ifelse(o2_device == 4,1,0), #hfnc reference category
         niv =ifelse(o2_device == 5,1,0),
         R1 = ifelse(region == "Midwest",1,0),
         R2 = ifelse(region == "West",1,0),
         R3 = ifelse(region == "South",1,0),
         R4 = ifelse(region == "Unknown",1,0),
         #G4 = ifelse(anchor_year_group == "2017-2019",1,0) # reference category
         black = ifelse(ethnicity == "African American", 1, 0),
         asian = ifelse(ethnicity == "Asian", 1, 0),
         hispanic = ifelse(ethnicity == "Hispanic",1, 0),
         # white reference category
        cardiac_icu = ifelse(unittype %in% c(
          "CTICU","Cardiac ICU", "CCU-CTICU","CSICU"),
          1,0),
        neuro_trauma_icu = ifelse(unittype== "Neuro ICU", 1, 0)
  ) %>%
  select(-o2_device, -region, -unittype) %>%
  mutate(imv_time = ifelse(state == 2, time, NA),
         dc_time = ifelse(state == 3, time, NA)) %>%
  select(patientunitstayid, state, time, heart_rate:pressor, 
         ra:niv, black:neuro_trauma_icu, age:dementia,
         R1:R4, imv_time, dc_time) %>%
  mutate(gender = ifelse((gender == "Male" | gender == "Unknown"), 1, 0), # female reference category
         teachingstatus = ifelse(teachingstatus == TRUE,1,0))

saveRDS(ariv_wide_stan, "ariv_eicu_stan.rds")

################################################################################

# put the two databases together

ariv_eicu  <- readRDS("ariv_eicu_stan.rds")
ariv_mimic <- readRDS("ariv_mimic_stan.rds")

# fix names

names(ariv_mimic)
names(ariv_eicu)

ariv_mimic <- 
  ariv_mimic %>%
  rename(id = stay_id,
         gcs = GCS,
         age = anchor_age) %>%
  mutate(teachingstatus = 1,
         R1 = 0,
         R2 = 0,
         R3 = 0,
         R4 = 0,
         eicu = 0) %>%
  relocate(id, time, state, black:hispanic,
           age:G3,
           cardiac_icu, neuro_trauma_icu,
           teachingstatus:eicu,
           heart_rate:gcs,
           pressor:niv,
           imv_time, dc_time,
           wob, medicaid, otherins) %>%
  select(-wob, -medicaid, -otherins)

ariv_eicu <- 
  ariv_eicu %>%
  rename(id = patientunitstayid) %>%
  mutate(G1 = 0,
         G2 = 0,
         G3 = 1,
         eicu = 1) %>%
  relocate(id, time, state, 
           black:hispanic,
           age, gender, 
           chf:dementia,
           G1:G3,
           cardiac_icu, neuro_trauma_icu,
           teachingstatus,
           R1:R4, eicu,
           heart_rate, resp_rate, spo2, fio2, gcs, pressor, ra:niv, 
           imv_time, dc_time) %>%
    select(-ethnicity)

# combine

ariv_joint <- bind_rows(ariv_mimic, ariv_eicu)

# scale age by mean and sd
meanage <- mean(ariv_joint$age) # 66.6
sdage <- sd(ariv_joint$age) # 15.46

ariv_joint$age <- (ariv_joint$age - 66.6)/15.46

# center and scale continuous predictors

cont_vars <- c("spo2","heart_rate","resp_rate","fio2", "gcs")
cont_var_num <- which(names(ariv_joint) %in% cont_vars)
mean_sd <- data.frame(variable = names(ariv_joint)[cont_var_num],
                      mean = apply(ariv_joint[,cont_var_num], 2, mean),
                      sd = apply(ariv_joint[,cont_var_num], 2, sd))

saveRDS(mean_sd, file = "CreateCohort/joint_mean_sd.rds")

for (i in 1:length(cont_var_num)){
  ariv_joint[,cont_var_num[i]] = 
    (ariv_joint[,cont_var_num[i]] - mean_sd$mean[i])/mean_sd$sd[i]
}

saveRDS(ariv_joint, "ariv_joint.rds")
ariv_joint <- readRDS("ariv_joint.rds")
# need to change to cause-specific hazards format

state_to_state = function(df, state1, state2){
  df %>%
    group_by(id) %>%
    arrange(id, time) %>%
    rename(from = state,
           start = time) %>%
    mutate(to = lead(from),
           stop = lead(start)) %>%
    mutate(status = ifelse(from == state1 & to == state2, 1, 0)) %>% 
    filter(from == state1) %>%
    mutate(to = state2) %>% 
    filter(!is.na(stop))
}

remove_timevarying = function(df, start){
  if(start == 1) {
    df %>%
    group_by(id) %>%
    summarise(start = first(start),
              stop = last(stop),
              status = max(status),
              black = first(black),
              asian = first(asian),
              hispanic = first(hispanic),
              age = first(age),
              gender = first(gender),
              chf = first(chf),
              copd = first(copd),
              cancer = first(cancer),
              dementia = first(dementia),
              G1 = first(G1),
              G2 = first(G2),
              G3 = first(G3),
              cardiac_icu = first(cardiac_icu),
              neuro_trauma_icu = first(neuro_trauma_icu),
              teachingstatus = first(teachingstatus),
              R1 = first(R1),
              R2 = first(R2),
              R3 = first(R3),
              R4 = first(R4),
              eicu = first(eicu)
    )
}
  else if(start == 2){
    df %>%
      group_by(id) %>%
      summarise(start = first(start),
                stop = last(stop),
                status = max(status),
                black = first(black),
                asian = first(asian),
                hispanic = first(hispanic),
                age = first(age),
                gender = first(gender),
                chf = first(chf),
                copd = first(copd),
                cancer = first(cancer),
                dementia = first(dementia),
                G1 = first(G1),
                G2 = first(G2),
                G3 = first(G3),
                cardiac_icu = first(cardiac_icu),
                neuro_trauma_icu = first(neuro_trauma_icu),
                teachingstatus = first(teachingstatus),
                R1 = first(R1),
                R2 = first(R2),
                R3 = first(R3),
                R4 = first(R4),
                eicu = first(eicu),
                imv_time = max(imv_time, na.rm = T))
  } else {
    df %>%
      group_by(id) %>%
      summarise(start = first(start),
                stop = last(stop),
                status = max(status),
                black = first(black),
                asian = first(asian),
                hispanic = first(hispanic),
                age = first(age),
                gender = first(gender),
                chf = first(chf),
                copd = first(copd),
                cancer = first(cancer),
                dementia = first(dementia),
                G1 = first(G1),
                G2 = first(G2),
                G3 = first(G3),
                cardiac_icu = first(cardiac_icu),
                neuro_trauma_icu = first(neuro_trauma_icu),
                teachingstatus = first(teachingstatus),
                R1 = first(R1),
                R2 = first(R2),
                R3 = first(R3),
                R4 = first(R4),
                eicu = first(eicu),
                dc_time = max(dc_time, na.rm = T))
  }
}

# arive_xy.rds is the data for state x to state y

ariv_joint %>%
  state_to_state(1,2) %>% #3809 events so keep all 11 timevarying and 22 baseline
  select(-from, -to, -imv_time, -dc_time) %>%
  saveRDS("arive_joint_12.rds")

# non timevarying version
ariv_joint %>%
  state_to_state(1,2) %>% #3809 events so keep all 11 timevarying and 22 baseline
  remove_timevarying(1) %>%
  saveRDS("arive_joint_12_notv.rds")


ariv_joint %>%
  state_to_state(1,3) %>% 
  #remove_timevarying(1) %>% # 34615 events so let's do timevarying
  select(-from, -to, -imv_time, -dc_time) %>%
  # remove rows with no duration ie start == stop
  filter(stop != start) %>% 
  saveRDS("arive_joint_13.rds")

ariv_joint %>%
  state_to_state(1,4) %>%  # 1416 so let's do timevarying 
  select(-from, -to, -imv_time, -dc_time) %>%
  # remove rows with no duration ie start == stop
  filter(stop != start) %>% 
  saveRDS("arive_joint_14.rds")

ariv_joint %>%
  state_to_state(2,3) %>%
  remove_timevarying(2) %>%  # 1984, keep all 22 baseline
  mutate(imv_time = log(imv_time)) %>%
  saveRDS("arive_joint_23.rds")

ariv_joint %>%
  state_to_state(2,4) %>%
  remove_timevarying(2) %>%  #260 so take 9, 260/9 = 29 events per variable
#  select(id:gender, eicu, cardiac_icu, neuro_trauma_icu, imv_time) %>%
  mutate(imv_time = log(imv_time)) %>%  
  saveRDS("arive_joint_24.rds")

ariv_joint %>%
    state_to_state(3,4) %>%
    remove_timevarying(3) %>% # 1393 so keep all 22 baseline
  mutate(dc_time = log(dc_time)) %>%  
  saveRDS("arive_joint_34.rds")
  

###################################################

# for table 1:
# oxygen devices used ever

tbl1_means <- readRDS("data/arive_joint_12.rds") %>%
  group_by(id) %>%
  mutate(hfnc = ifelse(niv == 0 &
                         nrb == 0 &
                         fm == 0 &
                         np == 0 &
                         ra == 0,1,0),
         white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0),
         R0 = ifelse(R1 == 0 & R2 == 0 & R3 == 0 & R4 == 0, 1, 0),
         G0 = ifelse(G1 == 0 & G2 == 0 & G3 == 0, 1, 0)) %>%
  pivot_longer(cols = c(black, asian, hispanic,white)) %>%
  filter(value == 1) %>%
  summarise(
    age = first(age)*15.46 + 66.6,
    sex = 1-max(gender),
    chf = max(chf),
    copd = max(copd),
    cancer = max(cancer),
    dementia = max(dementia),
    G0 = max(G0),
    G1 = max(G1),
    G2 = max(G2),
    G3 = max(G3),
    R0 = max(R0),
    R1 = max(R1),
    R2 = max(R2),
    R3 = max(R3),
    R4 = max(R4),
    teachingstatus = max(teachingstatus),
    eicu = max(eicu),
    cardiac = max(cardiac_icu),
    neurotrauma = max(neuro_trauma_icu),
    niv = max(niv),
    hfnc = max(hfnc),
    nrb = max(nrb),
    fm = max(fm),
    np = max(np),
    ra = max(ra),
    race = first(name)) %>%
  ungroup() %>%
  group_by(race) %>%
  summarise(across(age:ra, mean))  #%>%
  #write.csv("table1_means.csv")

#tbl1_mean_wc <- 
readRDS("data/arive_joint_12.rds") %>%
  group_by(id) %>%
  mutate(hfnc = ifelse(niv == 0 &
                         nrb == 0 &
                         fm == 0 &
                         np == 0 &
                         ra == 0,1,0),
         white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0),
         R0 = ifelse(R1 == 0 & R2 == 0 & R3 == 0 & R4 == 0, 1, 0),
         G0 = ifelse(G1 == 0 & G2 == 0 & G3 == 0, 1, 0)) %>%
#  pivot_longer(cols = c(black, asian, hispanic,white)) %>%
#  filter(value == 1) %>%
  summarise(
    age = first(age)*15.46 + 66.6,
    sex = 1-max(gender),
    chf = max(chf),
    copd = max(copd),
    cancer = max(cancer),
    dementia = max(dementia),
    G0 = max(G0),
    G1 = max(G1),
    G2 = max(G2),
    G3 = max(G3),
    R0 = max(R0),
    R1 = max(R1),
    R2 = max(R2),
    R3 = max(R3),
    R4 = max(R4),
    teachingstatus = max(teachingstatus),
    eicu = max(eicu),
    cardiac = max(cardiac_icu),
    neurotrauma = max(neuro_trauma_icu),
    niv = max(niv),
    hfnc = max(hfnc),
    nrb = max(nrb),
    fm = max(fm),
    np = max(np),
    ra = max(ra)) %>%
  ungroup() %>%
#  group_by(race) %>%
  summarise(across(age:ra, mean))  %>%
write.csv("tbl1_means_wc.csv")


tbl1_sums <- readRDS("data/arive_joint_12.rds") %>%
  group_by(id) %>%
  mutate(hfnc = ifelse(niv == 0 &
                         nrb == 0 &
                         fm == 0 &
                         np == 0 &
                         ra == 0,1,0),
         white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0),
         R0 = ifelse(R1 == 0 & R2 == 0 & R3 == 0 & R4 == 0, 1, 0),
         G0 = ifelse(G1 == 0 & G2 == 0 & G3 == 0, 1, 0)) %>%
  pivot_longer(cols = c(black, asian, hispanic,white)) %>%
  filter(value == 1) %>%
  summarise(
    age = first(age)*15.46 + 66.6,
    sex = 1-max(gender),
    chf = max(chf),
    copd = max(copd),
    cancer = max(cancer),
    dementia = max(dementia),
    G0 = max(G0),
    G1 = max(G1),
    G2 = max(G2),
    G3 = max(G3),
    R0 = max(R0),
    R1 = max(R1),
    R2 = max(R2),
    R3 = max(R3),
    R4 = max(R4),
    teachingstatus = max(teachingstatus),
    eicu = max(eicu),
    cardiac = max(cardiac_icu),
    neurotrauma = max(neuro_trauma_icu),
    niv = max(niv),
    hfnc = max(hfnc),
    nrb = max(nrb),
    fm = max(fm),
    np = max(np),
    ra = max(ra),
    race = first(name)) %>%
  ungroup() %>%
  group_by(race) %>%
  summarise(across(age:ra, sum)) 

readRDS("data/arive_joint_12.rds")%>%
  select(id, age) %>%
  group_by(id) %>%
  summarise(age = first(age)) %>%
  ungroup() %>%
  mutate(age = (age*15.4 + 66.6)) %>%
  summarise(median(age),
            quantile(age, 0.25),
            quantile(age, 0.75))

#tbl1_sums_wc <- 
  readRDS("data/arive_joint_12.rds") %>%
  group_by(id) %>%
  mutate(hfnc = ifelse(niv == 0 &
                         nrb == 0 &
                         fm == 0 &
                         np == 0 &
                         ra == 0,1,0),
         white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0),
         R0 = ifelse(R1 == 0 & R2 == 0 & R3 == 0 & R4 == 0, 1, 0),
         G0 = ifelse(G1 == 0 & G2 == 0 & G3 == 0, 1, 0)) %>%
  summarise(
    age = first(age)*15.46 + 66.6,
    sex = 1-max(gender),
    chf = max(chf),
    copd = max(copd),
    cancer = max(cancer),
    dementia = max(dementia),
    G0 = max(G0),
    G1 = max(G1),
    G2 = max(G2),
    G3 = max(G3),
    R0 = max(R0),
    R1 = max(R1),
    R2 = max(R2),
    R3 = max(R3),
    R4 = max(R4),
    teachingstatus = max(teachingstatus),
    eicu = max(eicu),
    cardiac = max(cardiac_icu),
    neurotrauma = max(neuro_trauma_icu),
    niv = max(niv),
    hfnc = max(hfnc),
    nrb = max(nrb),
    fm = max(fm),
    np = max(np),
    ra = max(ra)) %>%
  ungroup() %>%
  summarise(across(age:ra, sum)) %>%
    write_csv("tbl1_sums_wc.csv")


tbl1 <- 
  bind_rows(tbl1_means,tbl1_sums, .id = "name") %>%
  mutate(name = factor(name, levels = 1:2, labels = c("mean","sum"))) %>%
  pivot_longer(cols = age:ra,
               names_to = "variable",
               values_to = "value") %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(label = paste0(sum, " (", 
                          round(100*mean),
                          ")")) %>%
    select(race, variable, label) %>%
    pivot_wider(names_from = race, values_from = label) 

tbl1 %>% write.csv("tbl1_proportions.csv")

# age
  readRDS("data/arive_joint_12.rds") %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  distinct(id, age, white, black, hispanic, asian) %>%
  pivot_longer(cols = c(white, black, hispanic, asian)) %>%
    filter(value == 1) %>%
  mutate(age = age*15.4+66.6) %>%
    group_by(name) %>%
    summarise(median = round(median(age)),
              iqr25 = round(quantile(age, 0.25)),
              iqr75 = round(quantile(age, 0.75)))
  
  # A tibble: 4 x 4
##  name     median iqr25 iqr75
##  <chr>     <dbl> <dbl> <dbl>
#  1 asian        67    55    79
#  2 black        62    51    74
#  3 hispanic     68    53    80
#  4 white        69    58    80
  
  
# observed IMV and death probabilities
  
death14 <- 
  readRDS("data/arive_joint_14.rds") %>%
    filter(status == 1) %>%
    mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
    select(white, black, hispanic, asian) %>%
    pivot_longer(cols = c(white, black, hispanic, asian))  %>%
    ungroup() %>%
    filter(value == 1) %>%
    group_by(name) %>%
    summarise(count = sum(value))

death24 <-   readRDS("data/arive_joint_24.rds") %>%
  filter(status == 1) %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  select(white, black, hispanic, asian) %>%
  pivot_longer(cols = c(white, black, hispanic, asian))  %>%
  ungroup() %>%
  filter(value == 1) %>%
  group_by(name) %>%
  summarise(count = sum(value))

death34 <-   readRDS("data/arive_joint_34.rds") %>%
  filter(status == 1) %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  select(white, black, hispanic, asian) %>%
  pivot_longer(cols = c(white, black, hispanic, asian))  %>%
  ungroup() %>%
  filter(value == 1) %>%
  group_by(name) %>%
  summarise(count = sum(value))


bind_rows(death14, death24, death34, .id = "state") %>%
  group_by(name) %>%
  summarise(sum = sum(count)) %>%
  mutate(totals = c(854, 4255, 2050, 31923)) %>%
  mutate(prop = (totals-sum)/totals)

# imv

readRDS("data/arive_joint_12.rds") %>%
  filter(status == 1) %>%
  mutate(white = ifelse(black == 0 & hispanic == 0 & asian == 0, 1, 0)) %>%
  select(white, black, hispanic, asian) %>%
  pivot_longer(cols = c(white, black, hispanic, asian))  %>%
  ungroup() %>%
  filter(value == 1) %>%
  group_by(name) %>%
  summarise(count = sum(value)) %>%
  mutate(totals = c(854, 4255, 2050, 31923)) %>%
  mutate(prop = count/totals)
