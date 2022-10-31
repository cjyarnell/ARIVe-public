# functions for piecewise 4state Stan msm fit



#emod <- cmdstan_model("ARIVe_exponential_msm.stan")

fit_exponential_transition <- function(states){
  state1 = states[1]
  state2 = states[2]
  filename = paste0("arive_joint_", state1,
                    state2, ".rds")
  temp <- readRDS(filename)
  
  temp <- temp %>%
    mutate(intercept = 1,
           duration = stop-start) %>%
    ungroup() %>%
    select(-id, -start, -stop) %>%
    relocate(duration, status, intercept)
  
  stanlist = list(temp)
  
  standat <- list(
    P = ncol(stanlist[[1]])-2,
    N_cen = sum(stanlist[[1]]$status == 0),
    N_obs = sum(stanlist[[1]]$status == 1),
    X_cen = select(filter(stanlist[[1]], status == 0), -status, -duration),
    X_obs = select(filter(stanlist[[1]], status == 1), -status, -duration),
    y_cen = stanlist[[1]]$duration[stanlist[[1]]$status == 0],
    y_obs = stanlist[[1]]$duration[stanlist[[1]]$status == 1]
  )
  
  post <- emod$sample(data = standat,
                      chains = 4,
                      iter_warmup = 750,
                      iter_sampling = 250,
                      adapt_delta = 0.9)
  
  post$save_object(paste0("fits/post_joint",
                          state1, 
                          state2,
                          ".rds"))
  
}

# function to generate the 6 covariate sets
generate_X <- function(X){
  covlist <- list()
  
  covlist[[1]] <- X %>% 
    select(intercept:niv)
  
  covlist[[2]] <- X %>% 
    select(intercept:niv)
  
  covlist[[3]] <- X %>% 
    select(intercept:niv)
  
  covlist[[4]] <- X %>%
    select(intercept:eicu)
  
  covlist[[5]] <- X %>%
    select(intercept:eicu)
  
  covlist[[6]] <- X %>%
    select(intercept:eicu)
    
  names(covlist) <- c("X_pred12",
                     "X_pred13",
                     "X_pred14",
                     "X_pred23",
                     "X_pred24",
                     "X_pred34")
    
  covlist$P12 = ncol(covlist[[1]])
  covlist$P13 = ncol(covlist[[2]])
  covlist$P14 = ncol(covlist[[3]])
  covlist$P23 = ncol(covlist[[4]])+1
  covlist$P24 = ncol(covlist[[5]])+1
  covlist$P34 = ncol(covlist[[6]])+1

  covlist
}

# read models in 

read_models <- function(){
  models <- list()
  
    models[[1]] = readRDS("fits/post_joint12.rds")$draws(variables = "beta",
                                       format = "matrix")
    
    models[[2]] = 
      readRDS("fits/post_joint13.rds")$draws(variables = "beta",
                                       format = "matrix")
    models[[3]] = 
      readRDS("fits/post_joint14.rds")$draws(variables = "beta",
                                       format = "matrix")
    
    models[[4]] = 
      readRDS("fits/post_joint23.rds")$draws(variables = "beta",
                                       format = "matrix")
    
    models[[5]] = 
      readRDS("fits/post_joint24.rds")$draws(variables = "beta",
                                       format = "matrix")
    
    models[[6]] = 
      readRDS("fits/post_joint34.rds")$draws(variables = "beta",
                                       format = "matrix")
    
  models
}

# pick out the ith iter of each model

one_iter <- function(betas, i){
  list(
    as.numeric(betas[[1]][i,]),
    as.numeric(betas[[2]][i,]),
    as.numeric(betas[[3]][i,]),
    as.numeric(betas[[4]][i,]),
    as.numeric(betas[[5]][i,]),
    as.numeric(betas[[6]][i,])  )
}

# make Xbeta vector

make_Xbeta <- function(X, betas){
  Xbeta = rep(0, 6)
  X[[4]]$imv_time = 0
  X[[5]]$imv_time = 0
  X[[6]]$dc_time  = 0
  for(i in 1:6){
    Xbeta[i] = betas[[i]] %*% t(X[[i]])
    if(i < 4){Xbeta[i] = exp(Xbeta[i])}
  }
  Xbeta
}

make_Xbeta_counterfactual <- function(X, betas){
  Xbeta = rep(0, 6)
  X[[4]]$imv_time = 0
  X[[5]]$imv_time = 0
  X[[6]]$dc_time  = 0
  
  X[[1]] <- X[[1]] %>%
    mutate(black = 0,
           asian = 0,
           hispanic = 0)
  
  for(i in 1:6){
    Xbeta[i] = betas[[i]] %*% t(X[[i]])
    if(i < 4){Xbeta[i] = exp(Xbeta[i])}
  }
  Xbeta
}

# function to simulate a patient's path through the model

simulate_path <- function(Xbeta, imv23beta, 
                          imv24beta,
                          dc34beta){
  
  from = list()
  from[[1]] = 1
  times = list()
  to = list()
  status = list()
  
  max_time = 28
  
  # betas is a list of the 6 models
  # X is a list of the 6 covariate sets
  
  # First transition
  
  t1s <- c(rexp(1, rate = Xbeta[1]),
           rexp(1, rate = Xbeta[2]),
           rexp(1, rate = Xbeta[3]))
  
  if(min(t1s) > max_time){
    to[[1]] = 1
    status[[1]] = 0 
    times[[1]] = max_time
    return(bind_cols(from = as.numeric(from), 
                     to = as.numeric(to), 
                     times = as.numeric(times), 
                     status = as.numeric(status)))
  } else {
    status[[1]]= 1
    to[[1]]    = 1 + which.min(t1s)
    times[[1]] = min(t1s)
  }
  
  # Second transition (if needed)
  
  if(to[[1]] == 2){
    
    Xbeta[4] = exp(Xbeta[4] + imv23beta*log(times[[1]]))
    Xbeta[5] = exp(Xbeta[5] + imv24beta*log(times[[1]]))
    
    from[[2]] = to[[1]]
    t2s <- c(rexp(1, rate = Xbeta[4]),
             rexp(1, rate = Xbeta[5]))
    if(min(t2s) + times[[1]] > max_time){
      to[[2]] = 2
      status[[2]] = 0
      times[[2]] = max_time - times[[1]]
      return(bind_cols(from = as.numeric(from), 
                       to = as.numeric(to), 
                       times = as.numeric(times), 
                       status = as.numeric(status)))
    } else {
      status[[2]]= 1
      to[[2]]    = 2 + which.min(t2s)
      times[[2]] = min(t2s)
    }
  } else if(to[[1]] == 3){
    
    Xbeta[6] = exp(Xbeta[6] + dc34beta*log(times[[1]]))
    
    from[[2]] = to[[1]]
    t3s <- rexp(1, rate = Xbeta[6])
    if(min(t3s) + times[[1]] > max_time){
      to[[2]] = 3
      status[[2]] = 0
      times[[2]] = max_time - times[[1]]
      return(bind_cols(from = as.numeric(from), 
                       to = as.numeric(to), 
                       times = as.numeric(times), 
                       status = as.numeric(status)))
    } else {
      status[[2]] = 1
      to[[2]]     = 4
      times[[2]]  = min(t3s)
      return(bind_cols(from = as.numeric(from), 
                       to = as.numeric(to), 
                       times = as.numeric(times), 
                       status = as.numeric(status)))
    }
  } else { #to[[1]] == 4 or death
    return(bind_cols(from = as.numeric(from), 
                     to = as.numeric(to), 
                     times = as.numeric(times), 
                     status = as.numeric(status)))
  }
  
  # third transition if needed
  if(to[[2]] == 3){
    
    Xbeta[6] = exp(Xbeta[6] + dc34beta*log(times[[2]]))
    
    from[[3]] = to[[2]]
    t3s <- rexp(1, rate = Xbeta[6])
    if(min(t3s) + times[[1]] + times[[2]] > max_time){
      to[[3]] = 3
      status[[3]] = 0
      times[[3]] = max_time - times[[1]] - times[[2]]
      return(bind_cols(from = as.numeric(from), 
                       to = as.numeric(to), 
                       times = as.numeric(times), 
                       status = as.numeric(status)))
    } else {
      status[[3]] = 1
      to[[3]]     = 4
      times[[3]]  = min(t3s)
      return(bind_cols(from = as.numeric(from), 
                       to = as.numeric(to), 
                       times = as.numeric(times), 
                       status = as.numeric(status)))
    }
  } else { #to[[2]] == 4 or death
    return(bind_cols(from = as.numeric(from), 
                     to = as.numeric(to), 
                     times = as.numeric(times), 
                     status = as.numeric(status)))
  }
}

# generate vector of reference covariates

reference_covariates <- function(){
    X <- data.frame(
      intercept = 1,
      black = 0,
      asian = 0,
      hispanic = 0,
      age = 0,
      gender = 0,
      chf = 0,
      copd = 0,
      cancer = 0,
      dementia = 0,
      G1 = 0,
      G2 = 0,
      G3 = 0,
      cardiac_icu = 0,
      neuro_trauma_icu = 0,
      teachingstatus = 1,
      R1 = 0,
      R2 = 0,
      R3 = 0,
      R4 = 0,
      eicu = 0,
      heart_rate = 0.79, #100 BPM
      resp_rate = 0.90, # 25 breaths per minute
      spo2 = -0.76, # 94% 
      fio2 = 1.66, # 80%
      GCS = 0,
      pressor = 0,
      ra = 0,
      np = 0,  # hfnc
      fm = 0,
      nrb = 0,
      niv = 0) 
  X
}

# P(alive at 28 days | A = a)

wholecohort_oneiter <- function(ptlist, models, i, traj){
  
    temp = list()
    for (n in 1:nrow(ptlist)){
      temp[[n]] <- generate_X(ptlist[n,-1]) %>%
        onept_oneiter(models, i, traj)
    }
    bind_rows(temp) %>%
      group_by(quantity, iter) %>%
      summarise(across(survival:imv, mean))
}

onept_oneiter <- function(X, models, i, traj){
    betas <- one_iter(models, i)
    Xbeta <- make_Xbeta(X, betas)
    Xbeta_counterfactual <- make_Xbeta_counterfactual(X, betas)
    temp1 <- list()
    temp2 <- list()
    for (j in 1:traj){
      temp1[[j]] <- bind_cols(data.frame(quantity = "total effect",
                                         iter = i, traj = j), 
                             simulate_path(Xbeta, 
                                           imv23beta = betas[[4]][22],
                                           imv24beta = betas[[5]][22],
                                           dc34beta  = betas[[6]][22]))
      temp2[[j]] <- bind_cols(data.frame(quantity = "counterfactual",
                                         iter = i, traj = j), 
                              simulate_path(Xbeta_counterfactual, 
                                            imv23beta = betas[[4]][22],
                                            imv24beta = betas[[5]][22],
                                            dc34beta  = betas[[6]][22]))
      
    }
    output <- 
      bind_rows(temp1,temp2) %>% 
      group_by(quantity, iter, traj) %>%
      summarise(death = last(to) == 4,
                state2 = max(to == 2)) %>%
      group_by(quantity, iter) %>%
      summarise(survival = 1-mean(death),
                imv = mean(state2))
    output
}


counterfactual_effect28days <- function(X, models, iter, traj){
  
  output <- list()
  for (i in 1:iter){
    betas <- one_iter(models, i)
    temp <- list()
    for (j in 1:traj){
      temp[[j]] <- bind_cols(data.frame(iter = i, traj = j), 
                             simulate_path_counterfactual(betas, X))
    }
    output[[i]] <- 
      bind_rows(temp) %>% 
      group_by(iter, traj) %>%
      summarise(death = last(to) == 4) %>%
      group_by(iter) %>%
      summarise(survival = 1-mean(death))
  }
  bind_rows(output)
}


timeweighted_mean <- function(df){
  
  ptlist <- readRDS(df) %>% 
    mutate(intercept = 1,
           duration = stop-start) %>%
    ungroup() %>%
    select(-start, -stop) %>%
    relocate(id, duration, status, intercept)
  
  ptlist[,25:ncol(ptlist)] = ptlist[,25:ncol(ptlist)] * ptlist$duration
  
  ptlist <- 
    group_by(ptlist, id) %>%
    summarise(duration = sum(duration),
              across(intercept:eicu, mean),
              across(heart_rate:niv, sum)) %>%
    mutate(heart_rate = heart_rate/duration,
           resp_rate = resp_rate/duration,
           spo2 = spo2/duration,
           fio2 = fio2/duration,
           gcs = gcs/duration,
           pressor = pressor/duration,
           ra = ra/duration,
           np = np/duration,
           fm = fm/duration,
           nrb = nrb/duration,
           niv = niv/duration
    )
}

timeweighted_mean_fromdf <- function(df){
  
  ptlist = df
  
  ptlist[,24:ncol(ptlist)] = ptlist[,24:ncol(ptlist)] * ptlist$duration
  
  ptlist <- 
    ptlist %>%
    summarise(duration = sum(duration),
              across(intercept:eicu, mean),
              across(heart_rate:niv, sum)) %>%
    mutate(heart_rate = heart_rate/duration,
           resp_rate = resp_rate/duration,
           spo2 = spo2/duration,
           fio2 = fio2/duration,
           gcs = gcs/duration,
           pressor = pressor/duration,
           ra = ra/duration,
           np = np/duration,
           fm = fm/duration,
           nrb = nrb/duration,
           niv = niv/duration
    )
}

ready_for_gq <- function(df){
  
    bind_rows(white = df,
                    black = df,
                    asian = df,
                    hispanic = df) %>%
    mutate(rn = row_number()) %>%
    mutate(black = ifelse(rn == 2, 1, 0),
           asian = ifelse(rn == 3, 1, 0),
           hispanic = ifelse(rn == 4, 1, 0))
  
}


state_by_time <- function(df, times){

cdf <- list()

for (i in 1:length(times)){
  cdf[[i]] <- 
  df %>% 
    filter(start < times[i]) %>%
    group_by(id) %>%
    summarise(state = last(state)) %>%
    group_by(state) %>%
    count() %>%
    mutate(time = times[i])
}
bind_rows(cdf)
}

invgammaprior <- function(n, shape, scale){
  ig <- rinvgamma(n, shape = shape, scale = scale)
  plot(density(ig))
  c(mean(ig < 0.5),
  mean(ig>2),
  mean(ig))
}

simweibull <- function(n, shape, scale){
  rw <- rweibull(n = n, 
                 shape = shape,
                 scale = scale)
  delta <- rep(1, n)
  plot(survfit(Surv(rw, delta) ~ 1 ), col=c('black','blue'),
       xlab='Time',ylab='Survival Probability', conf.int=T)
}

lognormalprior <- function(n, logmean, logsd){
  rn <- rnorm(n, mean = logmean, sd = logsd)
  plot(density(exp(rn)))
  c(mean(mean(exp(rn) < 0.5)), mean(exp(rn) > 2), mean(exp(rn)))
}

make_full_stanlist <- function(stanlist){

    full_stanlist <- list(
    P12 = ncol(stanlist[[1]])-2,
  N_cen12 = sum(stanlist[[1]]$status == 0),
  N_obs12 = sum(stanlist[[1]]$status == 1),
  X_cen12 = select(filter(stanlist[[1]], status == 0), -status, -duration),
  X_obs12 = select(filter(stanlist[[1]], status == 1), -status, -duration),
  y_cen12 = stanlist[[1]]$duration[stanlist[[1]]$status == 0],
  y_obs12 = stanlist[[1]]$duration[stanlist[[1]]$status == 1],
  
  P13 = ncol(stanlist[[2]])-2,
  N_cen13 = sum(stanlist[[2]]$status == 0),
  N_obs13 = sum(stanlist[[2]]$status == 1),
  X_cen13 = select(filter(stanlist[[2]], status == 0), -status, -duration),
  X_obs13 = select(filter(stanlist[[2]], status == 1), -status, -duration),
  y_cen13 = stanlist[[2]]$duration[stanlist[[2]]$status == 0],
  y_obs13 = stanlist[[2]]$duration[stanlist[[2]]$status == 1],
  
  P14 = ncol(stanlist[[3]])-2,
  N_cen14 = sum(stanlist[[3]]$status == 0),
  N_obs14 = sum(stanlist[[3]]$status == 1),
  X_cen14 = select(filter(stanlist[[3]], status == 0), -status, -duration),
  X_obs14 = select(filter(stanlist[[3]], status == 1), -status, -duration),
  y_cen14 = stanlist[[3]]$duration[stanlist[[3]]$status == 0],
  y_obs14 = stanlist[[3]]$duration[stanlist[[3]]$status == 1],
  
  P23 = ncol(stanlist[[4]])-2,
  N_cen23 = sum(stanlist[[4]]$status == 0),
  N_obs23 = sum(stanlist[[4]]$status == 1),
  X_cen23 = select(filter(stanlist[[4]], status == 0), -status, -duration),
  X_obs23 = select(filter(stanlist[[4]], status == 1), -status, -duration),
  y_cen23 = stanlist[[4]]$duration[stanlist[[4]]$status == 0],
  y_obs23 = stanlist[[4]]$duration[stanlist[[4]]$status == 1],
  
  P24 = ncol(stanlist[[5]])-2,
  N_cen24 = sum(stanlist[[5]]$status == 0),
  N_obs24 = sum(stanlist[[5]]$status == 1),
  X_cen24 = select(filter(stanlist[[5]], status == 0), -status, -duration),
  X_obs24 = select(filter(stanlist[[5]], status == 1), -status, -duration),
  y_cen24 = stanlist[[5]]$duration[stanlist[[5]]$status == 0],
  y_obs24 = stanlist[[5]]$duration[stanlist[[5]]$status == 1],
  
  P34 = ncol(stanlist[[6]])-2,
  N_cen34 = sum(stanlist[[6]]$status == 0),
  N_obs34 = sum(stanlist[[6]]$status == 1),
  X_cen34 = select(filter(stanlist[[6]], status == 0), -status, -duration),
  X_obs34 = select(filter(stanlist[[6]], status == 1), -status, -duration),
  y_cen34 = stanlist[[6]]$duration[stanlist[[6]]$status == 0],
  y_obs34 = stanlist[[6]]$duration[stanlist[[6]]$status == 1]
  )
    
    full_stanlist
  
}

read_data <- function(df){
  df %>%
    mutate(intercept = 1,
           duration = stop-start) %>%
    ungroup() %>%
    select(-id, -start, -stop) %>%
    relocate(duration, status, intercept)
}

output_hazard_ratio <- function(model, transition, names){
    draws <- model$draws(variables = transition)
    data.frame(variable = names,
               mean = apply(draws, 3, mean),
               lb95 = apply(draws, 3, quantile, 0.025),
               ub95 = apply(draws, 3, quantile, 0.975))
}

gqmod <- cmdstan_model("models/ARIVe_4state_msm_gq.stan")

single_patient_gq <- function(X, name){

    gqdat <- generate_X(X)
    
    gqdat$N_pred <- nrow(gqdat$X_pred12)
    gqdat$N_rep <- 10000

gq <- gqmod$generate_quantities(fitted_params = post,
                                         data = gqdat)

gq$save_object(file = paste0(name, ".rds"))

# imv 

imv <- gq$draws(variables = "P_imv28")
imvdf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(imv, 3, mean),
           lb95 = apply(imv, 3, quantile, 0.025),
           ub95 = apply(imv, 3, quantile, 0.975))

ttimv <- gq$draws(variables = "T_imv28")
ttimvdf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(ttimv, 3, mean),
           lb95 = apply(ttimv, 3, quantile, 0.025),
           ub95 = apply(ttimv, 3, quantile, 0.975))
# causal mediation

totaleffect <- gq$draws(variables = "P_surv28")
tedf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(totaleffect, 3, mean),
           lb95 = apply(totaleffect, 3, quantile, 0.025),
           ub95 = apply(totaleffect, 3, quantile, 0.975))

indirecteffect <- gq$draws(variables = "P_cntr28") - 
                  gq$draws(variables = c("P_surv28[2]", "P_surv28[3]","P_surv28[4]"))
iedf <- data.frame(subgroup = c("Black","Asian","Hispanic"),
           mean = apply(indirecteffect, 3, mean),
           lb95 = apply(indirecteffect, 3, quantile, 0.025),
           ub95 = apply(indirecteffect, 3, quantile, 0.975))
    
    list(imv = imvdf,
        ttimv = ttimvdf,
        te = tedf,
        ie = iedf)
}

package_gq <- function(gq){
    # imv 

imv <- gq$draws(variables = "P_imv28")
imvdf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(imv, 3, mean),
           lb95 = apply(imv, 3, quantile, 0.025),
           ub95 = apply(imv, 3, quantile, 0.975))

ttimv <- gq$draws(variables = "T_imv28")
ttimvdf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(ttimv, 3, mean),
           lb95 = apply(ttimv, 3, quantile, 0.025),
           ub95 = apply(ttimv, 3, quantile, 0.975))
# causal mediation

totaleffect <- gq$draws(variables = "P_surv28")
tedf <- data.frame(subgroup = c("White","Black","Asian","Hispanic"),
           mean = apply(totaleffect, 3, mean),
           lb95 = apply(totaleffect, 3, quantile, 0.025),
           ub95 = apply(totaleffect, 3, quantile, 0.975))

indirecteffect <- gq$draws(variables = "P_cntr28") - 
                  gq$draws(variables = c("P_surv28[2]", "P_surv28[3]","P_surv28[4]"))
iedf <- data.frame(subgroup = c("Black","Asian","Hispanic"),
           mean = apply(indirecteffect, 3, mean),
           lb95 = apply(indirecteffect, 3, quantile, 0.025),
           ub95 = apply(indirecteffect, 3, quantile, 0.975))

imvdiff <- bind_rows(data.frame(subgroup = "Black",
                      value = -as.numeric(gq$draws(variables = "P_imv28[2]") - 
                                gq$draws(variables = "P_imv28[1]"))),
                      data.frame(subgroup = "Asian",
                      value = -as.numeric(gq$draws(variables = "P_imv28[3]") - 
                                gq$draws(variables = "P_imv28[1]"))),
                     data.frame(subgroup = "Hispanic",
                      value = -as.numeric(gq$draws(variables = "P_imv28[4]") - 
                                gq$draws(variables = "P_imv28[1]"))))

imvdiffdf <- 
    imvdiff %>%
        group_by(subgroup) %>%
        summarise(mean = mean(value),
                  lb95 = quantile(value, 0.025),
                  ub95 = quantile(value, 0.975),
                  p_less_imv = mean(value < 0))

p_harm <- data.frame(
            subgroup = c("Black","Asian","Hispanic"),
            prob = apply(indirecteffect>0, 3, mean))
    
    list(imv = imvdf,
        ttimv = ttimvdf,
        te = tedf,
        ie = iedf,
        imvdiff = imvdiffdf,
        p_harm = p_harm)
}


outcomes_by_race <- function(subset){

    # black race/ethnicity
    sb <- subset %>%
    filter(black == 1)
    sb <- generate_X(sb)
    sb$N_pred <- nrow(sb$X_pred12)
    sb$N_rep <- 10000

    postb <- gqmod$generate_quantities(fitted_params = postcsv,
                                     data = sb,
                                      parallel_chains = 4)
    
    # asian race/ethnicity
    sa <- subset %>%
    filter(asian == 1)
    sa <- generate_X(sa)
    sa$N_pred <- nrow(sa$X_pred12)
    sa$N_rep <- 10000

    posta <- gqmod$generate_quantities(fitted_params = postcsv,
                                     data = sa,
                                      parallel_chains = 4)
    

    # hispanic race/ethnicity
    sh <- subset %>%
    filter(hispanic == 1)
    sh <- generate_X(sh)
    sh$N_pred <- nrow(sh$X_pred12)
    sh$N_rep <- 10000

    posth <- gqmod$generate_quantities(fitted_params = postcsv,
                                     data = sh,
                                      parallel_chains = 4)
    

    # white race/ethnicity
    sw <- subset %>%
    filter(black == 0,
           hispanic == 0,
           asian == 0)
    sw <- generate_X(sw)
    sw$N_pred <- nrow(sw$X_pred12)
    sw$N_rep <- 10000

    postw <- gqmod$generate_quantities(fitted_params = postcsv,
                                     data = sw,
                                      parallel_chains = 4)
    
    
    outcomes <- list( outcomes = data.frame(
        subgroup = rep(c("white","black","asian","hispanic"), n = 9),
        outcome = rep(rep(c("P_imv28","P_mort28","P_ie28"), each = 12)),
        quantity = rep(rep(c("mean","lb95","ub95"), each = 4), n = 3),
        value = c(
            package_gq(postw)$imv$mean[1],
            package_gq(postb)$imv$mean[2],
            package_gq(posta)$imv$mean[3],
            package_gq(posth)$imv$mean[4],

            package_gq(postw)$imv$lb95[1],
            package_gq(postb)$imv$lb95[2],
            package_gq(posta)$imv$lb95[3],
            package_gq(posth)$imv$lb95[4],

            package_gq(postw)$imv$ub95[1],
            package_gq(postb)$imv$ub95[2],
            package_gq(posta)$imv$ub95[3],
            package_gq(posth)$imv$ub95[4],

            package_gq(postw)$te$mean[1],
            package_gq(postb)$te$mean[2],
            package_gq(posta)$te$mean[3],
            package_gq(posth)$te$mean[4],

            package_gq(postw)$te$lb95[1],
            package_gq(postb)$te$lb95[2],
            package_gq(posta)$te$lb95[3],
            package_gq(posth)$te$lb95[4],

            package_gq(postw)$te$ub95[1],
            package_gq(postb)$te$ub95[2],
            package_gq(posta)$te$ub95[3],
            package_gq(posth)$te$ub95[4],

            0,
            package_gq(postb)$ie$mean[1],
            package_gq(posta)$ie$mean[2],
            package_gq(posth)$ie$mean[3],

            0,
            package_gq(postb)$ie$lb95[1],
            package_gq(posta)$ie$lb95[2],
            package_gq(posth)$ie$lb95[3],

            0,
            package_gq(postb)$ie$ub95[1],
            package_gq(posta)$ie$ub95[2],
            package_gq(posth)$ie$ub95[3]
                    )),
        imvdiff = data.frame(
            subgroup = c("Black","Asian","Hispanic"),
            mean = c(as.numeric(filter(package_gq(postb)$imvdiff, subgroup == "Black")[2]),
                    as.numeric(filter(package_gq(posta)$imvdiff, subgroup == "Asian")[2]),
                    as.numeric(filter(package_gq(posth)$imvdiff, subgroup == "Hispanic")[2])),
            lb95 = c(as.numeric(filter(package_gq(postb)$imvdiff, subgroup == "Black")[3]),
                    as.numeric(filter(package_gq(posta)$imvdiff, subgroup == "Asian")[3]),
                    as.numeric(filter(package_gq(posth)$imvdiff, subgroup == "Hispanic")[3])),
            ub95 = c(as.numeric(filter(package_gq(postb)$imvdiff, subgroup == "Black")[4]),
                    as.numeric(filter(package_gq(posta)$imvdiff, subgroup == "Asian")[4]),
                    as.numeric(filter(package_gq(posth)$imvdiff, subgroup == "Hispanic")[4])),
            p_less_imv = c(as.numeric(filter(package_gq(postb)$imvdiff, subgroup == "Black")[5]),
                    as.numeric(filter(package_gq(posta)$imvdiff, subgroup == "Asian")[5]),
                    as.numeric(filter(package_gq(posth)$imvdiff, subgroup == "Hispanic")[5]))),

        p_harm = data.frame(
            subgroup = c("Black","Asian","Hispanic"),
            mean = c(package_gq(postb)$p_harm$prob[1],
                    package_gq(posta)$p_harm$prob[2],
                    package_gq(posth)$p_harm$prob[3])

        )
    )
    
    outcomes
}

forest1 <- function(df){
  
  reference <- data.frame(
    variable = c("white",
                 "G4",
                 "medical_surgical_icu",
                 "R0",
                 "hfnc"
    ),
    mean = 1,
    lb95 = 1,
    ub95 = 1)
  
  
  df <- bind_rows(df, reference)
  
  df <- mutate(df,
               Variable = factor(
                 variable,
                 levels = rev(c(
                   "white",
                   "asian",
                   "black",
                   "hispanic",
                   "age",
                   "gender",
                   "copd", 
                   "cancer", 
                   "chf",
                   "dementia",
                   "G1", 
                   "G2",
                   "G3",
                   "G4",
                   "cardiac_icu",
                   "neuro_trauma_icu",
                   "medical_surgical_icu",
                   "teachingstatus",
                   "R0",
                   "R1",
                   "R2",
                   "R3",
                   "R4",
                   "eicu",
                   "gcs",
                   "heart_rate",
                   "resp_rate",
                   "spo2",
                   "fio2",
                   "pressor",
                   "ra",
                   "np",
                   "fm",
                   "niv",
                   "nrb",
                   "hfnc"
                 )),
                 labels = rev(c("Race/ethnicity: White",
                                "Asian",
                                "Black",
                                "Hispanic",
                                "Age",
                                "Sex: Male",
                                "COPD",
                                "Cancer",
                                "CHF",
                                "Dementia",
                                "Admission year: 2008 - 2010",
                                "2011 - 2013",
                                "2014 - 2016",
                                "2017 - 2019",
                                "ICU type: Cardiac",
                                "Neuro-trauma",
                                "Medical-surgical",
                                "Hospital type: Teaching",
                                "Region: Northeast",
                                "Midwest",
                                "West",
                                "South",
                                "Unknown",
                                "Database: eICU",
                                "Glasgow Coma Scale",
                                "Heart rate",
                                "Respiratory rate",
                                "Peripheral oxygen saturation",
                                "Inspired oxygen fraction",
                                "Vasopressor use",
                                "Oxygen delivery: None",
                                "Nasal prongs",
                                "Face mask",
                                "Non-invasive ventilation",
                                "Non-rebreather mask",
                                "High flow nasal cannula"))))
  
  df <- df %>%
    mutate(label = paste0(
      round(mean, 2), " (",
      round(lb95, 2), " to ",
      round(ub95, 2), ")"
    ))
  
  df$label[-c(1:31)] <- "Reference"
  
  textsize = 3
  
  ggplot(df,
         aes(x = Variable,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
    geom_hline(yintercept = 1, color = "grey80", size = 1) +
    geom_pointrange() +
    theme_minimal() +
    coord_flip() +
    geom_text(aes(x = Variable,
                  label = label, y = 0.25),
              hjust = 0,vjust = 0.5,
              size = textsize) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.5, 0.7, 1, 1.4, 2),
                       limits = c(0.25, 2.5),
                       minor_breaks = NULL) +
    theme(panel.grid.major.y = element_blank()) +
    labs(y = "Posterior hazard ratio",
         x = "")
  
}

forest2 <- function(df){
  
  reference <- data.frame(
    variable = c("white",
                 "G4",
                 "medical_surgical_icu",
                 "R0"
    ),
    mean = 1,
    lb95 = 1,
    ub95 = 1)
  
  
  df <- bind_rows(df, reference)
  
  df <- mutate(df,
               Variable = factor(
                 variable,
                 levels = rev(c(
                   "white",
                   "asian",
                   "black",
                   "hispanic",
                   "age",
                   "gender",
                   "copd", 
                   "cancer", 
                   "chf",
                   "dementia",
                   "G1", 
                   "G2",
                   "G3",
                   "G4",
                   "cardiac_icu",
                   "neuro_trauma_icu",
                   "medical_surgical_icu",
                   "teachingstatus",
                   "R0",
                   "R1",
                   "R2",
                   "R3",
                   "R4",
                   "eicu",
                   "imv_time"
                 )),
                 labels = rev(c("Race/ethnicity: White",
                                "Asian",
                                "Black",
                                "Hispanic",
                                "Age",
                                "Sex: Male",
                                "COPD",
                                "Cancer",
                                "CHF",
                                "Dementia",
                                "Admission year: 2008 - 2010",
                                "2011 - 2013",
                                "2014 - 2016",
                                "2017 - 2019",
                                "ICU type: Cardiac",
                                "Neuro-trauma",
                                "Medical-surgical",
                                "Hospital type: Teaching",
                                "Region: Northeast",
                                "Midwest",
                                "West",
                                "South",
                                "Unknown",
                                "Database: eICU",
                                "Time to invasive ventilation"))))
  
  df <- df %>%
    mutate(label = paste0(
      round(mean, 2), " (",
      round(lb95, 2), " to ",
      round(ub95, 2), ")"
    ))
  
  df$label[-c(1:21)] <- "Reference"
  
  textsize = 3
  
  ggplot(df,
         aes(x = Variable,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
    geom_hline(yintercept = 1, color = "grey80", size = 1) +
    geom_pointrange() +
    theme_minimal() +
    coord_flip() +
    geom_text(aes(x = Variable,
                  label = label, y = 0.25),
              hjust = 0,vjust = 0.5,
              size = textsize) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.5, 0.7, 1, 1.4, 2, 4),
                       limits = c(0.25, 4),
                       minor_breaks = NULL) +
    theme(panel.grid.major.y = element_blank()) +
    labs(y = "Posterior hazard ratio",
         x = "")
  
}

forest3 <- function(df){
  
  reference <- data.frame(
    variable = c("white",
                 "G4",
                 "medical_surgical_icu",
                 "R0"
    ),
    mean = 1,
    lb95 = 1,
    ub95 = 1)
  
  
  df <- bind_rows(df, reference)
  
  df <- mutate(df,
               Variable = factor(
                 variable,
                 levels = rev(c(
                   "white",
                   "asian",
                   "black",
                   "hispanic",
                   "age",
                   "gender",
                   "copd", 
                   "cancer", 
                   "chf",
                   "dementia",
                   "G1", 
                   "G2",
                   "G3",
                   "G4",
                   "cardiac_icu",
                   "neuro_trauma_icu",
                   "medical_surgical_icu",
                   "teachingstatus",
                   "R0",
                   "R1",
                   "R2",
                   "R3",
                   "R4",
                   "eicu",
                   "dc_time"
                 )),
                 labels = rev(c("Race/ethnicity: White",
                                "Asian",
                                "Black",
                                "Hispanic",
                                "Age",
                                "Sex: Male",
                                "COPD",
                                "Cancer",
                                "CHF",
                                "Dementia",
                                "Admission year: 2008 - 2010",
                                "2011 - 2013",
                                "2014 - 2016",
                                "2017 - 2019",
                                "ICU type: Cardiac",
                                "Neuro-trauma",
                                "Medical-surgical",
                                "Hospital type: Teaching",
                                "Region: Northeast",
                                "Midwest",
                                "West",
                                "South",
                                "Unknown",
                                "Database: eICU",
                                "Time to ICU discharge"))))
  
  df <- df %>%
    mutate(label = paste0(
      round(mean, 2), " (",
      round(lb95, 2), " to ",
      round(ub95, 2), ")"
    ))
  
  df$label[-c(1:21)] <- "Reference"
  
  textsize = 3
  
  ggplot(df,
         aes(x = Variable,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
    geom_hline(yintercept = 1, color = "grey80", size = 1) +
    geom_pointrange() +
    theme_minimal() +
    coord_flip() +
    geom_text(aes(x = Variable,
                  label = label, y = 0.25),
              hjust = 0,vjust = 0.5,
              size = textsize) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.5, 0.7, 1, 1.4, 2, 4),
                       limits = c(0.25, 5),
                       minor_breaks = NULL) +
    theme(panel.grid.major.y = element_blank()) +
    labs(y = "Posterior hazard ratio",
         x = "")
  
}

forest_eicu <- function(df){
  
  reference <- data.frame(
    variable = c("white",
                 "medical_surgical_icu",
                 "R0",
                 "hfnc"
    ),
    mean = 1,
    lb95 = 1,
    ub95 = 1)
  
  
  df <- bind_rows(df, reference)
  
  df <- mutate(df,
               Variable = factor(
                 variable,
                 levels = rev(c(
                   "white",
                   "asian",
                   "black",
                   "hispanic",
                   "age",
                   "gender",
                   "copd", 
                   "cancer", 
                   "chf",
                   "dementia",
                   "cardiac_icu",
                   "neuro_trauma_icu",
                   "medical_surgical_icu",
                   "teachingstatus",
                   "R0",
                   "R1",
                   "R2",
                   "R3",
                   "R4",
                   "prop_nonwhite",
                   "gcs",
                   "heart_rate",
                   "resp_rate",
                   "spo2",
                   "fio2",
                   "pressor",
                   "ra",
                   "np",
                   "fm",
                   "niv",
                   "nrb",
                   "hfnc"
                 )),
                 labels = rev(c("Race/ethnicity: White",
                                "Asian",
                                "Black",
                                "Hispanic",
                                "Age",
                                "Sex: Male",
                                "COPD",
                                "Cancer",
                                "CHF",
                                "Dementia",
                                "ICU type: Cardiac",
                                "Neuro-trauma",
                                "Medical-surgical",
                                "Hospital type: Teaching",
                                "Region: Northeast",
                                "Midwest",
                                "West",
                                "South",
                                "Unknown",
                                "Higher proportion of non-white patients",
                                "Glasgow Coma Scale",
                                "Heart rate",
                                "Respiratory rate",
                                "Peripheral oxygen saturation",
                                "Inspired oxygen fraction",
                                "Vasopressor use",
                                "Oxygen delivery: None",
                                "Nasal prongs",
                                "Face mask",
                                "Non-invasive ventilation",
                                "Non-rebreather mask",
                                "High flow nasal cannula"))))
  
  df <- df %>%
    mutate(label = paste0(
      round(mean, 2), " (",
      round(lb95, 2), " to ",
      round(ub95, 2), ")"
    ))
  
  df$label[-c(1:28)] <- "Reference"
  
  textsize = 3
  
  ggplot(df,
         aes(x = Variable,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
    geom_hline(yintercept = 1, color = "grey80", size = 1) +
    geom_pointrange() +
    theme_minimal() +
    coord_flip() +
    geom_text(aes(x = Variable,
                  label = label, y = 0.25),
              hjust = 0,vjust = 0.5,
              size = textsize) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.5, 0.7, 1, 1.4, 2),
                       limits = c(0.25, 2.5),
                       minor_breaks = NULL) +
    theme(panel.grid.major.y = element_blank()) +
    labs(y = "Posterior hazard ratio",
         x = "")
  
}


forest_eicu2 <- function(df){
  
  reference <- data.frame(
    variable = c("white",
                 "medical_surgical_icu",
                 "R0",
                 "hfnc"
    ),
    mean = 1,
    lb95 = 1,
    ub95 = 1)
  
  
  df <- bind_rows(df, reference)
  
  df <- mutate(df,
               Variable = factor(
                 variable,
                 levels = rev(c(
                   "white",
                   "asian",
                   "black",
                   "hispanic",
                   "age",
                   "gender",
                   "copd", 
                   "cancer", 
                   "chf",
                   "dementia",
                   "cardiac_icu",
                   "neuro_trauma_icu",
                   "medical_surgical_icu",
                   "teachingstatus",
                   "R0",
                   "R1",
                   "R2",
                   "R3",
                   "R4",
                   "gcs",
                   "heart_rate",
                   "resp_rate",
                   "spo2",
                   "fio2",
                   "pressor",
                   "ra",
                   "np",
                   "fm",
                   "niv",
                   "nrb",
                   "hfnc"
                 )),
                 labels = rev(c("Race/ethnicity: White",
                                "Asian",
                                "Black",
                                "Hispanic",
                                "Age",
                                "Sex: Male",
                                "COPD",
                                "Cancer",
                                "CHF",
                                "Dementia",
                                "ICU type: Cardiac",
                                "Neuro-trauma",
                                "Medical-surgical",
                                "Hospital type: Teaching",
                                "Region: Northeast",
                                "Midwest",
                                "West",
                                "South",
                                "Unknown",
                                "Glasgow Coma Scale",
                                "Heart rate",
                                "Respiratory rate",
                                "Peripheral oxygen saturation",
                                "Inspired oxygen fraction",
                                "Vasopressor use",
                                "Oxygen delivery: None",
                                "Nasal prongs",
                                "Face mask",
                                "Non-invasive ventilation",
                                "Non-rebreather mask",
                                "High flow nasal cannula"))))
  
  df <- df %>%
    mutate(label = paste0(
      round(mean, 2), " (",
      round(lb95, 2), " to ",
      round(ub95, 2), ")"
    ))
  
  df$label[-c(1:28)] <- "Reference"
  
  textsize = 3
  
  ggplot(df,
         aes(x = Variable,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
    geom_hline(yintercept = 1, color = "grey80", size = 1) +
    geom_pointrange() +
    theme_minimal() +
    coord_flip() +
    geom_text(aes(x = Variable,
                  label = label, y = 0.25),
              hjust = 0,vjust = 0.5,
              size = textsize) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.5, 0.7, 1, 1.4, 2),
                       limits = c(0.25, 2.5),
                       minor_breaks = NULL) +
    theme(panel.grid.major.y = element_blank()) +
    labs(y = "Posterior hazard ratio",
         x = "")
  
}

makelab <- function(x){paste0(round(x[1]*100,1),
                              " (",
                              round(x[2]*100, 1),
                              " to ",
                              round(x[3]*100, 1),
                              ")")}

tableready <- function(outcomes){
  data.frame(
    white = c(makelab(filter(outcomes$outcomes,subgroup == "white", outcome == "P_imv28")$value),
              NA,
              makelab(filter(outcomes$outcomes,subgroup == "white", outcome == "P_mort28")$value),
              NA),
    black = c(makelab(filter(outcomes$outcomes,subgroup == "black", outcome == "P_imv28")$value),
              makelab(-filter(outcomes$imvdiff,subgroup == "Black")[c(2,4,3)]),
              makelab(filter(outcomes$outcomes,subgroup == "black", outcome == "P_mort28")$value),
              makelab(filter(outcomes$outcomes,subgroup == "black", outcome == "P_ie28")$value)),
    hispanic = c(makelab(filter(outcomes$outcomes,subgroup == "hispanic", outcome == "P_imv28")$value),
                 makelab(-filter(outcomes$imvdiff,subgroup == "Hispanic")[c(2,4,3)]),
                 makelab(filter(outcomes$outcomes,subgroup == "hispanic", outcome == "P_mort28")$value),
                 makelab(filter(outcomes$outcomes,subgroup == "hispanic", outcome == "P_ie28")$value)),
    asian = c(makelab(filter(outcomes$outcomes,subgroup == "asian", outcome == "P_imv28")$value),
              makelab(-filter(outcomes$imvdiff,subgroup == "Asian")[c(2,4,3)]),
              makelab(filter(outcomes$outcomes,subgroup == "asian", outcome == "P_mort28")$value),
              makelab(filter(outcomes$outcomes,subgroup == "asian", outcome == "P_ie28")$value)))
}
