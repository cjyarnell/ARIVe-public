// ARIVe project
// Stan code for Bayesian multistate model 

functions {

    real partial_sum_exponentiallcdf(array[] real y_slice,
                   int start, int end,
                   matrix x,
                   vector beta, 
                   int P) {
    return exponential_lcdf(y_slice | exp(x[start:end,] * beta));
  }
  
    real partial_sum_exponentiallccdf(array[] real y_slice,
                   int start, int end,
                   matrix x,
                   vector beta, 
                   int P) {
    return exponential_lccdf(y_slice | exp(x[start:end,] * beta));
  }

    real partial_sum_weibulllcdf(array[] real y_slice,
                   int start, int end,
                   matrix x,
                   vector beta,
                   real logalpha,
                   int P) {
           real alpha = exp(logalpha);
    return weibull_lcdf(y_slice | alpha, exp((x[start:end,] * beta)*alpha));
  }
  
    real partial_sum_weibulllccdf(array[] real y_slice,
                   int start, int end,
                   matrix x,
                   vector beta,
                   real logalpha,
                   int P) {        
           real alpha = exp(logalpha);
    return weibull_lccdf(y_slice | alpha, exp((x[start:end,] * beta)*alpha));
  }


}

data {

/*
   states as follows:
   1 = O2
   2 = IMV
   3 = ICU DC
   4 = Death
  
   Transition matrix Q
    (*, 1, 1, 1)
    (0, *, 1, 1)
    (0, 0, *, 1)
    (0, 0, 0, *)
  
   so we include the following 6 cause-specific models:
   
   1-2, 1-3, 1-4 (timevarying covariates, time unit = days)
   2-3, 2-4, 3-4 (time-to-transition, baseline, time unit = 28 days)
  
  Stan doesn't accept "lists" where each element has different lengths
  so for 6 separate models we will just code each separately for clarity
  using the state numbers to denote the model
  ie 12 = model for state 1 to state 2 transition
  
*/

// State 1 - State 2
  int<lower=0> P12; // number of beta parameters including intercept

  int<lower=0> N_cen12; // total number of censored observations
  int<lower=0> N_obs12; // total number of observed transitions
  
  // covariates 
  matrix[N_cen12, P12] X_cen12; 
  matrix[N_obs12, P12] X_obs12; 
  
  // observed times
  array[N_cen12] real y_cen12;
  array[N_obs12] real y_obs12;

// State 1 - State 3
  
  int<lower=0> P13; // number of beta parameters including intercept

  int<lower=0> N_cen13; // total number of censored observations
  int<lower=0> N_obs13; // total number of observed transitions
  
  // covariates 
  matrix[N_cen13, P13] X_cen13; 
  matrix[N_obs13, P13] X_obs13; 
  
  // observed times
  array[N_cen13] real y_cen13;
  array[N_obs13] real y_obs13;

// State 1 - State 4
  
    int<lower=0> P14; // number of beta parameters including intercept

  int<lower=0> N_cen14; // total number of censored observations
  int<lower=0> N_obs14; // total number of observed transitions
  
  // covariates 
  matrix[N_cen14, P14] X_cen14; 
  matrix[N_obs14, P14] X_obs14; 
  
  // observed times
  array[N_cen14] real y_cen14;
  array[N_obs14] real y_obs14;

// State 2 - State 3
  
    int<lower=0> P23; // number of beta parameters including intercept

  int<lower=0> N_cen23; // total number of censored observations
  int<lower=0> N_obs23; // total number of observed transitions
  
  // covariates 
  matrix[N_cen23, P23] X_cen23; 
  matrix[N_obs23, P23] X_obs23; 
  
  // observed times
  array[N_cen23] real y_cen23;
  array[N_obs23] real y_obs23;

// State 2 - State 4
  
    int<lower=0> P24; // number of beta parameters including intercept

  int<lower=0> N_cen24; // total number of censored observations
  int<lower=0> N_obs24; // total number of observed transitions
  
  // covariates 
  matrix[N_cen24, P24] X_cen24; 
  matrix[N_obs24, P24] X_obs24; 
  
  // observed times
  array[N_cen24] real y_cen24;
  array[N_obs24] real y_obs24;

// State 3 - State 4
  
    int<lower=0> P34; // number of beta parameters including intercept

  int<lower=0> N_cen34; // total number of censored observations
  int<lower=0> N_obs34; // total number of observed transitions
  
  // covariates 
  matrix[N_cen34, P34] X_cen34; 
  matrix[N_obs34, P34] X_obs34; 
  
  // observed times
  array[N_cen34] real y_cen34;
  array[N_obs34] real y_obs34;
  
/*  // for generated quantities
  
  // for causal mediation
  
  int N_rep; // number of reps for the causal mediation
  
  vector[P12] X_pred12;
  vector[P13] X_pred13;
  vector[P14] X_pred14;
  vector[P23] X_pred23;
  vector[P24] X_pred24;
  vector[P34] X_pred34;
  
  int N_pred; // predictions per iteration
  
  // for state-by-time plot
  
  int N_times; // number of times
  vector[N_times] times;
*/
}


parameters {
    vector[P12] beta12;
    
    vector[P13] beta13;

    vector[P14] beta14;

    vector[P23] beta23;
    real logalpha23;

    vector[P24] beta24;
    real logalpha24;

    vector[P34] beta34;
    real logalpha34;

}

transformed parameters{
  
  // model Weibull rate as function of covariates
  vector[N_cen23] lambda_cen23;
  vector[N_obs23] lambda_obs23;
  
  real alpha23 = exp(logalpha23);
  
  // standard weibull AFT re-parameterization
  lambda_cen23 = exp((X_cen23*beta23)*alpha23);
  lambda_obs23 = exp((X_obs23*beta23)*alpha23);
  
  vector[N_cen24] lambda_cen24;
  vector[N_obs24] lambda_obs24;
  
  real alpha24 = exp(logalpha24);
  
  // standard weibull AFT re-parameterization
  lambda_cen24 = exp((X_cen24*beta24)*alpha24);
  lambda_obs24 = exp((X_obs24*beta24)*alpha24);
  
  vector[N_cen34] lambda_cen34;
  vector[N_obs34] lambda_obs34;
  
  real alpha34 = exp(logalpha34);
  
  // standard weibull AFT re-parameterization
  lambda_cen34 = exp((X_cen34*beta34)*alpha34);
  lambda_obs34 = exp((X_obs34*beta34)*alpha34);
}


model {
  // priors
  beta12 ~ normal(0,0.1);
  
  beta13 ~ normal(0,0.1);

  beta14 ~ normal(0,0.1);

  beta23 ~ normal(0,0.3);
  logalpha23~ normal(0,0.4);

  beta24 ~ normal(0,0.3);
  logalpha24~ normal(0,0.4);

  beta34 ~ normal(0,0.3);
  logalpha34~ normal(0,0.4);
  
  int grainsize = 1;

  // likelihood function

  target += reduce_sum(partial_sum_exponentiallcdf, 
                       y_obs12, grainsize, X_obs12, beta12, P12);
  target += reduce_sum(partial_sum_exponentiallccdf, 
                       y_cen12, grainsize, X_cen12, beta12, P12);

  target += reduce_sum(partial_sum_exponentiallcdf, 
                       y_obs13, grainsize, X_obs13, beta13, P13);
  target += reduce_sum(partial_sum_exponentiallccdf, 
                       y_cen13, grainsize, X_cen13, beta13, P13);

  target += reduce_sum(partial_sum_exponentiallcdf, 
                       y_obs14, grainsize, X_obs14, beta14, P14);
  target += reduce_sum(partial_sum_exponentiallccdf, 
                       y_cen14, grainsize, X_cen14, beta14, P14);

  target += reduce_sum(partial_sum_weibulllcdf, 
                       y_obs23, grainsize, X_obs23, beta23, logalpha23, P23);
  target += reduce_sum(partial_sum_weibulllccdf, 
                       y_cen23, grainsize, X_cen23, beta23, logalpha23, P23);

  target += reduce_sum(partial_sum_weibulllcdf, 
                       y_obs24, grainsize, X_obs24, beta24, logalpha24, P24);
  target += reduce_sum(partial_sum_weibulllccdf, 
                       y_cen24, grainsize, X_cen24, beta24, logalpha24, P24);

  target += reduce_sum(partial_sum_weibulllcdf, 
                       y_obs34, grainsize, X_obs34, beta34, logalpha34, P34);
  target += reduce_sum(partial_sum_weibulllccdf, 
                       y_cen34, grainsize, X_cen34, beta34, logalpha34, P34);

}

generated quantities {}

