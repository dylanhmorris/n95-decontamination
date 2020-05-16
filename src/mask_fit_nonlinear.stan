/* name: mask_fit_nonlinear.stan
   author: Dylan Morris <dhmorris@princeton.edu>

   description: Stan code for inferring
   loss of mask integrity with a nonlinear
   (spline-like) model.
*/

functions {

  real integrity_rng(real mean_logit_obs,
                     real sd_logit_obs){
    return inv_logit(normal_rng(mean_logit_obs,
                                sd_logit_obs));
  }
  
  real obs_fit_factor_rng(real mean_logit_obs,
                          real sd_logit_obs) {
    real ff_rounded = 0;
    
    real ff_raw = 1 / (1 - integrity_rng(mean_logit_obs,
                                         sd_logit_obs));
    if(ff_raw > 200)
      ff_rounded = 200;
    else
      ff_rounded = round(ff_raw);
    return ff_rounded; 
  }
}

data {
  
  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 1> n_subjects;
  int<lower = 1> n_treatments;
  int<lower = 1> n_masks;
  int<lower = 1> max_n_decons;
  int<lower = 1> n_datapoints;
  int<lower = 0, upper = n_datapoints> n_observations;

  int<lower = 1, upper = n_subjects> subject_id[n_datapoints];
  int<lower = 1, upper = n_treatments> treatment_id[n_datapoints];
  int<lower = 1, upper = n_masks> mask_id[n_datapoints];
  int<lower = 1, upper = n_treatments> mask_treatment_id[n_masks];
  int<lower = 0, upper = max_n_decons> decon_number[n_datapoints];

  // parameters for censoring
  real upper_ff_limit;
  real ff_interval;

  // observations
  vector<lower = 1, upper = upper_ff_limit>[n_observations] obs_fit;

  // whether to print evaluations for debugging
  int debug; 
  
  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  // in logit space
  real fit_decay_mean;
  real<lower = 0> fit_decay_sd;
  
  // in logit space
  real<lower = 0> treatment_effect_mean;
  real<lower = 0> treatment_effect_sd;

  // in probability space
  real<lower = 0, upper = 1> mask_intercept_filt_mean_mean;
  real<lower = 0> mask_intercept_filt_mean_sd;

  // in logit space
  real<lower = 0> logit_mask_intercept_filt_sd_mean;
  real<lower = 0> logit_mask_intercept_filt_sd_sd;

  // on (0, \infty]

  real<lower = 0> sd_log_decay_mean;
  real<lower = 0> sd_log_decay_sd;
  real<lower = 0> sd_logit_obs_mean;
  real<lower = 0> sd_logit_obs_sd;
  
  // all these should be less than or equal to 0.5
  // to prevent squashing of filtrations
  // to 0 and 1 (due to the nonlinearity
  // of logit space. 0.5 is a good rough
  // choice for a weakly informative SD
  // prior

}

transformed data {

  real upper_filt_limit = 1 - 1.0 / upper_ff_limit;
}


parameters {

  // hierarchical decay rates (in logit space)
  matrix<lower = 0>[n_treatments, max_n_decons] treatment_effect;

  // intercept (in probability space),
  // for the filtration = 1 - (1 / fit factor)
  real<lower = 0, upper = 1> mask_intercept_filt_mean;
  // sd about the logit of that mean, in logit space
  real<lower = 0> logit_mask_intercept_filt_sd;
  
  vector[n_masks] mask_intercept_errors;
  matrix[n_masks, max_n_decons] log_decay_errors;

  // error terms for filtration, in logit space
  matrix<lower = 0>[n_treatments, max_n_decons] sd_log_decay;
  real<lower = 0> sd_logit_obs;

}

transformed parameters {
  real logit_mask_intercept_filt_mean;
  vector[n_masks] mask_intercept_filt;
  vector[n_masks] logit_mask_intercept_filt; // filtration intercepts in logit space

  // underlying value of p = 1 - 1/ff
  matrix[n_masks, max_n_decons] logit_underlying_integrity;
  matrix<lower = 0, upper = 1>[n_masks, max_n_decons]
    underlying_integrity;
  
  // calculate intercepts in logit space
  logit_mask_intercept_filt_mean = logit(mask_intercept_filt_mean);
  logit_mask_intercept_filt = logit_mask_intercept_filt_mean +
    logit_mask_intercept_filt_sd * mask_intercept_errors;
  
  // calculate values after decontamination
  for(i_mask in 1:n_masks){
    int mt_id = mask_treatment_id[i_mask];

    for(decon_no in 1:max_n_decons){
      real previous_level = 0;
      real decay = 0;
      real true_decay = 0;
      // initialize first entry by decaying from intercept;
      // subsequent entries decay from previous entry
      if(decon_no < 2)
        previous_level = logit_mask_intercept_filt[i_mask];
      else
        previous_level = logit_underlying_integrity[i_mask, decon_no - 1];

      // decay rate depends on which decontamination it is
      decay = treatment_effect[mt_id, decon_no];

      if(mt_id > 1)
        decay += treatment_effect[1, decon_no];

      true_decay = exp(log(decay) + sd_log_decay[mt_id, decon_no] *
                       log_decay_errors[i_mask, decon_no]);
      
      logit_underlying_integrity[i_mask, decon_no] =
        previous_level - true_decay;
    }
  }

  // calculate probability-scale values
  underlying_integrity = inv_logit(logit_underlying_integrity);
  mask_intercept_filt = inv_logit(logit_mask_intercept_filt);
  
} // close transformed parameter block

model {

  // priors for params

  // null decay rate
  treatment_effect[1] ~ normal(fit_decay_mean,
                               fit_decay_sd);

  // non null decay rates
  for(i_treat in 2:n_treatments){
    treatment_effect[i_treat] ~ normal(treatment_effect_mean,
                                       treatment_effect_sd);
  }
  mask_intercept_filt_mean ~ normal(mask_intercept_filt_mean_mean,
                                    mask_intercept_filt_mean_sd);
  
  logit_mask_intercept_filt_sd ~ normal(logit_mask_intercept_filt_sd_mean,
                                        logit_mask_intercept_filt_sd_sd);

  // non centered intercept errors in logit space
  mask_intercept_errors ~ normal(0, 1);
  for(i_mask in 1:n_masks){
    log_decay_errors[i_mask] ~ normal(0, 1);
  }

  for(t_id in 1:n_treatments)
    sd_log_decay[t_id] ~ normal(sd_log_decay_mean,
                                sd_log_decay_sd);
  
  sd_logit_obs ~ normal(sd_logit_obs_mean,
                        sd_logit_obs_sd);

  if(debug){
        print("null effect: ", treatment_effect[1]);
        print("treatment_effect: ", treatment_effect);
        print("mask_intercept_filt: ", mask_intercept_filt);
        print("total log prob: ", target());
  }
  
  // observation process: interval censoring
  for (i_obs in 1:n_observations) {

    int decon_no = decon_number[i_obs];
    int i_mask = mask_id[i_obs];
    real obs = obs_fit[i_obs];
    real mu = 0; // declare; define below
    real sigma = sd_logit_obs;

    real upper_bound = logit(1.0 - 1.0 / (obs + ff_interval / 2));
    real lower_bound = logit(1.0 - 1.0 / (obs - ff_interval / 2));

    // calculate what the mean should be
    if(decon_no < 1){
      mu = logit_mask_intercept_filt[i_mask];
    }
    else {
      mu = logit_underlying_integrity[i_mask, decon_no];
    }

    // now calculate log prob
    if (obs >= upper_ff_limit) {
      
      target += normal_lccdf(logit(upper_filt_limit) |  mu, sigma);

      if(debug){
        print("mean: ", mu);
        print("sd: ", sigma);
        print("observed logit filtration: ", logit(upper_filt_limit));
        print("log likelihood: ", normal_lccdf(logit(upper_filt_limit) |  mu, sigma));
        print("total log prob: ", target());
      }

    } else if (upper_bound - lower_bound > 0) {
      
      target +=
        log_diff_exp(normal_lcdf(upper_bound | mu, sigma),
                     normal_lcdf(lower_bound | mu, sigma));
    } else {
            
      target += normal_lpdf(logit(1.0 - 1.0/obs) | mu, sigma);

      if(debug){
        print("mean", mu);
        print("sd", sigma);
        print("observed logit filtration", logit(1.0 - 1.0/obs));
        print("likelihood", normal_lpdf(logit(1.0 - 1.0/obs) | mu, sigma));
        print("total log prob: ", target());
      }
    } // close if
  } // close loop over observations

  if(debug){
    print("log prob :", target());
  }
  
} // close model block

generated quantities {

  // predictive checks

  // posterior
  vector<lower = 0, upper = 1>[n_datapoints] post_pred_obs_integrity;

  // predicted mask distributions
  // at each point for each treatment
  matrix<lower = 0, upper = 1>[n_treatments, max_n_decons]
    predicted_integrity;

  matrix[n_treatments, max_n_decons] logit_predicted_integrity;
  vector<lower = 0, upper = 1>[n_treatments] predicted_intercept;
  vector[n_treatments] logit_predicted_intercept;

  
  for(t_id in 1:n_treatments){
    real previous_level = 0;
    real decay = 0;
    real true_decay = 0;

    previous_level = normal_rng(logit_mask_intercept_filt_mean,
                                logit_mask_intercept_filt_sd);

    logit_predicted_intercept[t_id] = previous_level;

    for(decon_no in 1:max_n_decons){
      decay = treatment_effect[t_id, decon_no];

      if(t_id > 1)
        decay += treatment_effect[1, decon_no];

      // first decay is rel to pristine condition
      // others relative to previous
      if(decon_no > 1)
        previous_level = logit_predicted_integrity[t_id, decon_no - 1];

      true_decay = exp(normal_rng(log(decay),
                                  sd_log_decay[t_id, decon_no]));
      logit_predicted_integrity[t_id, decon_no] =
        previous_level - true_decay;
    }
  }

  // convert integrities back to probability space
  predicted_integrity = inv_logit(logit_predicted_integrity);
  predicted_intercept = inv_logit(logit_predicted_intercept);

  // predictive checks
  for (i_obs in 1:n_datapoints) {
    
    int decon_no = decon_number[i_obs];
    int i_mask = mask_id[i_obs];
    real mu = 0; // declare; define below
    
    // calculate what the mean should be
    if(decon_no < 1){
      mu = logit_mask_intercept_filt[i_mask];
    }
      else {
        mu = logit_underlying_integrity[i_mask, decon_no];
      }

    post_pred_obs_integrity[i_obs] =
      integrity_rng(mu, sd_logit_obs);
  }
}
