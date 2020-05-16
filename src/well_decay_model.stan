functions {

  int well_rng(real true_log10_titer,
               real dilution) {
    int hit = 0;
    real virions = log(2) * 10^(true_log10_titer - dilution);
    if(virions > 1e9){
      hit = 1;
    } else if(virions > 0) {
      hit = poisson_rng(virions) > 0;
    } else {
      hit = 0;
    }
    return hit;    
  }

  real obs_titer_rng(real true_log10_titer,
                     int n_wells,
                     int first_dilution,
                     int last_dilution,
                     int debug) {

    real result = -0.5 + first_dilution; // RML style titer
    int hit = 0;
    real virions = 0;
    
    for (dil in first_dilution:last_dilution) {
      for(well in 1:n_wells){
        virions = log(2) * 10^(true_log10_titer - dil);
        if(virions > 1e9){
          hit = 1;
        } else if(virions > 0) {
          hit = poisson_rng(virions) > 0;
        } else {
          hit = 0;
        }
        result += hit / (n_wells * 1.0);

        if (debug) {
          print("true titer: ", true_log10_titer);
          print("dilution: ", dil);
          print("virions: ", virions);
          print("hit: ", hit);
        }
      }
    }

    if (debug) {
      print("measured titer:", result);
    }

    return (result);
  }

}

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 1> n_total_datapoints;
  int<lower = 1> n_experiments;
  int<lower = 1> n_replicates;
  int<lower = 1> n_titers;
  int<lower = 0, upper = 1> well_status [n_total_datapoints];
  int dilution [n_total_datapoints];
  
  vector<lower = 0>[n_total_datapoints] times;
  int<lower = 1> replicate_id[n_total_datapoints];
  int<lower = 1> experiment_id[n_total_datapoints];

  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real intercept_prior_mean;
  real<lower = 0> intercept_prior_sd;

  real decay_rate_prior_mean;
  real<lower = 0> decay_rate_prior_sd;
  
  real lower_lim_decay_rate;
  
  int<lower = 0, upper = 1> debug;

}

transformed data {
  real lld;
  lld = lower_lim_decay_rate;
}


parameters{
  matrix[n_experiments, n_replicates] intercept;
  vector<lower = lld>[n_experiments] decay_rate;
  vector<lower = lld>[n_experiments] decay_prior_pred; 
}

transformed parameters {
}

model {

  // observation process: poisson single hit
  for (i_dat in 1:n_total_datapoints) {
    int i_exp = experiment_id[i_dat];
    int i_repl = replicate_id[i_dat];

    real ith_predicted_titer = intercept[i_exp, i_repl] -
      decay_rate[i_exp] * times[i_dat];
    real dilute_dose = ith_predicted_titer + dilution[i_dat];

    if(well_status[i_dat] == 0) {
      target += poisson_lpmf(0 | log(2) * pow(10, dilute_dose));

      // use debug to check log prob getting properly incremented
      if(debug) {
        print("dose: ", dilute_dose);
        print("lpmf: ", poisson_lpmf(0 | log(2) * pow(10, dilute_dose)));
      }
      
    } else if (well_status[i_dat] == 1) {
      target += poisson_lccdf(0 | log(2) * pow(10, dilute_dose));

      // use debug to check log prob getting properly incremented
      if(debug) {
        print("dose: ", dilute_dose);
        print("lccdf: ", poisson_lccdf(0 | log(2) * pow(10, dilute_dose)));
      }
      
    } else {
      // well status must be negative (0) or positive (1)
      reject("well_status data must be one or zero, found",
             well_status[i_dat]);
      
    }
  }

  // priors
  for(i_exp in 1:n_experiments)
    intercept[i_exp] ~ normal(intercept_prior_mean,
                              intercept_prior_sd);
  
  decay_rate ~ normal(decay_rate_prior_mean,
                      decay_rate_prior_sd);
  decay_prior_pred ~ normal(decay_rate_prior_mean,
                            decay_rate_prior_sd);

}

generated quantities {

  vector[n_total_datapoints] prior_predictive_wells;
  vector[n_total_datapoints] posterior_predictive_wells;
  
  matrix[n_experiments, n_replicates] intercept_prior_pred; 
  
  for(i_exp in 1:n_experiments){
    for(i_repl in 1:n_replicates){
      intercept_prior_pred[i_exp, i_repl] =
        normal_rng(intercept_prior_mean,
                   intercept_prior_sd);
    }
      }
  
  for (i_well in 1:n_total_datapoints) {
    int i_exp = experiment_id[i_well];
    int i_repl = replicate_id[i_well];
    real well_time = times[i_well];
    real dil = dilution[i_well];
    
    real ith_post_titer = intercept[i_exp, i_repl] -
      decay_rate[i_exp] * well_time;
    
    real ith_prior_titer = intercept_prior_pred[i_exp, i_repl] -
      decay_prior_pred[i_exp] * well_time;

    posterior_predictive_wells[i_well] =
      well_rng(ith_post_titer,
               dil);
    
    prior_predictive_wells[i_well] =
      well_rng(ith_prior_titer,
               dil);
  }
}
