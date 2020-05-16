functions {

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
  int<lower = 1> control_index[n_total_datapoints];
  vector<lower = 0>[n_total_datapoints] treatment_end_time;

  int<lower = 1> n_wells;
  int<lower = 1> titer_experiment_id[n_titers];
  int<lower = 1> titer_replicate_id[n_titers];
  vector<lower = 0>[n_titers] titer_treatment_end_time;
  int<lower = 1> titer_control_index[n_titers];
  vector<lower = 0>[n_titers] titer_time;
  
  int first_dilution;
  int<lower = first_dilution> last_dilution;

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

transformed parameters {}

model {

  // observation process: poisson single hit
  real ith_predicted_titer;
  real dilute_dose;

  for (i_dat in 1:n_total_datapoints) {
    int i_exp = experiment_id[i_dat];
    int i_repl = replicate_id[i_dat];
    
    real time = times[i_dat];
    real t_end = treatment_end_time[i_dat];
    
    if(time < t_end)
      ith_predicted_titer = intercept[i_exp, i_repl] -
        decay_rate[i_exp] * time;
    else
      ith_predicted_titer = intercept[i_exp, i_repl] -
        decay_rate[i_exp] * treatment_end_time[i_dat] -
        decay_rate[control_index[i_dat]] *
        (time - treatment_end_time[i_dat]);
    
    dilute_dose = ith_predicted_titer + dilution[i_dat];

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

  /* vector[n_titers] prior_predictive_titers; */
  /* vector[n_titers] posterior_predictive_titers; */
  
  /* matrix[n_experiments, n_replicates] intercept_prior_pred;  */
  
  /* for(i_exp in 1:n_experiments){ */
  /*   for(i_repl in 1:n_replicates){ */
  /*     intercept_prior_pred[i_exp, i_repl] = */
  /*       normal_rng(intercept_prior_mean, */
  /*                  intercept_prior_sd); */
  /*   } */
  /*     } */
  
  /* for (i_titer in 1:n_titers) { */
  /*   int i_exp = titer_experiment_id[i_titer]; */
  /*   int i_repl = titer_replicate_id[i_titer]; */
  /*   real t_time = titer_time[i_titer]; */
  /*   real tt_end = titer_treatment_end_time[i_titer]; */
  /*   real ith_post_titer; */
  /*   real ith_prior_titer; */
    
  /*   if(t_time < titer_treatment_end_time[i_titer]) */
  /*     ith_post_titer = intercept[i_exp, i_repl] - */
  /*       decay_rate[i_exp] * t_time; */
  /*   else */
  /*     ith_post_titer = intercept[i_exp, i_repl] - */
  /*       decay_rate[i_exp] *  tt_end - */
  /*       decay_rate[titer_control_index[i_titer]] * */
  /*       (t_time - tt_end); */
    
  /*   ith_prior_titer = intercept_prior_pred[i_exp, i_repl] - */
  /*     decay_prior_pred[i_exp] * t_time; */

  /*   posterior_predictive_titers[i_titer] = */
  /*     obs_titer_rng(ith_post_titer, */
  /*                   n_wells, */
  /*                   first_dilution, */
  /*                   last_dilution, */
  /*                   debug); */
  /*   prior_predictive_titers[i_titer] = */
  /*     obs_titer_rng(ith_prior_titer, */
  /*                   n_wells, */
  /*                   first_dilution, */
  /*                   last_dilution, */
  /*                   debug); */
  /* } */
}
