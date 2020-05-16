/* name: well_titer_means.stan
   author: Dylan Morris <dhmorris@princeton.edu>

   description: Stan code for inferring 
   mean titers across multiple replicates
*/


functions {
}

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 0> n_total_datapoints;
  int<lower = 0> n_titers;
  int<lower = 0> n_replicate_sets;
  int<lower = 0, upper = 1> well_status [n_total_datapoints];
  int dilution [n_total_datapoints];
  
  int<lower = 0, upper = n_titers> titer_id[n_total_datapoints];
    
  int<lower = 0, upper = n_replicate_sets>
    titer_replicate_set_id[n_titers];

  int<lower = 0, upper = 1> debug;

  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real titer_mean_prior_mean;
  real<lower = 0> titer_mean_prior_sd;
  
  real titer_sd_prior_mean;
  real<lower = 0> titer_sd_prior_sd;

}

transformed data {
}


parameters{
  vector[n_titers] titer_errors;
  vector[n_replicate_sets] set_mean_log10_tcid50;
  real<lower = 0> titer_sd;
}

transformed parameters {
  vector[n_titers] log10_tcid50;
  // non centered hierarchical parameterization
  for (t_id in 1:n_titers) {
    int set_id = titer_replicate_set_id[t_id];
    log10_tcid50[t_id] = set_mean_log10_tcid50[set_id] +
      titer_sd * titer_errors[t_id];
  }
}

model {

  // observation process: poisson single hit
  for (obs_id in 1:n_total_datapoints) {
    int t_id = titer_id[obs_id];
    int set_id = titer_replicate_set_id[t_id];
    real ith_titer = log10_tcid50[t_id];
    real dilute_dose = ith_titer + dilution[obs_id];

    if(well_status[obs_id] == 0) {
      target += poisson_lpmf(0 | pow(10, dilute_dose));
      
    } else if (well_status[obs_id] == 1) {
      target += poisson_lccdf(0 | pow(10, dilute_dose));
      
    } else {
      // well status must be negative (0) or positive (1)
      reject("well_status data must be one or zero, found",
             well_status[obs_id]);
    }
  }

  // non centered normal errors
  titer_errors ~ normal(0, 1);
  
  // priors
  set_mean_log10_tcid50 ~ normal(titer_mean_prior_mean,
                                 titer_mean_prior_sd);

  titer_sd ~ normal(titer_sd_prior_mean,
                    titer_sd_prior_sd);

  if(debug) {
    print("transition successful!");
    print("log prob: ", target());
    print("mean titers:   ", set_mean_log10_tcid50);
    print("sd:   ", titer_sd);
  }
}

generated quantities {

}
