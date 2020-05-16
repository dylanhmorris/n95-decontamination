functions {
}

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 0> n_total_datapoints;
  int<lower = 0> n_titers;
  int<lower = 0, upper = 1> well_status [n_total_datapoints];
  int dilution [n_total_datapoints];
  
  int<lower = 0> titer_id[n_total_datapoints];

  int<lower = 0, upper = 1> debug;

  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real titer_prior_mean;
  real<lower = 0> titer_prior_sd;

}

transformed data {
}


parameters{
  real log10_tcid50[n_titers];
}

transformed parameters {
}

model {

  // observation process: poisson single hit
  for (i_dat in 1:n_total_datapoints) {
    real ith_titer = log10_tcid50[titer_id[i_dat]];
    real dilute_dose = ith_titer + dilution[i_dat];

    if(well_status[i_dat] == 0) {
      target += poisson_lpmf(0 | pow(10, dilute_dose));

      // use debug to check log prob getting properly incremented
      if(debug) {
        print("dose: ", dilute_dose);
        print("lpmf: ", poisson_lpmf(0 | pow(10, dilute_dose)));
      }
      
    } else if (well_status[i_dat] == 1) {
      target += poisson_lccdf(0 | pow(10, dilute_dose));

      // use debug to check log prob getting properly incremented
      if(debug) {
        print("dose: ", dilute_dose);
        print("lccdf: ", poisson_lccdf(0 | pow(10, dilute_dose)));
      }
      
    } else {
      // well status must be negative (0) or positive (1)
      reject("well_status data must be one or zero, found",
             well_status[i_dat]);
      
    }
  }

  // priors
  log10_tcid50 ~ normal(titer_prior_mean,
                        titer_prior_sd);
  
}

generated quantities {

}
