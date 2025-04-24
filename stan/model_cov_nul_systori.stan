data {
  int<lower=0> N;
  int<lower=0> p;
  int <lower=1> sp;
  int <lower=0> so;
  int<lower=1,upper=p> plot[N];
  int <lower=0,upper=sp> species [N];
  int <lower=0,upper=so> systori [N];
  vector <lower=0> [N] H;
  vector <lower=0> [N] dbh;
}


parameters {
  vector <lower=0,upper=100> [so] alpha;
  vector <lower=0,upper=100> [so] beta;
  vector <lower=0> [p] gamma_plot;
  real<lower=0> sigma_plot;
  vector <lower=0> [sp] gamma_sp;
  real <lower=0> sigma_sp;
  real<lower=0> sigma;
}
 
model {
  real mu [N];
  for (i in 1:N) {
    mu[i] = gamma_sp[species[i]]*gamma_plot[plot[i]]*(alpha[systori[i]] * dbh[i])/ 
            (beta[systori[i]]+dbh[i]);
  }
  alpha~normal(40,10);
  beta~normal(40,10);
  gamma_sp~lognormal(0,sigma_sp);
  gamma_plot ~ lognormal(0, sigma_plot);
  H~lognormal(log(mu),sigma);
}

generated quantities {
  vector[N] log_lik; // Log-likelihood for each observation
  for (i in 1:N) {
    real mu_i = gamma_sp[species[i]] * gamma_plot[plot[i]] * (alpha[systori[i]] * dbh[i]) / 
                (beta[systori[i]] + dbh[i]);
    log_lik[i] = lognormal_lpdf(H[i] | log(mu_i), sigma);
  }
}
