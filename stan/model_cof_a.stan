data {
  int<lower=0> N;
  int<lower=0> p;
  int <lower=1> sp;
  int <lower=0> cof;
  int<lower=1,upper=p> plot[N];
  int <lower=0,upper=sp> species [N];
  int <lower=0,upper=cof> cofactor [N];
  vector <lower=0> [N] H;
  vector <lower=0> [N] dbh;
}


parameters {
  vector <lower=0,upper=100> [cof] alpha;
  real <lower=0,upper=100> beta_0;
  // plot random effect
  vector <lower=0> [p] gamma_plot;
  real<lower=0> sigma_plot;
  // species random effect
  vector <lower=0> [sp] gamma_sp;
  real <lower=0> sigma_sp;
  real<lower=0> sigma;
}
 
model {
  real mu [N];
  for (i in 1:N) {
    mu[i] = gamma_sp[species[i]]*gamma_plot[plot[i]]*(alpha[cofactor[i]] * dbh[i])/ 
            (beta_0+dbh[i]);
  }
  gamma_sp~lognormal(0,sigma_sp);
  gamma_plot ~ lognormal(0, sigma_plot);
  // Likelihood part of Bayesian inference
  H~lognormal(log(mu),sigma);
}
