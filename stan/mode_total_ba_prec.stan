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
  vector <lower=0> [N] ba;
  vector <lower=0> [N] precmin;
}


parameters {
  real <lower=20,upper=100> alpha_0;
  vector <lower=20,upper=100> [so] alpha;
  real <lower=0> sigma_a;
  real <lower=0,upper=100> beta_0;
  vector <lower=0,upper=100> [so] beta;
  real <lower=0> sigma_b;
  real <lower=-3,upper=3> beta_ba;
  real <lower=-3,upper=3> beta_precmin;
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
  real beta_i [N];
  for (i in 1:N) {
    beta_i[i]=beta[systori[i]]*pow(ba[i],beta_ba)*pow(precmin[i],beta_precmin);
    mu[i] = gamma_sp[species[i]]*gamma_plot[plot[i]]*(alpha[systori[i]] * dbh[i])/ 
            (beta_i[i]+dbh[i]);
  }
  alpha~lognormal(log(alpha_0),sigma_a);
  beta~lognormal(log(beta_0),sigma_b);
  gamma_sp~lognormal(0,sigma_sp);
  gamma_plot ~ lognormal(0, sigma_plot);
  // Likelihood part of Bayesian inference
  H~lognormal(log(mu),sigma);
}
