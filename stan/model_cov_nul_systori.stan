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
  real <lower=20,upper=60> alpha_0;
  vector <lower=20,upper=60> [so] alpha;
  real <lower=0> sigma_a;
  real <lower=0,upper=80> beta_0;
  vector <lower=0,upper=80> [so] beta;
  real <lower=0> sigma_b;
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
  alpha~lognormal(log(alpha_0),sigma_a);
  beta~lognormal(log(beta_0),sigma_b);
  gamma_sp~lognormal(0,sigma_sp);
  gamma_plot ~ lognormal(0, sigma_plot);
  H~lognormal(log(mu),sigma);
}
