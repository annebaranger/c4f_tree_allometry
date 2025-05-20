data {
  int<lower=0> N;
  int<lower=0> p;
  int <lower=1> sp;
  int <lower=0> ncof;
  int<lower=1,upper=p> plot[N];
  int <lower=0,upper=sp> species [N];
  int <lower=0,upper=ncof> cof [N];
  vector <lower=0> [N] H;
  vector <lower=0> [N] dbh;
  vector<lower=0>[N] weight;
}


parameters {
  vector <lower=0,upper=100> [ncof] alpha;
  vector <lower=0,upper=100> [ncof] beta;
  vector <lower=0> [p] gamma_plot;
  real<lower=0> sigma_plot;
  vector <lower=0> [sp] gamma_sp;
  real <lower=0> sigma_sp;
  real<lower=0> sigma;
}
 
model {
  real mu [N];
  for (i in 1:N) {
    mu[i] = gamma_sp[species[i]]*gamma_plot[plot[i]]*(alpha[cof[i]] * dbh[i])/ 
            (beta[cof[i]]+dbh[i]);
  }
  alpha~normal(40,10);
  beta~normal(40,10);
  gamma_sp~lognormal(0,sigma_sp);
  gamma_plot ~ lognormal(0, sigma_plot);

  for (i in 1:N) {
    target += weight[i] * lognormal_lpdf(H[i] | log(mu[i]), sigma);
    }
}
