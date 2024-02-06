data {
  int <lower=0> N;
  vector [N] H ;
  vector [N] dbh ;
}

parameters {
real  alpha;
real <lower=0> beta;
real <lower=0> sigma;
}

model {
vector [N] height;
height = exp(alpha+beta*log(dbh)) ;
H ~ lognormal(log(height), sigma);
}
