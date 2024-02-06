data {
  int <lower=0> N;
  vector [N] H ;
  vector [N] dbh ;
}

parameters {real <lower=0> alpha;real <lower=0> beta;real <lower=0> sigma;}
model {vector [N] height;height = alpha*(1-exp(-dbh/beta)) ;
H ~ lognormal(log(height), sigma);
}