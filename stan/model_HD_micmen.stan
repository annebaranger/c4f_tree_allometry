data {
  int <lower=0> N;
  vector [N] H ;
  vector [N] dbh ;
}

parameters {real <lower=0> alpha;real <lower=0> beta;real <lower=0> sigma;}
model {vector [N] height;height = (alpha*dbh)./((1/beta)+dbh) ;
H ~ lognormal(log(height), sigma);
}