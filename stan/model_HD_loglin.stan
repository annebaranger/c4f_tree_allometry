data {
  int <lower=0> N;
  vector [N] H ;
  vector [N] dbh ;
}

parameters {

H ~ lognormal(log(height), sigma);
}