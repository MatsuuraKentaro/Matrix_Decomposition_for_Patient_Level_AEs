data {
  int<lower=1> E;       // num of AE types
  int<lower=1> K;       // num of patterns for treatments
  matrix[E,K] Phi_est;  // estimated phi
  int<lower=0> Y[E];    // AE data vector
  real Cyc0;            // num of cycles
}

transformed data {
  vector[E] Cyc = rep_vector(Cyc0, E);
}

parameters {
  vector<lower=0>[K] theta;
}

model {
  vector[E] mu = Cyc .* (Phi_est * theta);
  theta ~ exponential(1);
  Y  ~ poisson(mu);
}
