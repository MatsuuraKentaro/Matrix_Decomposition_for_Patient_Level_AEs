library(dplyr)
library(rstan)

## read our phi estimated using all data
d_phi_AC <- read.csv(file = 'input/phi.csv', stringsAsFactors=FALSE)
AE_types <- d_phi_AC %>% select(AE_type, AE_type_j) %>% unique() %>% arrange(AE_type_j) %>% .$AE_type

phi_AC <- d_phi_AC %>%   
  filter(Treatment == 'AC') %>% 
  select(Pattern_k, AE_type_j, Phi_kj) %>% 
  tidyr::spread(key = AE_type_j, value = Phi_kj) %>% 
  select(-Pattern_k)

## AE data example (i' = 1 at the end of the second cycle)  
d_AE_eg1 <- read.csv('input/AE_data_example.csv', stringsAsFactors=FALSE) %>% 
  filter(i_prime == 1, Cycle == 2)

K <- nrow(phi_AC)
E <- ncol(phi_AC)
Y <- integer(E)
Y[d_AE_eg1$AE_type_j] <- d_AE_eg1$Count
Cyc <- d_AE_eg1$Cycle[1]

stanmodel <- stan_model(file='model/model_theta.stan')

## Compare 100 results of MLE
compa <- sapply(1:100, function(seedID) {
  fit_tmp <- NULL
  try(fit_tmp <- optimizing(stanmodel, data=list(E=E, K=K, Phi_est=t(phi_AC), Y=Y, Cyc0=Cyc), seed=seedID, iter=10000, as_vector=FALSE))
  if (is.null(fit_tmp)) {
    modes <- rep(NA, K)
  } else {
    modes <- fit_tmp$par$theta
  }
  c(lp=fit_tmp$value, mode=modes)
}) %>% t() %>% as.data.frame()
seed_best <- which.max(compa$lp)

theta <- unlist(compa[seed_best, -1])
d_theta <- data.frame(Pattern_k=1:K, Theta_k=theta)

lambda <- t(phi_AC) %*% theta %>% as.vector()
d_lambda <- data.frame(AE_type=AE_types, Lambda_j=lambda) %>% arrange(desc(lambda))
