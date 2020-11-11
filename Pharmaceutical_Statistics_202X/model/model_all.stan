data {
  int<lower=1> E;       // num of AE types
  int<lower=1> E_GX;    // num of AE types that have >=2 grades
  int<lower=1> E_GX0;   // num of AE types that have the lowest grade
  int<lower=1> K;       // num of patterns for treatments
  int<lower=1> L;       // num of patterns for study
  int<lower=1, upper=E> GX[E_GX, 2];  // array of (start, end) indices of AE types that have the same AE name
  int<lower=1, upper=E> GX0[E_GX0];   // array of index of AE types that have the lowest grade
  int<lower=1> N_A1;      // num of patients of Arm 1, baseline
  int<lower=1> N_A1_AC;   // num of patients of Arm 1, first treatment
  int<lower=1> N_A1_T;    // num of patients of Arm 1, second treatment
  int<lower=1> N_A2_FAC;  // num of patients of Arm 2
  int<lower=1> N_A3;      // num of patients of Arm 3, baseline
  int<lower=1> N_A3_A;    // num of patients of Arm 3, first treatment
  int<lower=1> N_A3_CMF;  // num of patients of Arm 3, second treatment
  int<lower=1> N_A4;      // num of patients of Arm 4, baseline
  int<lower=1> N_A4_AC;   // num of patients of Arm 4, first treatment
  int<lower=1> N_A4_CMF;  // num of patients of Arm 4, second treatment
  int Nx_A1_Re[N_A1_T];     // indices of patients of Arm 1, second treatment
  int Nx_A3_Re[N_A3_CMF];   // indices of patients of Arm 3, second treatment
  int Nx_A4_Re[N_A4_CMF];   // indices of patients of Arm 4, second treatment
  int<lower=0> Y_A1_S1 [E, N_A1];     // AE data matrix of Arm 1, baseline
  int<lower=0> Y_A1_AC [E, N_A1_AC];  // AE data matrix of Arm 1, first treatment
  int<lower=0> Y_A1_T  [E, N_A1_T];   // AE data matrix of Arm 1, second treatment
  int<lower=0> Y_A2_FAC[E, N_A2_FAC]; // AE data matrix of Arm 2
  int<lower=0> Y_A3_S3 [E, N_A3];     // AE data matrix of Arm 3, baseline
  int<lower=0> Y_A3_A  [E, N_A3_A];   // AE data matrix of Arm 3, first treatment
  int<lower=0> Y_A3_CMF[E, N_A3_CMF]; // AE data matrix of Arm 3, second treatment
  int<lower=0> Y_A4_S3 [E, N_A4];     // AE data matrix of Arm 4, baseline
  int<lower=0> Y_A4_AC [E, N_A4_AC];  // AE data matrix of Arm 4, first treatment
  int<lower=0> Y_A4_CMF[E, N_A4_CMF]; // AE data matrix of Arm 4, second treatment
  vector<lower=0>[N_A1_AC]  Cyc_A1_AC;  // num of cycles of Arm 1, first treatment
  vector<lower=0>[N_A1_T]   Cyc_A1_T;   // num of cycles of Arm 1, second treatment
  vector<lower=0>[N_A2_FAC] Cyc_A2_FAC; // num of cycles of Arm 2
  vector<lower=0>[N_A3_A]   Cyc_A3_A;   // num of cycles of Arm 3, first treatment
  vector<lower=0>[N_A3_CMF] Cyc_A3_CMF; // num of cycles of Arm 3, second treatment
  vector<lower=0>[N_A4_AC]  Cyc_A4_AC;  // num of cycles of Arm 4, first treatment
  vector<lower=0>[N_A4_CMF] Cyc_A4_CMF; // num of cycles of Arm 4, second treatment
  real<lower=0> Hyp1;  // Hyperparameter 1
  real<lower=0> Hyp2;  // Hyperparameter 2
}

transformed data {
  matrix[N_A1_AC,  L+K] cyc_A1_AC  = rep_matrix(Cyc_A1_AC,  L+K);
  matrix[N_A1_T,   L+K] cyc_A1_T   = rep_matrix(Cyc_A1_T,   L+K);
  matrix[N_A2_FAC,   K] cyc_A2_FAC = rep_matrix(Cyc_A2_FAC,   K); // because of no baseline data
  matrix[N_A3_A,   L+K] cyc_A3_A   = rep_matrix(Cyc_A3_A,   L+K);
  matrix[N_A3_CMF, L+K] cyc_A3_CMF = rep_matrix(Cyc_A3_CMF, L+K);
  matrix[N_A4_AC,  L+K] cyc_A4_AC  = rep_matrix(Cyc_A4_AC,  L+K);
  matrix[N_A4_CMF, L+K] cyc_A4_CMF = rep_matrix(Cyc_A4_CMF, L+K);
}

parameters {
  matrix<lower=0>[N_A1,     L] th_A1_S1;
  matrix<lower=0>[N_A1_AC,  K] th_A1_AC;
  matrix<lower=0>[N_A1_T,   K] th_A1_T;
  matrix<lower=0>[N_A2_FAC, K] th_A2_FAC;
  matrix<lower=0>[N_A3,     L] th_A3_S3;
  matrix<lower=0>[N_A3_A,   K] th_A3_A;
  matrix<lower=0>[N_A3_CMF, K] th_A3_CMF;
  matrix<lower=0>[N_A4,     L] th_A4_S3;
  matrix<lower=0>[N_A4_AC,  K] th_A4_AC;
  matrix<lower=0>[N_A4_CMF, K] th_A4_CMF;
  matrix[E,L] phi_S1_raw;
  matrix[E,L] phi_S3_raw;
  matrix[E,K] phi_Me_DNA;
  matrix[E,K] phi_Me_TAX;
  matrix[E,K] phi_r_AC;
  matrix[E,K] phi_r_FAC;
  matrix[E,K] phi_r_A;
  matrix[E,K] phi_r_CMF;
}

model {
  matrix[L,E] phi_S1;
  matrix[L,E] phi_S3;
  matrix[K,E] phi_AC;
  matrix[K,E] phi_FAC;
  matrix[K,E] phi_A;
  matrix[K,E] phi_CMF;
  matrix[K,E] phi_T;
  matrix[N_A1,     E] mu_A1_S1;
  matrix[N_A1_AC,  E] mu_A1_AC;
  matrix[N_A1_T,   E] mu_A1_T;
  matrix[N_A2_FAC, E] mu_A2_FAC;
  matrix[N_A3,     E] mu_A3_S3;
  matrix[N_A3_A,   E] mu_A3_A;
  matrix[N_A3_CMF, E] mu_A3_CMF;
  matrix[N_A4,     E] mu_A4_S3;
  matrix[N_A4_AC,  E] mu_A4_AC;
  matrix[N_A4_CMF, E] mu_A4_CMF;

  for (l in 1:L) {
    vector[E] sx_S1 = softmax(phi_S1_raw[,l]);
    vector[E] sx_S3 = softmax(phi_S3_raw[,l]);
    for (e in 1:E) {
      phi_S1[l,e] = sx_S1[e];
      phi_S3[l,e] = sx_S3[e];
    }
    for (e in 1:E_GX) {
      phi_S1_raw[(GX[e,1]+1):GX[e,2], l] ~ normal(phi_S1_raw[GX[e,1]:(GX[e,2]-1), l], Hyp2);
      phi_S3_raw[(GX[e,1]+1):GX[e,2], l] ~ normal(phi_S3_raw[GX[e,1]:(GX[e,2]-1), l], Hyp2);
    }
  }

  for (k in 1:K) {
    vector[E] sx_AC  = softmax(phi_Me_DNA[,k] + phi_r_AC [,k]);
    vector[E] sx_FAC = softmax(phi_Me_DNA[,k] + phi_r_FAC[,k]);
    vector[E] sx_A   = softmax(phi_Me_DNA[,k] + phi_r_A  [,k]);
    vector[E] sx_CMF = softmax(phi_Me_DNA[,k] + phi_r_CMF[,k]);
    vector[E] sx_T   = softmax(phi_Me_TAX[,k]);
    for (e in 1:E) {
      phi_AC [k,e] = sx_AC [e];
      phi_FAC[k,e] = sx_FAC[e];
      phi_A  [k,e] = sx_A  [e];
      phi_CMF[k,e] = sx_CMF[e];
      phi_T  [k,e] = sx_T  [e];
    }
    for (e in 1:E_GX) {
      phi_Me_DNA[(GX[e,1]+1):GX[e,2], k] ~ normal(phi_Me_DNA[GX[e,1]:(GX[e,2]-1), k], Hyp2);
      phi_Me_TAX[(GX[e,1]+1):GX[e,2], k] ~ normal(phi_Me_TAX[GX[e,1]:(GX[e,2]-1), k], Hyp2);
    }
  }

  mu_A1_S1  =                          th_A1_S1                       *            phi_S1;
  mu_A1_AC  = cyc_A1_AC  .* append_col(th_A1_S1,           th_A1_AC)  * append_row(phi_S1, phi_AC);
  mu_A1_T   = cyc_A1_T   .* append_col(th_A1_S1[Nx_A1_Re], th_A1_T)   * append_row(phi_S1, phi_T);
  mu_A2_FAC = cyc_A2_FAC .*                                th_A2_FAC  *                    phi_FAC;
  mu_A3_S3  =                          th_A3_S3                       *            phi_S3;
  mu_A3_A   = cyc_A3_A   .* append_col(th_A3_S3,           th_A3_A)   * append_row(phi_S3, phi_A);
  mu_A3_CMF = cyc_A3_CMF .* append_col(th_A3_S3[Nx_A3_Re], th_A3_CMF) * append_row(phi_S3, phi_CMF);
  mu_A4_S3  =                          th_A4_S3                       *            phi_S3;
  mu_A4_AC  = cyc_A4_AC  .* append_col(th_A4_S3,           th_A4_AC)  * append_row(phi_S3, phi_AC);
  mu_A4_CMF = cyc_A4_CMF .* append_col(th_A4_S3[Nx_A4_Re], th_A4_CMF) * append_row(phi_S3, phi_CMF);

  to_vector(phi_S1_raw[GX0, ]) ~ normal(0, 5);
  to_vector(phi_S3_raw[GX0, ]) ~ normal(0, 5);
  to_vector(phi_Me_DNA[GX0, ]) ~ normal(0, 5);
  to_vector(phi_Me_TAX[GX0, ]) ~ normal(0, 5);
  to_vector(phi_r_AC)  ~ normal(0, Hyp1);
  to_vector(phi_r_FAC) ~ normal(0, Hyp1);
  to_vector(phi_r_A)   ~ normal(0, Hyp1);
  to_vector(phi_r_CMF) ~ normal(0, Hyp1);

  to_array_1d(Y_A1_S1)  ~ poisson(to_vector(mu_A1_S1));
  to_array_1d(Y_A1_AC)  ~ poisson(to_vector(mu_A1_AC));
  to_array_1d(Y_A1_T)   ~ poisson(to_vector(mu_A1_T));
  to_array_1d(Y_A2_FAC) ~ poisson(to_vector(mu_A2_FAC));
  to_array_1d(Y_A3_S3)  ~ poisson(to_vector(mu_A3_S3));
  to_array_1d(Y_A3_A)   ~ poisson(to_vector(mu_A3_A));
  to_array_1d(Y_A3_CMF) ~ poisson(to_vector(mu_A3_CMF));
  to_array_1d(Y_A4_S3)  ~ poisson(to_vector(mu_A4_S3));
  to_array_1d(Y_A4_AC)  ~ poisson(to_vector(mu_A4_AC));
  to_array_1d(Y_A4_CMF) ~ poisson(to_vector(mu_A4_CMF));
}
