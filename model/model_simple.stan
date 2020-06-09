data{
  int<lower=0> N; // number of sample
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  matrix[N, K] x; // tree-level predictor
  row_vector[L] u[J]; // sp-level predictor
  int<lower=0,upper=1> suv[N]; // 1 or 0
  int<lower=1,upper=J> sp[N]; // integer
}

parameters{
  matrix[L, K] gamma;
  vector[K] beta[J];
  vector<lower=0>[K] L_sigma;
  cholesky_factor_corr[K] L_Omega;
}

model {
  vector[N] p;
  row_vector[K] u_gamma[J];
  // Hyper-priors
  L_Omega ~ lkj_corr_cholesky(2); // uniform of L_Omega * L_Omega'
  L_sigma ~ cauchy(0, 2.5);
  // transformed parameters
  for (j in 1:J)
    u_gamma[j] = u[j] * gamma;
  for (n in 1:N)
    p[n] = x[n] * beta[sp[n]];
  // Prior
  to_vector(gamma) ~ normal(0, 2.5);
  beta ~ multi_normal_cholesky(u_gamma, diag_pre_multiply(L_sigma, L_Omega));
  // Likelihood
  suv ~ bernoulli_logit(p);
}

generated quantities {
  corr_matrix[K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

//generated quantities {
//  vector[N] log_lik;
//  vector[N] p;
//  row_vector[K] u_gamma[J];
//  // transformed parameters
//  for (j in 1:J)
//    u_gamma[j] = u[j] * gamma;
//  for (n in 1:N) {
//    p[n] = x[n] * beta[sp[n]] + phi[plot[n]] + tau[census[n]];
//    log_lik[n] = bernoulli_logit_lpmf(suv[n] | p[n]);
//  }
// }

