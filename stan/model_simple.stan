data{
  int<lower=0> N; // number of sample
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> J; // number of sp
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  int<lower=1,upper=J> sp[N]; // integer: sp group for trees
  matrix[N, K] x; // tree-level predictor
  //row_vector[L] u[J]; // sp-level predictor
  matrix[J, L] u;
  int<lower=0,upper=1> suv[N]; // 1 or 0
}

parameters{
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
  matrix[L, K] gamma;
}

transformed parameters{
  matrix[J, K] beta;
  vector<lower=0>[K] tau;
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]); // implies tau ~ cauchy(0, 2.5)
  beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)';
}

model {
//  vector[N] p;
  // Hyper-priors
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2); // uniform of L_Omega * L_Omega'
  // transformed parameters
  // Prior
  to_vector(gamma) ~ normal(0, 5);
  // Likelihood
 // suv ~ bernoulli_logit(p);
  suv ~ bernoulli_logit(rows_dot_product(beta[sp] , x));
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

