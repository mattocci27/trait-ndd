data {
  int<lower=1> N; // number of samples
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditors (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictors (i.e., interecept, SLA,...)
  int<lower=1> S; // number of sites
  int<lower=1> T; // number of censuses
  matrix[N, K] x; // tree-level predictor
  matrix[L, J] u; // sp-level predictor
  array[N] int<lower=0, upper=1> suv; // 1 or 0
  array[N] int<lower=1, upper=J> sp; // integer
  array[N] int<lower=1, upper=S> plot; // integer
  array[N] int<lower=1, upper=T> census; // integer
}

parameters {
  matrix[K, L] gamma;
  matrix[K, J] z;
  vector[S] phi_raw;
  vector[T] xi_raw;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0, upper=pi() / 2>[K] tau_unif;
  vector<lower=0, upper=pi() / 2>[2] sig_unif;
}

transformed parameters {
  matrix[K, J] beta;
  vector<lower=0>[K] tau;
  vector<lower=0>[2] sig;
  vector[S] phi;
  vector[T] xi;
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]);
  for (i in 1:2) sig[i] = 2.5 * tan(sig_unif[i]);
  beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;
  phi = phi_raw * sig[1];
  xi = xi_raw * sig[2];
}

model {
  vector[N] mu;
  to_vector(z) ~ std_normal();
  to_vector(phi_raw) ~ std_normal();
  to_vector(xi_raw) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 2.5);
  for (n in 1:N) {
    mu[n] = x[n, ] * beta[, sp[n]];
  }
  suv ~ bernoulli_logit(mu + phi[plot] + xi[census]);
}

generated quantities {
  vector[N] log_lik;
  corr_matrix[K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(suv[n] | x[n, ] * beta[, sp[n]] +
      phi[plot[n]] + xi[census[n]]);
  }
}
