data{
  int<lower=0> N; // number of sample
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  int<lower=1> M; // number of seedling individuals (tag)
  int<lower=1> S; // number of site
  int<lower=1> T; // number of census
  matrix[N, K] x; // tree-level predictor
  matrix[L, J] u; // sp-level predictor
  array[N] int<lower=0,upper=1> suv; // 1 or 0
  array[N] int<lower=1,upper=J> sp; // integer
  array[N] int<lower=1,upper=S> plot; // integer
  array[N] int<lower=1,upper=T> census; // integer
  array[N] int<lower=1> tag; // integer
}

parameters{
  matrix[K, L] gamma;
  matrix[K, J] z;
  vector[S] phi_raw;
  vector[T] xi_raw;
  vector[M] psi_raw;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
  vector<lower=0,upper=pi()/2>[3] sig_unif;
}

transformed parameters{
  matrix[K, J] beta;
  vector<lower=0>[K] tau;
  vector<lower=0>[3] sig;
  vector[S] phi;
  vector[T] xi;
  vector[M] psi;
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]);
  for (i in 1:3) sig[i] = 2.5 * tan(sig_unif[i]);
  beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;
  phi = phi_raw * sig[1];
  xi = xi_raw * sig[2];
  psi = psi_raw * sig[3];
}

model {
  vector[N] mu;
  to_vector(z) ~ std_normal();
  to_vector(phi_raw) ~ std_normal();
  to_vector(xi_raw) ~ std_normal();
  to_vector(psi_raw) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 5);
  for (n in 1:N) {
    mu[n] = x[n, ] * beta[, sp[n]];
  }
  suv ~ bernoulli_logit(mu + phi[plot] + xi[census] + psi[tag]);
}


