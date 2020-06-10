data{
  int<lower=0> N; // number of sample
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  int<lower=1> M; // number of seedling individuals (tag)
  int<lower=1> S; // number of site
  int<lower=1> T; // number of census
  matrix[N, K] x; // tree-level predictor
  matrix[J, L] u; // sp-level predictor
  int<lower=0,upper=1> suv[N]; // 1 or 0
  int<lower=1> sp[N]; // integer
  int<lower=1> plot[N]; // integer
  int<lower=1> census[N]; // integer
  int<lower=1> tag[N]; // integer
}

parameters{
  matrix[K, J] z;
  vector[S] phi_raw;
  vector[M] psi_raw;
  vector[T] xi;
  matrix[L, K] gamma;
  vector<lower=0>[2] sigma;
  vector<lower=0>[K] L_sigma;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
  vector<lower=0,upper=pi()/2>[2] sig_unif;
}

transformed parameters{
  matrix[J, K] beta;
  vector<lower=0>[K] tau;
  vector<lower=0>[2] sig;
  vector[S] phi;
  vector[M] psi;
  for (k in 1:K) tau[k] = 2.5 * tan(tau_unif[k]); // implies tau ~ cauchy(0, 2.5)
  for (i in 1:2) sig[i] = 2.5 * tan(sig_unif[i]); // implies sig ~ cauchy(0, 2.5)
  beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)';
  phi = phi_raw * sig[1];
  psi = psi_raw * sig[2];
}

model {
  // Hyper-priors
  to_vector(z) ~ std_normal();
  to_vector(phi_raw) ~ std_normal();
  to_vector(psi_raw) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2); // uniform of L_Omega * L_Omega'
  // Priors
  to_vector(gamma) ~ normal(0, 5);
  xi ~ normal(0, 5);
  // Likelihood
  suv ~ bernoulli_logit(rows_dot_product(beta[sp] , x) + phi[plot] + xi[census] + psi[tag]);
}

generated quantities {
  corr_matrix[K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}


