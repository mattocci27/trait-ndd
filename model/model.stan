data{
  int<lower=1> N; // number of sample
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  int<lower=1> S; // number of site
  int<lower=1> T; // number of census
  matrix[N, K] x; // tree-level predictor
  row_vector[L] u[J]; // sp-level predictor
  int<lower=0> suv[N]; // 1 or 0
  int<lower=1> sp[N]; // integer
  int<lower=1> plot[N]; // integer
  int<lower=1> census[N]; // integer
}

parameters{
  matrix[L, K] gamma;
  vector[K] beta[J];
  vector[S] phi;
  vector[T] tau;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_tau;
  vector<lower=0>[K] L_sigma;
  cholesky_factor_corr[K] L_Omega;
}

model {
  vector[N] p;
  row_vector[K] u_gamma[J];
  // priors
  sigma_phi ~ cauchy(0, 5);
  sigma_tau ~ cauchy(0, 5);
  L_Omega ~ lkj_corr_cholesky(2); // uniform of L_Omega * L_Omega'
  L_sigma ~ cauchy(0, 5);
  to_vector(gamma) ~ normal(0, 5);
  // transformed parameters
  for (j in 1:J)
    u_gamma[j] = u[j] * gamma;
  for (n in 1:N)
    p[n] = x[n] * beta[sp[n]] + phi[plot[n]] + tau[census[n]];
  // model
  beta ~ multi_normal_cholesky(u_gamma, diag_pre_multiply(L_sigma, L_Omega));
  phi ~ normal(0, sigma_phi);
  tau ~ normal(0, sigma_tau);
  suv ~ bernoulli_logit(p);
}

