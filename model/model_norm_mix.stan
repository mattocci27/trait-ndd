data{
  int<lower=0> N; // number of sample
  int<lower=1> J; // number of sp
  int<lower=1> K; // number of tree-level preditor (i.e, CONS, HETS,...)
  int<lower=1> L; // number of sp-level predictor (i.e., interecept and WP)
  int<lower=1> S; // number of site
  int<lower=1> T; // number of census
  matrix[N, K] x; // tree-level predictor
  row_vector[L] u[J]; // sp-level predictor
  vector[N] y; //
  int<lower=1,upper=J> sp[N]; // integer
  int<lower=1,upper=S> plot[N]; // integer
  int<lower=1,upper=T> census[N]; // integer 
}

parameters{
  matrix[L, K] gamma;
  vector[K] beta[J];
  vector<lower=0>[K] L_sigma;
  cholesky_factor_corr[K] L_Omega;
  real<lower=0> sig;
  vector[S] phi;
  vector[T] tau;
  vector<lower=0>[2] sigma;
}

model {
  vector[N] mu;
  row_vector[K] u_gamma[J];
  // priors
  L_Omega ~ lkj_corr_cholesky(2); // uniform of L_Omega * L_Omega'
  L_sigma ~ cauchy(0, 2.5);
  sig ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 2.5);
  to_vector(gamma) ~ normal(0, 5);
  // transformed parameters
  for (j in 1:J)
    u_gamma[j] = u[j] * gamma;
  for (n in 1:N)
    mu[n] = x[n] * beta[sp[n]] + phi[plot[n]] + tau[census[n]];
  // model
  beta ~ multi_normal_cholesky(u_gamma, diag_pre_multiply(L_sigma, L_Omega));
  phi ~ normal(0, sigma);
/ tau ~ normal(0, sigma[2]);
  y ~ normal(mu, sig);
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
