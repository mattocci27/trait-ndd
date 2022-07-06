functions {
  real lkj_cor_lpdf(matrix y, real eta) {
    return lkj_corr_lpdf(y | eta);
  }

  matrix lkj_cor_rng(int K, real eta) {
    return lkj_corr_rng(K, eta);
  }
  real cauchy__lpdf(real y, real mu, real sigma) {
    return cauchy_lpdf(y | mu, sigma);
  }
  real cauchy__rng(real mu, real sigma) {
    return cauchy_rng(mu, sigma);
  }
}
