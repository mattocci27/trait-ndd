functions {
  real cauchy_lpdf(real y, real mu, real sigma) {
    return cauchy_lpdf(y | mu, sigma);
  }
  real cauchy_rng(real mu, real sigma) {
    return cauchy_rng(mu, sigma);
  }
}
