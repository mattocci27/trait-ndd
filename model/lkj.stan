functions {
  real lkj_cor_lpdf(matrix y, real eta) {
    return lkj_corr_lpdf(y | eta);
  }
  matrix lkj_cor_rng(int K, real eta) {
    return lkj_corr_rng(K, eta);
  }
}
