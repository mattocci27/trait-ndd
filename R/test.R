targets::tar_load(dry_full)
targets::tar_load(dry_wd)
targets::tar_load(dry_pca)

str(dry_full)
str(dry_wd)
str(dry_pca)


model_code <-
  '
  functions {
    vector test() {
      return rows_dot_product(beta[sp] , x);
   }
  }
'
