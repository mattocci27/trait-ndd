#!/bin/bash

for SPAB in 1 10 20
  do
  echo "min. species abundance = ${SPAB}"
  export SPAB
  # n_iter n_warm n_thin n_chain adapt_delta n_ab
  #R --vanilla --slave --args 2 1 1 1 0.9 ${SPAB} < dry_stan.r > ./log/dry_stan_${SPAB}.log
  nohup R --vanilla --slave --args 4000 3000 1 4 0.9 ${SPAB} < dry_stan.r > ./log/dry_stan_${SPAB}.log &
done
