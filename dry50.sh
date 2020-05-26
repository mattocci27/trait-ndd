#!/bin/bash
SPAB=50
DRY=dry
for MODEL in model_ind model
  do 
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  export SPAB DRY MODEL
  #model  n_iter n_warm n_thin n_chain adapt_delta n_ab season
  R --vanilla --slave --args ${MODEL} 3000 2000 1 4 0.95 ${SPAB} ${DRY} < run_stan.r > ./log2/${DRY}_stan_${SPAB}_${MODEL}.log
  sleep 1
done
