#!/bin/bash
SPAB=30
DRY=dry
for MODEL in model model_ind
  do 
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  export SPAB DRY MODEL
  #model  n_iter n_warm n_thin n_chain adapt_delta n_ab season
  nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.99 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}.log &
  sleep 1
done
