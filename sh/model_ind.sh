#!/bin/bash
SPAB=50
MODEL=model_ind

for DRY in dry wet
do
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  export SPAB DRY MODEL
    # n_iter n_warm n_thin n_chain adapt_delta n_ab season
  nohup R --vanilla --slave --args ${MODEL} 2000 1000 1 4 0.95 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}.log &
sleep 1 # pause to be kind to the scheduler
done
