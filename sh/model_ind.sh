#!/bin/bash
SPAB=50
MODEL=model_ind
DRY=dry
for NEIGHBOR in full seedling adult
do
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  echo "Use ${NEIGHBOR} data"
  export SPAB DRY MODEL NEIGHBOR
    # n_iter n_warm n_thin n_chain adapt_delta n_ab season
  nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}.log &
  sleep 1 # pause to be kind to the scheduler
done
wait

DRY=wet
for NEIGHBOR in full seedling adult
do
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  echo "Use ${NEIGHBOR} data"
  export SPAB DRY MODEL NEIGHBOR
    # n_iter n_warm n_thin n_chain adapt_delta n_ab season
  nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}.log &
  sleep 1 # pause to be kind to the scheduler
done
wait
