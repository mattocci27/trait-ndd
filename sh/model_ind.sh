#!/bin/bash
SPAB=50
MODEL=model_ind_gq
for DRY in dry wet
do
  for TRAIT in C13 LT WP PCA
  do
    for NEIGHBOR in full seedling adult
    do
      echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
      echo "Run ${MODEL}.stan"
      echo "Use ${NEIGHBOR} data"
      echo "Use ${TRAIT} for sp-level"
      export SPAB DRY MODEL NEIGHBOR
        # n_iter n_warm n_thin n_chain adapt_delta n_ab season
      nohup R --vanilla --slave --args ${MODEL} 2000 1000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}.log &
      sleep 1 # pause to be kind to the scheduler
    done
    wait
  done
  wait
done
