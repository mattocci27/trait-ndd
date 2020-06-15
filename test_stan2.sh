#!/bin/bash
SPAB=50
DRY=dry
MODEL=model_ind_gq
NEIGHBOR=adult
echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
echo "Run ${MODEL}.stan"

export SPAB DRY MODEL
  # n_iter n_warm n_thin n_chain adapt_delta n_ab season
R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} < run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}.log
#R --vanilla --slave --args ${MODEL} 3000 2000 1 4 0.95 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}.log
