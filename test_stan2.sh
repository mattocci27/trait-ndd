#!/bin/bash
SPAB=50
DRY=dry
echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
export SPAB DRY
  # n_iter n_warm n_thin n_chain adapt_delta n_ab season
R --vanilla --slave --args 2 1 1 1 0.9 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}.log
