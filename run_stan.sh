#!/bin/bash
for SPAB in 1 10 20
  do 
  for DRY in dry wet
  do
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  export SPAB DRY
  # n_iter n_warm n_thin n_chain adapt_delta n_ab season
#  R --vanilla --slave --args 2 1 1 1 0.9 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}.log
  #nohup R --vanilla --slave --args 4000 3000 1 4 0.9 ${SPAB} < run_stan.r > ./log/run_stan_${SPAB}.log &
  R --vanilla --slave --args 4000 3000 1 4 0.9 ${SPAB} ${DRY} < run_stan.r > ./log/${DRY}_stan_${SPAB}.log
  sleep 1
  done
done

