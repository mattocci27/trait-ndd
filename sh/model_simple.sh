#!/bin/bash
SPAB=50
MODEL=model_inter
NEIGHBOR=full
TRAIT=WD
 #for TRAIT in PC2
HAB=Full
for DRY in dry rainy
do
  echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
  echo "Run ${MODEL}.stan"
  echo "Use ${NEIGHBOR} data"
  echo "Use ${TRAIT} for sp-level"
  echo "Habitat: ${HAB}"
  export SPAB DRY MODEL NEIGHBOR HAB TRAIT
    # n_iter n_warm n_thin n_chain adapt_delta n_ab season
  nohup R --vanilla --slave --args ${MODEL} 4000 2000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan_simple.R > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}_simple.log &
 # nohup R --vanilla --slave --args ${MODEL} 20 10 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan_simple.R > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}_simple.log &
  sleep 1 # pause to be kind to the scheduler
  #wait
done
