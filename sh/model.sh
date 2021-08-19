#!/bin/bash
SPAB=50
MODEL=model_inter
NEIGHBOR=full
for TRAIT in Full WD PCA
#for HAB in valley
do
  #for TRAIT in PC2
  for DRY in dry rainy
  do
    #for HAB in valley ridge slope all
    for HAB in Full valley ridge slope
    do
      echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
      echo "Run ${MODEL}.stan"
      echo "Use ${NEIGHBOR} data"
      echo "Use ${TRAIT} for sp-level"
      echo "Habitat: ${HAB}"
      export SPAB DRY MODEL NEIGHBOR HAB TRAIT
        # n_iter n_warm n_thin n_chain adapt_delta n_ab season
      nohup R --vanilla --slave --args ${MODEL} 4000 2000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
     # nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
      sleep 1 # pause to be kind to the scheduler
      #wait
    done
    wait
  done
  wait
done