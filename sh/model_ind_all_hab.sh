#!/bin/bash
SPAB=50
MODEL=model_ind
NEIGHBOR=full
for HAB in all valley ridge slope
#for HAB in valley
do
  for TRAIT in SLA_WD_TLP
  #for TRAIT in SLA_WD_TLP Full SDMC WD PC2 PC3
  #for TRAIT in PC2
  do
    #for HAB in valley ridge slope all
    for DRY in dry rainy
    do
      echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
      echo "Run ${MODEL}.stan"
      echo "Use ${NEIGHBOR} data"
      echo "Use ${TRAIT} for sp-level"
      echo "Habitat: ${HAB}"
      export SPAB DRY MODEL NEIGHBOR HAB TRAIT
        # n_iter n_warm n_thin n_chain adapt_delta n_ab season
      nohup R --vanilla --slave --args ${MODEL} 2000 1000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
     # nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
      sleep 1 # pause to be kind to the scheduler
      #wait
    done
    wait
  done
  wait
done


for HAB in all
do
  for TRAIT in Full SDMC WD PC2 PC3
  do
    #for HAB in valley ridge slope all
    for DRY in dry rainy
    do
      echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
      echo "Run ${MODEL}.stan"
      echo "Use ${NEIGHBOR} data"
      echo "Use ${TRAIT} for sp-level"
      echo "Habitat: ${HAB}"
      export SPAB DRY MODEL NEIGHBOR HAB TRAIT
        # n_iter n_warm n_thin n_chain adapt_delta n_ab season
      nohup R --vanilla --slave --args ${MODEL} 2000 1000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
     # nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_stan.r > ./log/${DRY}_stan_${SPAB}_${MODEL}_${NEIGHBOR}_${TRAIT}_${HAB}.log &
      sleep 1 # pause to be kind to the scheduler
      #wait
    done
    wait
  done
  wait
done