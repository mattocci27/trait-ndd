#!/bin/bash
SPAB=50
MODEL=model_inter
NEIGHBOR=full
TRAIT=WD
HAB=all
DRY=dry

echo "Model for ${DRY} season. Min. species abundance = ${SPAB}"
echo "Run ${MODEL}.stan"
echo "Use ${NEIGHBOR} data"
echo "Use ${TRAIT} for sp-level"
echo "Habitat: ${HAB}"
export SPAB DRY MODEL NEIGHBOR HAB TRAIT
  # n_iter n_warm n_thin n_chain adapt_delta n_ab season
nohup R --vanilla --slave --args ${MODEL} 2 1 1 1 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_cmdstan.r > ./log/test.log &
#nohup R --vanilla --slave --args ${MODEL} 4000 2000 1 4 0.95 ${SPAB} ${DRY} ${NEIGHBOR} ${TRAIT} ${HAB}< run_cmdstan.r > ./log/test.log &
wait