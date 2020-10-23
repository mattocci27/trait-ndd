#!/bin/bash
singularity exec $HOME/dockerfiles/singularity/rstan_4.0.2.sif R --args $n < ./sh/check_loo0.r
for n in {1..48}
  do 
  singularity exec $HOME/dockerfiles/singularity/rstan_4.0.2.sif R --args $n < ./sh/check_loo.r
done
