# seedling-stan

Research project code for seedling survival analyses using Stan.

## Clone this repository

```{bash}
git clone git@github.com:mattocci27/seedling_stan.git
cd seedling_stan
```

## Singularlity

To run MCMC

```
singularity shell ../dockerfiles/singularity/rstan_4.0.2.sif
# in singularity shell
sh sh/model_ind.sh
```

or 

```
singularity exec ../dockerfiles/singularity/rstan_4.0.2.sif sh/model_ind.sh
```


## Docker (not recommended)

It's not recommended for running MCMC on docker but it's convenient for checking html
files on the server.

```{bash}
docker pull mattocci/rstan:4.0.2

docker run --rm -v $(pwd):/home/rstudio/seedling-stan \ 
  -p 8787:8787 \
  -e PASSWORD=< your_password > \
  mattocci/rstan:4.0.2
```

or 

```{bash}
docker pull mattocci/myenv:4.0.2

docker run --rm -v $(pwd):/home/rstudio/seedling-stan \ 
  -p 8787:8787 \
  -e PASSWORD=< your_password > \
  mattocci/myenv:4.0.2

```

Then, go to http://210.72.93.96:8787/:8787/ or localhost:8787/ in your browser.

## Log

20200921
- fixed erros in `run_sran.r` 
- now using 8 cores

20200909

- started running sh/model_ind.sh for the following combinations
  - 2 seasons
  - 5 traits combinations (1. all traits, 2. except for SDMC, 3. except for WD, 4. use PC1 and PC2, 5. use PC1, PC2, PC3)
  - 3 habitats (valley, ridge, slope)
- we may need to run models for all habitat together later 

2020 July-Augst

- New models with this form `tau[k] = 2.5 * tan(tau_unif[k])` worked
- Waiting for new tlp values

20200402
- dry, min abund 30, cc - tlp
  - tlp on height has positive effect
  - 4599.68 seconds on n2s4
  - Rhat for  `lp__` is 1.13 and L_sigma[5] (ahet) is 1.103
  - will use 50 again

- dry, min abund 50, cc - tlp
  - tlp on height has positive effect
  - 5975.35 seconds on n2s4
  - convergence for L_sigma[4-5] (ahet and acon) still looks bad


20200401

- dry, min abund 50, cc, ind effect
    - `sigma[3]`` and `lp__` showed 1.1 < Rhat <1.2
- 8615.59 seconds on n2s4
- wp on cons adulut has negative effect
- I will try model without ind again



