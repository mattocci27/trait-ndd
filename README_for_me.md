# seedling_stan

Clone this repository.


build docker

```
docker build --add-host="api.github.com:140.82.113.6" -t mattocci/cmdstanr .
```

build sing

```
sudo SINGULARITY_NOHTTPS=true singularity build gpu_stan.sif gpu_stan.def

singularity exec --nv gpu_stan.sif bash sh/test.sh

```


```{bash}

kr pair

sudo fallocate -l 32G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

sudo echo "/swapfile swap swap defaults 0 0" >> /etc/fstab

git clone git@github.com:mattocci27/seedling_stan.git
cd seedling_stan

```


## Docker

```{bash}

sudo docker load < mattocci_rstan_3.6.3.tar

sudo docker pull mattocci/rstan

docker run -d -p 8787:8787 -v $(pwd):/home/rstudio -e PASSWORD=mogemoge --name rstudio mattocci/rstan:3.6.3


docker run -p 8787:8787 -v $(pwd):/home/rstudio/seedling-stan -e PASSWORD=mogemoge mattocci/myenv:3.6.3

docker run -it --rm -v $(pwd):/home/rstudio -u rstudio mattocci/myenv:3.6.3 /bin/bash


F85hPRItkcsaQ7lR6AHK

docker run --rm -v $(pwd):/home/rstudio/seedling-stan -p 8787:8787 -e PASSWORD=F85hPRItkcsaQ7lR6AHK mattocci/myenv:4.0.2

docker run --rm -v $(pwd):/home/rstudio/seedling-stan -p 59000:8787 -e PASSWORD=F85hPRItkcsaQ7lR6AHK mattocci/myenv:4.0.2

#sudo docker run -it -v $(pwd):/home/rstudio/seedling_stan -u rstudio mattocci/rstan R


docker run -p 8787:8787 -v $(pwd):/home/rstudio/seedling-stan -e PASSWORD=mogemoge mattocci/myenv:3.6.3


4cpu 16GB
$0.1664 ~ 1.18 RMB / h
```

Then, go to http://xxx.xxx.xx.xx:8787/ or localhost:8787/ in your browser.

## Docker X11

```

SOCK=/tmp/.X11-unix
XAUTH=$HOME/.Xauthority
xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge -
chmod 777 $XAUTH

docker-compose run stan /bin/bash

```

## Singularlity

```
singularity exec ../r-containers/rstan_4.0.2.sif sh/model_ind.sh

singularity shell ../r-containers/myenv_4.0.2.sif
```


## Local terminal or terminal inside the docker.

To see if the model can be compiled:

```{bash}
sh ./test_stan.sh
```

To run the model:

```{bash}
sh ./run_stan.sh

sh ./dry50.sh

```

## Code

- model_ind.stan
    - added function to calculate likelihood at each MCMC step
    - added random intercept for each tag
- `run_stan.r`
    - r code to run `model_ind.stan`.
    - this code works manually if you comment out [L12](https://github.com/mattocci27/seedling_stan/blob/2d065e240222943a0abc6b68df3839e6fa3eaef4/run_stan.r#L12) and specify [L13-18](https://github.com/mattocci27/seedling_stan/blob/2d065e240222943a0abc6b68df3839e6fa3eaef4/run_stan.r#L13-L18).
- `run_stan.sh` and `run_test.sh`
    - shell scrip to run `run_stan.r`.
    - logs will be stored in `/log/` and rda files will be stored in `/data`


```{r}


fs <- data.frame(summary(fit)$summary)

fs %>%
  mutate(para = rownames(.)) %>%
  DT::datatable(.)

mcmc_trace(fit, pars = c("gamma[2,1]", "gamma[2,2]","gamma[2,3]", "gamma[2,4]","gamma[2,5]", "gamma[2,6]"))

```

## log

20200909

- sh/model_ind.sh (for each habita)

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



