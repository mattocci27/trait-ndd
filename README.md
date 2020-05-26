# seedling_stan

Clone this repository.

```{bash}
git clone git@github.com:mattocci27/seedling_stan.git
cd seedling_stan
```


## Docker

```{bash}
sudo docker pull mattocci/rstan
sudo docker run -d -p 8787:8787 -v $(pwd):/home/rstudio -e PASSWORD=<your_password> mattocci/rstan
```

Then, go to http://xxx.xxx.xx.xx:8787/ or localhost:8787/ in your browser.


## Local terminal or terminal inside the docker.

To see if the model can be compiled:

```{bash}
sh ./test_stan.sh
```

To run the model:

```{bash}
sh ./run_stan.sh
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
