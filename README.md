# seedling_stan

## Docker

```{bash}
sudo pull mattocci/rstan
sudo docker run -d -p 8787:8787 -v $(pwd):/home/rstudio -e PASSWORD=<your_password> mattocci/r-debian
```

Then, go to http://xxx.xxx.xx.xx:8787/ or localhost:8787/ in your browser.


## Terminal

```{bash}
sh ./dry_stan.sh
```

## Code

- model.stan
    - changed
    - added function to calculate likelihood at each MCMC step
- dry_stan.r
    - r code to run `model.stan`.
    - this code works manually if you comment out [L12](https://github.com/mattocci27/seedling_stan/blob/a96dff6cb0d85c9826436141f912b9da96e4ca2d/dry_stan.r#L12) and specify [L13-18](https://github.com/mattocci27/seedling_stan/blob/a96dff6cb0d85c9826436141f912b9da96e4ca2d/dry_stan.r#L13-L18)
- dry_stan.sh
    - shell scrip to run `dry_stan`.