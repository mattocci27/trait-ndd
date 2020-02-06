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