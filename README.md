# Seedling survival model

Research project code for seedling survival analyses using Stan.

## Reproduce the results

Codes (R and STAN) and workflow are managed with the R package `targets` (https://github.com/ropensci/targets).
Analysis can be run locally or inside an Apptainer container.

### Creating an Apptainer container

The following command creates an Apptainer container that contains R, quarto and, cmdstan.

```
sudo apptainer build radian.sif radian.def
```

### Running code

Run `./make_renviron.sh` to create `.Renviron` first to use renv's cache in your home directory `$HOME/renv/`.

Run `./run.sh` to run the analysis.

```bash
# To install R packages for the first run
# Rscript -e "renv::restore()"

> ./run.sh

1) tar_make() on local
2) tar_make_clustermq() on local
3) tar_make() on Apptainer
4) tar_make_clustermq() on Apptainer
5) Enter in the Apptainer container
6) Enter in the Singularity container on HPC
Enter number：
```

`tar_make()` runs a pipeline with a single thread, and
`tar_make_clustermq()` runs a pipeline in parallel.

Requirements to run on your local computer:

- cmdstan 2.29.2
- quarto
- latexdiff
- R (4.2.1)
	- renv (`renv::restore()` will install all the R packages)

Requirements to run inside the Apptainer container:

- Apptainer (or Singularity)
- Linux
