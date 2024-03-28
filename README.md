[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Seedling survival model

Research project code for seedling survival analyses using Stan.

## Data

Data files required for this project can be downloaded from <https://doi.org/10.57760/sciencedb.02276>.

## Reproduce the results

Codes (R and STAN) and workflow are managed with the R package `targets` (https://github.com/ropensci/targets).
Analysis can be run locally or inside an Apptainer container.

### Apptainer Container Setup

To create an Apptainer container (which includes R, Quarto, and cmdstan), use the following command:

```
sudo apptainer build radian.sif radian.def
```

### Running code

Prior to running the analysis, execute `./make_renviron.sh` to generate a `.Renviron file`.
This is essential for utilizing renv's cache located in your home directory `$HOME/renv/`.

Next, initiate the analysis using the `./run.sh` command:

```
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

The function `tar_make()` executes a pipeline using a single thread, while `tar_make_clustermq()` runs the pipeline in parallel.

Requirements for Local Execution:

- cmdstan 2.33.1
- Quarto
- latexdiff
- R (4.3.1)
	- renv (`renv::restore()` will install all the R packages)

Requirements for Execution within the Apptainer Container:

- Apptainer (or Singularity)
- Linux OS
