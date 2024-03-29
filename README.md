[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Seedling Survival Model

Research project code for seedling survival analyses using Stan.

## Data

The required data files for this project are available for download at ScienceDB (https://doi.org/10.57760/sciencedb.02276).

## Reproduce the Results

The workflow and code (in R and Stan) for this analysis are managed with the `targets` R package, which can be found here: [ropensci/targets](https://github.com/ropensci/targets).
Analyses can be conducted either locally or within an container (Docker and Apptainer).

For instructions on using the containerized environment, refer to the [Reproducible R Project Template](https://github.com/mattocci27/rdv-template).

### Apptainer Container Setup

Create an Apptainer container with the following command, which includes R, Quarto, and cmdstan:

```sh
sudo apptainer build radian.sif radian.def
```

### Running the Code
Before starting the analysis, run `./make_renviron.sh` to create a `.Renviron` file.
This helps to manage package caching efficiently.

To initiate the analysis, use the `./run.sh` script:

```sh
# For the first run, to install R packages
# Rscript -e "renv::restore()"

> ./run.sh
1) tar_make() on local (or inside Docker)
2) tar_make() on Apptainer
3) Enter the Apptainer container
4) Enter the Singularity container on HPC
Enter number:
```

Use `tar_make()` for single-threaded execution and `tar_make_clustermq()` for multi-threaded, parallel processing.

Please update the following:
**Computational Intensity Warning:**
Due to the computationally intensive nature of some analyses, you may need to execute only a subset of the targets pipeline, especially when system memory (RAM) is a limiting factor.
Some slurm scripts for High-Performance Computing (HPC) are also available under `slurm`.

**Docker Usage for Less Intensive Analysis:**
For parts of the analysis that are less computationally demanding, Docker can be utilized, as explained in the [Reproducible R Project Template](https://github.com/mattocci27/rdv-template).

### Requirements

- Apptainer (or Singularity)
- Docker
- A Linux-based OS
