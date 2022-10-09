# Seedling survival model

Research project code for seedling survival analyses using Stan.

## Reproduce the results

Codes (R and STAN) and workflow are managed with the R package `targets` (https://github.com/ropensci/targets).

### Running code on local

To run analysis:

```bash
# To install R packages for the first run
# Rscript -e "renv::restore()"
Rscript run.R
```

To generate the method descripton:

```bash
make
```

Requirements:

- cmdstan 2.29.2
- pandoc
- pandoc-crossref
- latexdiff
- R (4.2.1)
	- renv (`renv::restore()` will install all the R packages)

### Running code in Apptainer (Linux)

First, change `RENV_PATHS_CACHE` in `radian.def` and `tinytex.def` to your path (i.e.,
`
RENV_PATHS_CACHE=<your_path>"
`
).

To build Apptainer containers:

```bash
sudo apptainer build radian.sif radian.def
```

To run analysis:

```bash
# To install R packages for the first run
# apptainer exec apptainer exec --env RENV_PATHS_CACHE=/home/mattocci/renv \
#		--env RENV_PATHS_PREFIX_AUTO=TRUE \
#	 radian.sif Rscript -e "renv::restore()"
apptainer exec --env RENV_PATHS_CACHE=/home/mattocci/renv \
	--env RENV_PATHS_PREFIX_AUTO=TRUE \
	radian.sif Rscript run.R
```

To generate the method descripton:

```bash
apptainer exec --env RENV_PATHS_CACHE=/home/mattocci/renv \
	--env RENV_PATHS_PREFIX_AUTO=TRUE \
	radian.sif make
```

Requirements:

- Apptainer (or Singularity)
