#!/usr/bin/env bash
set -e

menu() {
	echo "1) tar_make() on local"
	echo "2) tar_make() on Apptainer"
	echo "3) Enter in the Apptainer container"
	read -rp "Enter number：" menu_num
  case $menu_num in
  1)
    Rscript run.R
    ;;
  2)
 		# singularity exec --env RENV_PATHS_CACHE=/public1/home/${USER}/renv \
 		apptainer exec --env RENV_PATHS_CACHE=/home/${USER}/renv \
		--env RENV_PATHS_PREFIX_AUTO=TRUE \
 		radian.sif Rscript run.R
    ;;
  3)
 		# singularity shell --env RENV_PATHS_CACHE=/public1/home/${USER}/renv \
 		apptainer shell --env RENV_PATHS_CACHE=/home/${USER}/renv \
		--env RENV_PATHS_PREFIX_AUTO=TRUE \
 		radian.sif bash
    ;;
	*)
    echo "Type 1-3"
    ;;
  esac
}

menu "$@"
