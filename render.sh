#!/bin/bash

nohup Rscript -e 'library(rmarkdown); system.time(render("./rmd/model.rmd"))' &

nohup Rscript -e 'library(rmarkdown); system.time(render("./rmd/model_simple.rmd"))' &
