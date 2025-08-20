#!/bin/bash

Rscript -e  "install.packages(\"remotes\", repos=\"http://lib.stat.cmu.edu/R/CRAN/\")"
Rscript -e  "remotes::install_version(\"BiocManager\", \"1.30.22\", repos=\"http://lib.stat.cmu.edu/R/CRAN/\")"
Rscript -e  "install.packages(\"renv\", version = \"1.0.11\", repos=\"http://lib.stat.cmu.edu/R/CRAN/\")"
Rscript -e  "renv::init()"
cp renv_mac.lock renv.lock
Rscript -e  "renv::restore()"