# R-package-EcotoneFinder: Characterising and Locating Ecotones and Communities 
This repository contains the necssary files to build and use the *EcotoneFinder* package on R.

This package regroups analytical methods to locate and characterise ecotones, ecosystems and environmental patchiness along ecological gradients. Methods are implemented for both isolated sampling and space/time series. It includes Detrended Correspondence Analysis (Hill & Gauch (1980) <doi:10.1007/BF00048870>), fuzzy clustering (De Cáceres et al. (2010) <doi:10.1080/01621459.1963.10500845>), biodiversity indices (Jost (2006) <doi:10.1111/j.2006.0030-1299.14714.x>), and network analyses (Epskamp et al. (2012) <doi:10.18637/jss.v048.i04>) -- as well as tools to explore the number of clusters in the data. Functions to produce synthetic ecological datasets are also provided.

# Instalation:
*The updated EcotoneFinder package will be submitted to CRAN once vignette compatibility issues accross platforms are solved*

In the meantime:

From source file:

`# install.packages("path-to-file/EcotoneFinder_0.2.4.tar.gz", repos = NULL, type = "source")`

From GitHub:

`# library("devtools"); install_github("Ecotoni/EcotoneFinder",dependencies=TRUE)`

# Documentation:
Vignette for the package is provided: EcotoneFinder-vignette.Rmd & EcotoneFinder-vignette.pdf.
