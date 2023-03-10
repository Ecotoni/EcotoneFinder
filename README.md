# R-package-EcotoneFinder: Characterising and Locating Ecotones and Communities 
This repository contains the necssary files to build and use the *EcotoneFinder* package on R.

## Description:
This package regroups analytical methods to locate and characterise ecotones, ecosystems and environmental patchiness along ecological gradients.\
Methods are implemented for both isolated sampling and space/time series. It includes Detrended Correspondence Analysis [[1]](#1), fuzzy clustering [[2]](#2), biodiversity indices [[3]](#3), and network analyses [[4]](#4) – as well as tools to explore the number of clusters in the data.\
Functions to produce synthetic ecological datasets are also provided.

## Aim:
The initial motivation to put together the set of analyses implemented in this package was an attempt at addressing one of the shortcomings presented by Ries et al., 2017 [[5]](#5) in their review on ecotone research, namely that there was no consistent framework for the exploration of edge characteristics at the ecosystem/community level. Alongside with proposing such a framework, we tried to make the outputs of the *EcotoneFinder package* as compatible as possible with the single-species framework presented earlier, in 2004, by the same authors [[6]](#6).\
The basic idea of this framework is to delineate ecotonal regions through the two wdges of adjacent communities. Simple metrics can be extracted from the shapes of these edges (e.g. width of the ecotone, magnitude of change, sharpness and regularity of change), and provide standardised elements of comparison for studies on ecotones accross ecological systems.\
Within a singular study, the characteristics of these edges may easily be compared with other aspets of the ecosystem, e.g. the evolution of environmetal variables (to investigate the relative importance of biotic *versus* abiotic factors in shaping the ecotone), or the response of different types of organisms (typically at other trophic levels) to the ecotonal region.\
Multivariate definition of community edges also avoids most arbitrary decisions regarding ecotones positions and extents. The framework presented here is thus particularly relevant for systems where the delineations between community types are not obvious, either due to the size of the organisms (making visual assessments impossible), the apparent homogeneity of the environment, the existence of closely related groups (making the classification of near-boundary objects unsatisfactory), or the existence of important mismatches between pre-existing classifications and the reality (driving the necessity for re-evaluation). These issues often arise in aquatic environments and microbial ecosystems. Both these cases were of primary interest when designing this package.\

## Instalation:
*The updated EcotoneFinder package will be submitted to CRAN once vignette compatibility issues accross platforms are solved*

In the meantime:
From source file:\
`# install.packages("path-to-file/EcotoneFinder_0.2.4.tar.gz", repos = NULL, type = "source")`

From GitHub:\
`# library("devtools"); install_github("Ecotoni/EcotoneFinder",dependencies=TRUE)`

## Documentation:
Vignette for the package is provided: EcotoneFinder-vignette.Rmd & EcotoneFinder-vignette.pdf.

### References:
<a id="1">[1]</a> 
Hill MO, Gauch JJr HG (1980) Detrended correspondence analysis: An improved ordination technique. Vegetatio 42:47–58 *doi:10.1007/BF00048870*\
<a id="2">[2]</a> 
De Cáceres M, Font X, Oliva F (2010) The management of vegetation classifications with fuzzy clustering. Journal of Vegetation Science 21:1138–1151 *doi:10.1080/01621459.1963.10500845*\
<a id="3">[3]</a> 
Lou Jost (2006) Entropy and diversity. Oikos 113:363–375 *doi:10.1111/j.2006.0030-1299.14714.x*\
<a id="4">[4]</a> 
Epskamp S, Cramer AOJ, Waldorp LJ, Schmittmann VD, Borsboom D (2012) qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software 48:1–18 *doi:10.18637/jss.v048.i04*\
<a id="5">[5]</a> 
Ries L, Murphy SM, Wimp GM, Fletcher RJ (2017) Closing Persistent Gaps in Knowledge About Edge Ecology. Curr Landscape Ecol Rep:1–12 *doi:10.1007/s40823-017-0022-4*\
<a id="6">[6]</a> 
Ries L, Sisk TD (2004) A predictive model of edge effects. Ecology 85:2917–2926\
