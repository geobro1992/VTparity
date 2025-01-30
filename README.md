# VTparity
Supporting data and code for the manuscript: 
Brooks, G. C., Uyeda J. C., H. Conrad, N. Bone, C. Mull, and H. K. Kindsvater. 2024. Fundamental constraints on the evolution of vertebrate life histories.
[link](https://doi.org/10.1101/2024.01.23.576873)

## Contents
**datasets** - compiled vertebrate life-history data (body size, age at maturity, reproductive output, reproductive mode)\
**trees** - phylogenetic trees for each vertebrate clade\
**code** - R scripts to create the supertree and perform the discrete-state Markovian evolutionary model of vertebrate life-history transitions\
**output** - model output and figures

## Instructions
- download repository and unzip in local drive\
- download Rtools from cran.r-project.org/bin/windows/Rtools\
- open _code.Rproj_ file in the 'code' folder\
- install required packages: ape, rotl, castor, geiger, phytools, phylolm, rmarkdown, dplyr, plyr, Rphylopars, MCMCpack, Ternary, lazyeval, ComplexUpset, kableExtra, tricolore, UpSetR, viridis, treeplyr*

*NOTE: package treeplyr has been removed from CRAN and so needs to be installed from source using Rtools. We have provided the _tar.gz_ file in the code folder.\
Once Rtools has been installed and the .Rproj file has been opened, the following code should manually install the package from your local drive:\ 
install.packages("treeplyr_0.1.10.tar.gz", repos = NULL, type="source")

**If running from scratch (without calling saved objects already created in the 'output' folder), you will need to run scripts in this order:**
- _tree_merge.Rmd_\
- _figs.R_\
- _phylo_analysis_summary.Rmd_
