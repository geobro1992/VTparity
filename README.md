# VTparity
Supporting data and code for the manuscript: 
Brooks, G. C., Uyeda J. C., H. Conrad, N. Bone, C. Mull, and H. K. Kindsvater. 2024. Fundamental constraints on the evolution of vertebrate life histories.
[link](https://doi.org/10.1101/2024.01.23.576873)

## Contents
**datasets** - compiled vertebrate life-history data (body size, age at maturity, reproductive output, reproductive mode)\
**trees** - phylogenetic trees for each vertebrate clade\
**code** - R scripts to create the supertree and perform the discrete-state Markovian evolutionary model of vertebrate life-history transitions\
**output** - model output and figures

## Instructions\
- download repository and unzip in local drive\
- download Rtools from cran.r-project.org/bin/windows/Rtools\
- open _code.Rproj_ file in the 'code' folder\
- install required packages: lazyeval, ape, rotl, castor, geiger, phytools, phylolm, rmarkdown, dplyr, plyr, Rphylopars, MCMCpack, Ternary, treeplyr*

*NOTE: package treeplyr has been removed from CRAN and so needs to be installed from source using Rtools. We have provided the _tar.gz_ file in the code folder.\
Once Rtools has been installed and the .Rproj file has been opened, the following code should manually install the package from your local drive:\ 
install.packages("treeplyr_0.1.10.tar.gz", repos = NULL, type="source")

**Run scripts in this order:**\
run _tree_merge.Rmd_ to create the full phylogeny and match life-history data to tree tips\
run _figs.R_ to create the raw data figures and basic regression\
run _analysis.Rmd_ to fit the ternary model and plot model predictions
