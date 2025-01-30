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
download
unzip
open _code.Rproj_ file in code folder
install required packages: ape, treeplyr, rotl, castor, geiger, phytools, phylolm, rmarkdown (NOTE: package treeplyr may need to be installed from source)

run _tree_merge.Rmd_ to create the full phylogeny and match life-history data to tree tips
