#rm(list=ls(all=TRUE))
library(shiftPlot)
library(ape)
library(treeplyr)
library(rotl)
library(castor)
library(geiger)
library(phytools)
library(phylolm)
library(rphylopic)
library(Ternary)
library(Rphylopars)
library(optimParallel)
library(tricolore)
library(colorspace)
#dat <- read.csv("../datasets/vtparity.csv")
#td0 <- readRDS("../output/chordatedata_matched_new.rds")
#newmass <- readRDS("../output/NB_massEst_NoOutliers.RDS")
#mergemass <- left_join(td0$dat, newmass)
#mergemass <- mutate(mergemass, lnSize=ifelse(is.na(logMass)==TRUE, logEstMass/3, logMass/3))
#mergemass$lnSize[mergemass$class=="Actinopterygii" & is.na(mergemass$lnSize)] <- (3.0904 * mergemass$svl[mergemass$class=="Actinopterygii" & is.na(mergemass$lnSize)]^0.007)/3
#mergemass$species <- td0$phy$tip.label
#td <- make.treedata(td0$phy, mergemass)
getwd()
td <-readRDS("../output/chordatedata_matched_estSizeComplete.rds")


#Data patches
td$dat$repmode[td$dat$family=="Plethodontidae" & !is.na(td$dat$repmode) & td$dat$repmode==1] <- 0
td$dat$svl[td$phy$tip.label=="Plethodon_cinereus"] <- 112
#td$dat$Terr[td$phy$tip.label=="Desmognathus_ocoee"] <- NA
#td$dat$Terr[td$phy$tip.label=="Desmognathus_wrighti"] <- 1
#td$dat$Terr[td$phy$tip.label=="Desmognathus_organi"] <- 1
#td$dat$Terr[td$phy$tip.label=="Desmognathus_ochrophaeus"] <- 1
#td$dat$Terr[td$phy$tip.label=="Desmognathus_orestes"] <- 1
#td$dat$Terr[grep("Plethodon_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Oedipina_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Batrachoseps_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Bolitoglossa_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Aneides_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Thorius_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Dendrotriton_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Chiropterotriton_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Cryptotriton_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Nototriton_", td$phy$tip.label)] <- 1
#td$dat$Terr[grep("Pseudoeurycea_", td$phy$tip.label)] <- 1
td$dat$AdTerr[grep("Pipa_", td$phy$tip.label)] <- 0
td$dat$Terr[grep("Pipa_", td$phy$tip.label)] <- 0
td$dat$species <- td$phy$tip.label

source("./phylotern/phyloternFunctions.R")

completetd <- filter(td, !is.na(adMort), !is.na(lnSize), !is.na(juvMort))
completeRanges <- c(diff(quantile(td$dat$lnSize, c(0, 1), na.rm=TRUE)), 
                    diff(quantile(td$dat$adMort, c(0, 0.9), na.rm=TRUE)),
                    diff(quantile(td$dat$juvMort, c(0, 1), na.rm=TRUE)))


terntd <- td %>% filter(is.na(lnSize) + is.na(adMort) + is.na(juvMort) + is.na(repmode) < 3) %>% mutate(scSize = (lnSize-min(lnSize,na.rm=TRUE))/completeRanges[1], scadMort= ((adMort-min(adMort,na.rm=TRUE)))/completeRanges[2], scjMort=(juvMort-min(juvMort,na.rm=TRUE))/completeRanges[3]) %>% mutate(totT = scSize + scadMort + scjMort)
quantile(terntd$dat$scadMort, na.rm=TRUE, c(0, 0.01, 0.025, 0.5, 0.975, 0.99, 1))
quantile(terntd$dat$scSize, na.rm=TRUE, c(0, 0.01, 0.025, 0.5, 0.975, 0.99, 1))
quantile(terntd$dat$scjMort, na.rm=TRUE, c(0, 0.01, 0.025, 0.5, 0.975, 0.99, 1))


totDat <- data.frame("species"=terntd$phy$tip.label, totT = as.vector(terntd[['totT']]))
totFit <- Rphylopars::phylopars(totDat, force.ultrametric(terntd$phy, "extend"), model="OU")
terntd$dat$itotT <- totFit$anc_recon[1:nrow(totDat)]

terntd <- terntd %>% mutate(trAdMort = scadMort/itotT, trJMort = scjMort/itotT, trSize = scSize/itotT, totNew=trAdMort + trJMort + trSize)

terntd$dat <- recodeTraits(terntd$dat, c("trAdMort", "trJMort", "trSize"), bins=4)
terntd <- filter(terntd, ddat!="", ddat!="1&2&3&4&5&6&7&8&9&10&11&12&13&14&15&16")
terntd$dat$ddat_AqTerr <- layerTrait(terntd[['ddat']], terntd[['Terr']], bins=4)
terntd$dat$ddat_repmode <- layerTrait(terntd[['ddat']], terntd[['repmode']], bins=4)
completetd <- filter(terntd, !is.na(totT))
dplyr::select(terntd, lnSize, adMort, juvMort, scSize, scadMort, scjMort, itotT, totT, trSize, trAdMort, trJMort, repmode, ddat,ddat_AqTerr, ddat_repmode, totNew)
table(terntd$dat$ddat)
counts_complete <- round(count_ambiguous(completetd$dat$ddat, 16),0)
counts_ambiguous <- round(count_ambiguous(terntd$dat$ddat, 16),0)

saveRDS(terntd, "../output/terntd80822.RDS")

