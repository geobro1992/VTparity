rm(list = ls())

library(tidyverse)
library(lubridate)
library(forcats)
library(spdep)
library(reshape)
library(dplyr)
library(rfishbase)


setwd("I:/VTech/Holly/Data")

#############################################
# AmphiBIO, a global database for amphibians  
#Brunno Freire Oliveira et al.
############################################

dat = read.csv("AmphiBIO_v1.csv")
#dat = dat[complete.cases(dat),]

# calculate mean values for ranges provided
dat$AaM = ((dat$Age_at_maturity_min_y + dat$Age_at_maturity_max_y) / 2)
dat$Litter_size_mean = ((dat$Litter_size_min_n + dat$Litter_size_max_n) / 2)
dat$Offspring_size_mean = ((dat$Offspring_size_min_mm + dat$Offspring_size_max_mm) / 2)
dat$SS = ((dat$Age_at_maturity_max_y - dat$Age_at_maturity_min_y)) / dat$AaM

# calculate juvenile and adult mortality 
# lifetime fecundity and age at maturity as proxies respectively
dat$log.AaM = log(dat$AaM)
dat$juvMort = ifelse((dat$Longevity_max_y - dat$AaM) < 1, 
                     log(dat$Litter_size_mean * dat$Reproductive_output_y),
                     log((dat$Longevity_max_y - dat$AaM) * (dat$Litter_size_mean * dat$Reproductive_output_y)))

dat$adMort = 1/(1 + dat$AaM)

dat$Order = as.factor(dat$Order)


# create factor for reproductive mode
dat[which(dat[,"Dir"] == 0), "Dir"] = NA
dat[which(dat[,"Lar"] == 0), "Lar"] = NA
dat[which(dat[,"Viv"] == 0), "Viv"] = NA
dat[which(dat[,"Ped"] == 0), "Ped"] = NA


repmode = dat[,c("id","Dir", "Lar", "Viv", "Ped")] %>% gather(type, value, -id) %>% na.omit() %>% select(-value)%>% arrange(id)

dat = left_join(dat, repmode)
dat$RepMode = as.factor(dat$type)

summary(dat$RepMode)

# remove species I already have data for (I trust my numbers more than these)
#spp = c("Alytes obstetricans", "Arthroleptis poecilonotus", "Bombina bombina", "Bombina pachypus", "Bombina variegata", 
#        "Anaxyrus americanus", "Rhinella marina", "Anaxyrus quercicus", "Dendrobates auratus", "Nanorana parkeri", "Eleutherodactylus coqui",
#        "Osteopilus septentrionalis", "Pseudacris crucifer", "Hyla gratiosa", "Hyla versicolor", "Gastrophryne carolinensis", "Pelobates fuscus",
#        "Pelodytes punctatus", "Pipa pipa", "Xenopus laevis", "Pyxicephalus adspersus", "Lithobates catesbeianus", "Lithobates pipiens", 
#        "Lithobates clamitans", "Lithobates sylvaticus", "Rana aurora", "Rana latastei", "Scaphiopus couchii", "Scaphiopus holbrookii", 
#        "Ambystoma cingulatum", "Ambystoma gracile", "Ambystoma jeffersonianum", "Ambystoma macrodactylum", "Ambystoma maculatum", "Ambystoma opacum",
#        "Ambystoma talpoideum", "Amphiuma tridactylum", "Andrias japonicus", "Cryptobranchus alleganiensis", "Plethodon glutinosus", "Desmognathus quadramaculatus",
#        "Desmognathus ochrophaeus", "Necturus maculosus", "Proteus anguinus", "Rhyacotriton olympicus", "Notophthalmus viridescens",
#        "Taricha torosa", "Salamandra salamandra", "Salamandrina perspicillata", "Siren intermedia", "Siren lacertina", "Typhlonectes compressicauda")
        
#dat = dat %>%
#  filter(!Species %in% spp)

# extract useful columns
amphibians = dat %>%
  mutate(Class = rep("Amphibia", length.out = length(Order))) %>%
  select(class = "Class", order = "Order", family = "Family", species = "Species", AaM = "AaM", longevity = "Longevity_max_y", "juvMort", "adMort", mass = "Body_mass_g", svl = "Body_size_mm", fecundity = "Litter_size_mean", birth.svl = "Offspring_size_mean", repmode = "RepMode", Terr = "Terr", AdTerr = "AdTerr", IF = "IF")

amphibians$na_count <- apply(amphibians, 1, function(x) sum(is.na(x)))

rm(dat)
rm(repmode)
#rm(spp)

# length mass relationships
lm(log(mass) ~ log(svl) + order, data = amphibians)

amphibians[which(amphibians$order == "Anura"), "mass"] = exp(-7.692 + (2.497*log(amphibians[which(amphibians$order == "Anura"), "svl"]))) 
amphibians[which(amphibians$order == "Caudata"), "mass"] = exp(-2.525 -7.692 + (2.497*log(amphibians[which(amphibians$order == "Caudata"), "svl"])))
amphibians[which(amphibians$order == "Gymnophiona"), "mass"] = exp(-3.711 -7.692 + (2.497*log(amphibians[which(amphibians$order == "Gymnophiona"), "svl"])))


################################
# personally compiled amphibians
dat = read.csv("reptilesBIO.csv")

dat = as.data.frame(dat)
#dat = dat[complete.cases(dat),]

dat = dat %>%
  filter(Class == "Amphibia")

dat$repmode = as.factor(dat$repmode)
dat$log.AaM = log(dat$AaM)
dat$juvMort = ifelse((dat$Longevity - dat$AaM) < 1, 
                     log(dat$Clutches.per.year * dat$Clutch),
                     log((dat$Longevity - dat$AaM) * (dat$Clutches.per.year * dat$Clutch)))

dat$adMort = 1/(1 + dat$AaM)

dat$GS <- paste(dat$Genus, dat$Species)

amps = dat %>%
  select(class = "Class", order = "Order", family = "Family", species = "GS", "AaM", longevity = "Longevity", "juvMort", "adMort", svl = "Max.Size", repmode = "repmode", Terr = "Terr", IF = "IF")

amps$na_count <- apply(amps, 1, function(x) sum(is.na(x)))

rm(dat)

##############################
# Personally compiled reptiles
##############################
dat = read.csv("reptilesBIO.csv")

dat = as.data.frame(dat)
#dat = dat[complete.cases(dat),]
dat = dat %>%
  filter(Class == "Reptilia")

dat$repmode = as.factor(dat$repmode)
dat$log.AaM = log(dat$AaM)

dat$juvMort = ifelse((dat$Longevity - dat$AaM) < 1, 
                     log(dat$Clutches.per.year * dat$Clutch),
                     log((dat$Longevity - dat$AaM) * (dat$Clutches.per.year * dat$Clutch)))

dat$adMort = 1/(1 + dat$AaM)

dat$GS <- paste(dat$Genus, dat$Species)


reptiles = dat %>%
  select(class = "Class", order = "Order", family = "Family", species = "GS", "AaM", longevity = "Longevity", "juvMort", "adMort", svl = "Max.Size", fecundity = "Clutch", birth.svl = "birth.svl", repmode = "repmode", Terr = "Terr", AdTerr = "AdTerr", IF = "IF")

reptiles$na_count <- apply(reptiles, 1, function(x) sum(is.na(x)))

rm(dat)

###################################################
# Meiri, S., 2018. Traits of lizards of the world: 
# Variation around a successful evolutionary design. Global ecology and biogeography, 27(10), pp.1168-1172.
dat = read.csv("lizards2.csv")
#dat = dat[complete.cases(dat),]

dat$maximum.SVL = as.numeric(as.character(dat$maximum.SVL))
dat$female.SVL = as.numeric(as.character(dat$female.SVL))
dat$hatchling.neonate.SVL = as.numeric(as.character(dat$hatchling.neonate.SVL))
dat$Longevity = as.numeric(dat$Longevity)

# if no female specific max size, use max size
dat[is.na(dat[,13]),13] = dat[is.na(dat[,13]),12]

# calculate average clutch size (using range if mean not reported)
dat$Clutch_mean = (dat$smallest.mean.clutch.size + dat$largest.mean.clutch.size) / 2
dat$Clutch_range = (dat$smallest.clutch + dat$largest.clutch) / 2
dat[is.na(dat[,"Clutch_mean"]),"Clutch_mean"] = dat[is.na(dat[,"Clutch_mean"]),"Clutch_range"]

# caculate age at maturity (in years)
dat$AaM = (dat$oldest.age.at.first.breeding..months. + dat$youngest.age.at.first.breeding..months.) / 24
dat$log.AaM = log(dat$AaM)
dat$adMort = 1/(1 + dat$AaM)

dat$SS = (dat$oldest.age.at.first.breeding..months. - dat$youngest.age.at.first.breeding..months.) / dat$AaM

dat$juvMort = ifelse((dat$Longevity - dat$AaM) < 1, 
                     log(dat$Clutch_mean),
                     log((dat$Longevity - dat$AaM) * (dat$Clutch_mean)))

dat$Order = as.factor(dat$Order)
dat$Species = paste(dat$Genus, dat$epithet, sep=" ")

lizards = dat %>%
  mutate(Terr = rep("1", length.out = length(AaM))) %>%
  mutate(AdTerr = rep("1", length.out = length(AaM))) %>%
  mutate(IF = rep("1",length.out = length(AaM))) %>%
  select(class = "X", order = "Order", family = "Family", species = "Species", longevity = "Longevity", "AaM", "adMort", "juvMort",  svl = "female.SVL", birth.svl = "hatchling.neonate.SVL", fecundity = "Clutch_mean", repmode = "reproductive.mode", Terr = "Terr", IF = "IF")

lizards$na_count <- apply(lizards, 1, function(x) sum(is.na(x)))

rm(dat)

####################################
# endotherms
# Nathan P. Myhrvold, Elita Baldridge, Benjamin Chan, Dhileep Sivam, Daniel L. Freeman, and S. K. Morgan Ernest. 2015. 
# An amniote life-history database to perform comparative analyses with birds, mammals, and reptiles. Ecology 96:3109.
dat = read.csv("Mammals.csv")

dat = dat %>%
  filter(class != "Reptilia") %>%
  droplevels()

dat[dat == -999] = NA
dat$AaM = dat$female_maturity_d/365
dat$longevity = rowMeans(dat[,c("longevity_y", "maximum_longevity_y")], na.rm = T)

dat$juvMort = ifelse((dat$longevity - dat$AaM) < 1, 
                     log(dat$litters_or_clutches_per_y * dat$litter_or_clutch_size_n),
                     log((dat$longevity - dat$AaM) * (dat$litters_or_clutches_per_y * dat$litter_or_clutch_size_n)))

dat$adMort = 1/(1 + dat$AaM)

dat$repmode = as.factor(as.numeric(as.factor(dat$class))-1)

dat$SS = (dat$female_maturity_d - dat$male_maturity_d) / dat$female_maturity_d
dat$Species = paste(dat$genus, dat$species, sep=" ")

endos = dat %>% 
  mutate(order = class) %>%
  mutate(IF = rep("1",length.out = length(AaM))) %>%
  select("order", "class", "family", species = "Species", "AaM", "longevity", "juvMort", "adMort", svl = "adult_svl_cm", birth.svl = "birth_or_hatching_svl_cm", mass = "adult_body_mass_g", birth.mass = "birth_or_hatching_weight_g", fecundity = "litter_or_clutch_size_n", repmode = "repmode", Terr = "Terr", AdTerr = "AdTerr", IF = "IF") %>%
  mutate(svl = svl * 10) %>%
  mutate(birth.svl = birth.svl * 10) # convert to mm for consistency

endos[which(endos$family == "Tachyglossidae"),"repmode"] = "0"
endos[which(endos$family == "Ornithorhynchidae"),"repmode"] = "0"

endos$na_count <- apply(endos, 1, function(x) sum(is.na(x)))

rm(dat)

#############################################
# FishTraits, FishBase
############################################

#dat = read.csv("fishdata2021.csv")

#species = unique(dat[,c(2,3,16,17,18)])

#fish = data.frame(species, 
#                  AgeMatMin = as.vector(tapply(dat$AgeMatMin, dat$SpecCode, mean, na.rm=TRUE)),
#                  AgeMatMin2 = as.vector(tapply(dat$AgeMatMin2, dat$SpecCode, mean, na.rm=TRUE)),
#                  FecundityMin = as.vector(tapply(dat$FecundityMin, dat$SpecCode, mean, na.rm=TRUE)),
#                  FecundityMax = as.vector(tapply(dat$FecundityMax, dat$SpecCode, mean, na.rm=TRUE)),
#                  longevity = as.vector(tapply(dat$tmax, dat$SpecCode, mean, na.rm=TRUE)),
#                  birth.svlMin = as.vector(tapply(dat$LengthOffspringMin, dat$SpecCode, mean, na.rm=TRUE)), 
#                  birth.svlMax = as.vector(tapply(dat$LengthOffspringMax, dat$SpecCode, mean, na.rm=TRUE)))
                  
#fish$AaM <- rowMeans(cbind(fish$AgeMatMin, fish$AgeMatMin2), na.rm = T)
#fish$fecundity <- rowMeans(cbind(fish$FecundityMin, fish$FecundityMax), na.rm = T)
#fish$birth.svl <- rowMeans(cbind(fish$birth.svlMin, fish$birth.svlMax), na.rm = T)

#fish$juvMort = ifelse((fish$longevity - fish$AaM) < 1, 
 #                     log(fish$fecundity),
 #                     log((fish$longevity - fish$AaM) * (fish$fecundity)))

#fish$log.AaM = log(fish$AaM)
#fish$adMort = 1/(1 + fish$AaM)

# body sizes from fishbase
#dat2 = rfishbase::species(version="21.04")

#dat2$svl <- rowMeans(cbind(as.numeric(dat2$LengthFemale), dat2$Length), na.rm = T)

#dat2 = dat2 %>% 
#  select(SpecCode, Species, svl, Weight)


# repmode and IF
#dat3 = read.csv("bony_fish_IF.csv")
#dat3 = dat3 %>% 
#  select(SpecCode = "X", "repmode", "fertmode")


# merge fish databases
#fish = merge(fish, dat2, by = "SpecCode",all.x = T) 
#fish = merge(fish, dat3, by = "SpecCode", all.x = T)

# select appropriate columns
#fish = fish %>%
#  mutate(Terr = rep("0", length.out = length(Class))) %>%
#  mutate(AdTerr = rep("0", length.out = length(Class))) %>%
#  select(class = "Class", order = "Order", family = "Family", species = "Species.x", "AaM", longevity = "longevity", "juvMort", "adMort", svl = "svl", fecundity, birth.svl, mass = "Weight", Terr = "Terr", AdTerr = "AdTerr", IF = "fertmode", "repmode")

##fish$na_count <- apply(fish, 1, function(x) sum(is.na(x)))

#rm(dat)
#rm(dat2)
#rm(dat3)
#rm(species)

#write.csv(fish, "bony_fish_JACOB.csv")

#### shark data ######

dat = read.csv("SharkData.csv")
#dat = dat[complete.cases(dat),]


# calculate juvenile and adult mortality 
# lifetime fecundity and age at maturity as proxies respectively
dat$log.AaM = log(dat$age.mat)

dat$juvMort = ifelse((dat$max.age - dat$age.mat) < 1, 
                     log(dat$litter.size * dat$interval),
                     log((dat$max.age - dat$age.mat) * (dat$litter.size * dat$interval)))

dat$adMort = 1/(1 + dat$age.mat)

dat$class = as.factor(dat$ï..superorder)

dat$max.size = dat$max.size * 10
dat$pup.size = dat$pup.size * 10

sharks = dat %>% mutate(Class = rep("Sharks", length.out = length(max.age))) %>%
  mutate(AdTerr = rep("0", length.out = length(max.age))) %>%
  mutate(Terr = rep("0", length.out = length(max.age))) %>%
  mutate(IF = rep("1", length.out = length(max.age))) %>%
  select(class = "Class", order = "Class", family = "family", "species", AaM = "age.mat", longevity = "max.age", "juvMort", "adMort", svl = "max.size", fecundity = "litter.size", birth.svl = "pup.size", repmode = "bear", Terr = "Terr", AdTerr = "AdTerr", IF = "IF")
  
sharks$na_count <- apply(sharks, 1, function(x) sum(is.na(x)))

rm(dat)


################
# bony fish final
dat = read.csv("final_fish.csv")

dat$juvMort = ifelse((dat$longevity - dat$AaM) < 1, 
                     log(dat$fecundity),
                     log((dat$longevity - dat$AaM) * (dat$fecundity)))

dat$adMort = 1/(1 + dat$AaM)

# select appropriate columns
fish = dat %>%
  mutate(Terr = rep("0", length.out = length(class))) %>%
  mutate(AdTerr = rep("0", length.out = length(class))) %>%
  select(class = "class", order = "order", family = "family", species = "species", "AaM", longevity = "longevity", "juvMort", "adMort", svl = "svl", fecundity, birth.svl, mass = "mass", Terr = "Terr", AdTerr = "AdTerr", IF = "fertmode", "repmode")

fish$na_count <- apply(fish, 1, function(x) sum(is.na(x)))

write.csv(x = fish, file = "vtparity_bonyfish.csv")

rm(dat)
################
# merge datasets
x = merge(endos, amps, all = TRUE)
x = merge(x, amphibians, all = TRUE)
x = merge(x, reptiles, all = TRUE)
x = merge(x, lizards, all = TRUE)
x = merge(x, fish, all = TRUE)
x = merge(x, sharks, all = TRUE)



x = x %>% 
  mutate(repmode = recode(repmode, 
                          `Dir`="0",
                          `Viv`="1",
                          `Lar`="0",
                          `Ped`="0",
                          `Mixed`="1",
                          `Oviparous`="0",
                          `unclear`="0",
                          `Viviparous`="1",
                          `live-bearing`="1",
                          `egg-laying`="0")) %>%
  droplevels()

x = x[!duplicated(x),]
x = x[!duplicated(x$species), ]

table(x$repmode, x$class)
table(x$adMort, x$class)

write.csv(x = x, file = "vtparity.csv")

dat = read.csv("vtparity.csv")
dat = dat[which(dat$class == "Actinopterygii"),]
summary(dat$repmode)

comp.x = x[which(x$juvMort != "NA" & x$adMort != "NA" & x$svl != "NA"),]

table(comp.x$class, comp.x$repmode)



#######################
# mass-length equations

# bony fish
dat = rfishbase::length_weight(version="21.04")

dat = dat %>%
  select(Species, a, b)

# lizards
dat2 = read.csv("lizard_mass_length_eqs.csv")

dat2 = dat2 %>%
  select(Species = "Binomial", a = "intercept", b = "slope")

# sharks
dat3 = read.csv("SharkData.csv")

dat3 = dat3 %>%
  select(Species = "species", a = "length.weight.a", b = "length.weight.b")


x2 = rbind(dat, dat2, dat3)
x2 = x2[complete.cases(x2),]
write.csv(x2, "VTparity_length_mass.csv")
