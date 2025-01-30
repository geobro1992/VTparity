source("./phylotern/phyloternFunctions_reduced.R")
library(parallel)
library(ape)
library(treeplyr)
library(castor)
library(geiger)
library(phytools)
library(phylolm)
library(Ternary)
library(Rphylopars)
library(optimParallel)
library(nloptr)

fits <- list()
terntd <- readRDS("../output/terntd80822.RDS")
completetd <- filter(terntd, !is.na(totT))

counts_complete <- round(count_ambiguous(completetd$dat$ddat, 16),0)
counts_ambiguous <- round(count_ambiguous(terntd$dat$ddat, 16),0)

pi <- plotTernaryIndex(4)
#terncols <- Tricolore(data.frame(pi$cells), p1="a", p2="b", p3="c")
#terncols <- setNames(terncols$rgb, pi$labels)
#terncols2 <- terncols[as.character(1:16)]

#########Global Fits###############
cl <- makeCluster(20, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pps <- list(c(-0.6, 0.7,  -5, 0.01))

lowerbound <- c(-30,-30,-30,  1e-8)
upperbound <- c( 10, 10, 10, 10)
upperbound.draw <- c( 10,  10, 10, 0.1)


for(i in 2:50){
  pps[[i]] <- runif(length(lowerbound), lowerbound, upperbound.draw)
}


tree0 <- completetd$phy
dat0 <- completetd[['ddat']]


lnLfx0 <- make.lnLTernaryGrad(tree0, dat0, bins=4)
lnLfx0(pps[[1]])

tree1 <- terntd$phy
dat1 <- terntd[['ddat']]
lnLfx1 <- make.lnLTernaryGrad(tree1, dat1, bins =4)
#lnLfx1(c(0,0,0,0,0.1))

fit0s <- list()
fit1s <- list()

for(i in 1:length(pps)){
  
  fit0s[[i]] <- nloptr(pps[[i]], lnLfx0, opts=opts, lb = lowerbound, ub=upperbound)
  fit1s[[i]] <- nloptr(fit0s[[i]]$solution, lnLfx1, opts=opts, lb = lowerbound, ub=upperbound)
  
  #fit0s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
  #                            lower=lowerbound, upper=upperbound, 
  #                            parallel=list(loginfo=TRUE))
  
  #fit1s[[i]] <- optimParallel(par=fit0s[[i]]$par, fn=lnLfx1, 
  #                            lower=lowerbound,
  #                            upper=upperbound, parallel=list(loginfo=TRUE))
  print(i)
}

o <-order(sapply(fit0s, function(x) x$objective))
o2 <-order(sapply(fit1s, function(x) x$objective))
sapply(fit0s, function(x) x$objective)[o]
sapply(fit1s, function(x) x$objective)[o2]


#f0Qs <- do.call(cbind, lapply(fit0s[o], function(x){qq <- r2Q(x$solution[1:3], x$solution[4], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit1s[o2], function(x){qq <- r2Q(x$solution[1:3], x$solution[4], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

bestmodel0 <- fit0s[o][[1]]
bestmodel1 <- fit1s[o2][[1]]

pdf('../output/GlobalTernary_121224.pdf')
q0 <- plotTernaryGradient(tree0, bestmodel0$solution[1:3], bestmodel0$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
.simcounts <- round(10*log(table(factor(q0$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q1 <- plotTernaryGradient(tree1, bestmodel1$solution[1:3], bestmodel1$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
.simcounts <- round(10*log(table(factor(q1$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
dev.off()

fits$global_complete <- bestmodel0
fits$global_all <- bestmodel1



#########Aquatic & Terrestrial###############
cl <- makeCluster(20, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster
bins <- 4
pps <- list(c(c(-0.6, 0.7,  -5, 0.01), c(-0.6, 0.7,  -5, 0.01), 0.05, 0.05))

lowerbound <- c(-30,-30,-30,  1e-8,-30,-30,-30,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10, 10,  10, 10, 10)
upperbound.draw <- c( 10,  10, 10, 0.1, 10,  10, 10, 0.1, 0.1, 0.1)

for(i in 2:50){
  pps[[i]] <- runif(length(lowerbound), lowerbound, upperbound.draw)
}


tree0 <- completetd$phy
dat0 <- completetd$dat$ddat_AqTerr #layerTrait(completetd[['ddat']], completetd[['Terr']], bins)
lnLfx0 <- make.lnL2xTernaryGrad(tree0, dat0, bins=4, model="ARD")
lnLfx0(pps[[1]])

tree1 <- terntd$phy
dat1 <- terntd$dat$ddat_AqTerr #layerTrait(terntd[['ddat']], terntd[['Terr']], bins)
lnLfx1 <- make.lnL2xTernaryGrad(tree1, dat1, bins =4, model="ARD")
lnLfx1(pps[[1]])

fit2s <- list()
fit3s <- list()

opts <- list("algorithm"="NLOPT_LN_SBPLX",
             "xtol_rel"=1.0e-6, 
             "maxeval"=1000)

for(i in 2:20){
  #if(i>length(fit2s)){
    fit2s[[i]] <- nloptr(pps[[i]], lnLfx0, opts=opts, lb = lowerbound, ub=upperbound)#optimParallel(pps[[i]], lnLfx0, 
  #}
                              #lower=lowerbound, upper=upperbound, 
                              #parallel=list(loginfo=TRUE))
  #fits2s[[length(pps)+i]] <- 
  fit3s[[i]] <- nloptr(pps2[[i]], lnLfx1, opts=opts, lb = lowerbound, ub=upperbound)#reran with ending pps because hit bounds
  #fit3s[[i]] <- nloptr(fit2s[[i]]$solution, lnLfx1, opts=opts, lb = lowerbound, ub=upperbound)
                                #optimParallel(par=fit2s[[i]]$par, fn=lnLfx1, 
  #                            lower=lowerbound,
  #                            upper=upperbound, parallel=list(loginfo=TRUE))
  print(i)
}
o <-order(sapply(fit2s, function(x) x$objective))
o2 <-order(sapply(fit3s, function(x) x$objective))
sapply(fit2s, function(x) x$objective)[o]
sapply(fit3s, function(x) x$objective)[o2]

#f0Qs <- do.call(cbind, lapply(fit2s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit3s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)
# 
bestmodel0 <- fit2s[o][[1]]
bestmodel1 <- fit3s[o2][[1]]

pdf('../output/AqTerrTernary_112724.pdf')
q2a <- plotTernaryGradient(tree0, bestmodel0$solution[1:3], bestmodel0$solution[4], bins=4, line.weight=50, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q2b <- plotTernaryGradient(tree0, bestmodel0$solution[5:7], bestmodel0$solution[8], bins=4, line.weight=50, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])


q3a <- plotTernaryGradient(tree1, bestmodel1$solution[1:3], bestmodel1$solution[4], bins=4, line.weight=10, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q3b <- plotTernaryGradient(tree1, bestmodel1$solution[5:7], bestmodel1$solution[8], bins=4, line.weight=10, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
dev.off()

fits$aqterr_complete <- bestmodel0
fits$aqterr_all <- bestmodel1


#########Aquatic & Terrestrial (ADULT ONLY)###############
cl <- makeCluster(20, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster
bins <- 4
pps <- list(c(c(-0.6, 0.7,  -5, 0.01), c(-0.6, 0.7,  -5, 0.01), 0.05, 0.05))


lowerbound <- c(-30,-30,-30,  1e-8,-30,-30,-30,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10, 10,  10, 10, 10)
upperbound.draw <- c( 10,  10, 10, 0.1, 10,  10, 10, 0.1, 0.1, 0.1)

for(i in 2:50){
  pps[[i]] <- runif(length(lowerbound), lowerbound, upperbound.draw)
}


tree0 <- completetd$phy
dat0 <- layerTrait(completetd[['ddat']], completetd[['AdTerr']], bins)
lnLfx0 <- make.lnL2xTernaryGrad(tree0, dat0, bins=4, model="ARD")
lnLfx0(pps[[1]])

tree1 <- terntd$phy
dat1 <- layerTrait(terntd[['ddat']], terntd[['AdTerr']], bins)
lnLfx1 <- make.lnL2xTernaryGrad(tree1, dat1, bins =4, model="ARD")
lnLfx1(pps[[1]])

fit2xs <- list()
fit3xs <- list()

opts <- list("algorithm"="NLOPT_LN_SBPLX",
             "xtol_rel"=1.0e-6, 
             "maxeval"=1000)

for(i in 1:length(pps)){
  #fit2xs[[i]] <- nloptr(pps[[i]], lnLfx0, opts=opts, lb = lowerbound, ub=upperbound)#optimParallel(pps[[i]], lnLfx0, 
  #lower=lowerbound, upper=upperbound, 
  #parallel=list(loginfo=TRUE))
  #fits2s[[length(pps)+i]] <- 
  
  fit3xs[[i]] <- nloptr(fit2s[[i]]$solution, lnLfx1, opts=opts, lb = lowerbound, ub=upperbound)
  
  print(i)
}

o <-order(sapply(fit2xs, function(x) x$objective))
o2 <-order(sapply(fit3xs[1:3], function(x) x$objective))
sapply(fit2xs, function(x) x$objective)[o]
sapply(fit3xs, function(x) x$objective)[o2]

#f0Qs <- do.call(cbind, lapply(fit2s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit3s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

bestmodel0 <- fit2xs[o][[1]]
bestmodel1 <- fit3xs[o2][[1]]
sapply(fit2xs, function(x) x$objective)[o]
sapply(fit3xs, function(x) x$objective)[o2]

par(mfrow=c(1,2))
pdf('../output/Aq_AdTerrTernary_12062024.pdf')
q2a <- plotTernaryGradient(tree0, bestmodel0$solution[1:3], bestmodel0$solution[4], bins=4, line.weight=50, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q2b <- plotTernaryGradient(tree0, bestmodel0$solution[5:7], bestmodel0$solution[8], bins=4, line.weight=50, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])


q3a <- plotTernaryGradient(tree1, bestmodel1$solution[1:3], bestmodel1$solution[4], bins=4, line.weight=10, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q3b <- plotTernaryGradient(tree1, bestmodel1$solution[5:7], bestmodel1$solution[8], bins=4, line.weight=10, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
dev.off()

fits$aqadterr_complete <- bestmodel0
fits$aqadterr_all <- bestmodel1


#########Parity Mode###############
bins <- 4
pps <- list(c(c(-0.6, 0.7, -5, 0.01), c(-0.6, 0.7,  -5, 0.01), 0.05, 0.05))



lowerbound <- c(-30,-30,-30,  1e-8,-30,-30,-30,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10, 10,  10, 10, 10)
upperbound.draw <- c( 10,  10, 10, 0.1, 10,  10, 10, 0.1, 0.1, 0.1)

for(i in 2:50){
  pps[[i]] <- runif(length(lowerbound), lowerbound, upperbound.draw)
}


tree0 <- completetd$phy
dat0 <- completetd$dat$ddat_repmode #layerTrait(completetd[['ddat']], completetd[['repmode']], bins)
lnLfx0 <- make.lnL2xTernaryGrad(tree0, dat0, bins=4, model="ARD")
lnLfx0(pps[[1]])

tree1 <- terntd$phy
dat1 <- terntd$dat$ddat_repmode #layerTrait(terntd[['ddat']], terntd[['repmode']], bins)
lnLfx1 <- make.lnL2xTernaryGrad(tree1, dat1, bins =4, model="ARD")
lnLfx1(pps[[1]])

fit4s <- list()
fit5s <- list()

for(i in 1:length(pps)){
  fit4s[[i]] <- nloptr(pps[[i]], lnLfx0, opts=opts, lb = lowerbound, ub=upperbound)#optimParallel(pps[[i]], lnLfx0, 
  fit5s[[i]] <- nloptr(pps2[[i]], lnLfx1, opts=opts, lb = lowerbound, ub=upperbound)#reran with ending pps because hit bounds
  print(i)
  #fit4s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
  #                            lower=lowerbound, upper=upperbound, 
  #                            parallel=list(loginfo=TRUE))
  
  #fit5s[[i]] <- optimParallel(par=fit4s[[i]]$par, fn=lnLfx1, 
   #                           lower=lowerbound,
  #                            upper=upperbound, parallel=list(loginfo=TRUE))
}

o <-order(sapply(fit4s, function(x) x$objective))
o2 <-order(sapply(fit5s, function(x) x$objective))

#f0Qs <- do.call(cbind, lapply(fit2s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit3s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

bestmodel0 <- fit4s[o][[1]]
bestmodel1 <- fit5s[o2][[1]]

par(mfrow=c(1,2))
pdf('../output/RepModeTernary_120624.pdf')
q4a <- plotTernaryGradient(tree0, bestmodel0$solution[1:3], bestmodel0$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=1000)
.simcounts <- round(10*log(table(factor(q4a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q4b <- plotTernaryGradient(tree0, bestmodel0$solution[5:7], bestmodel0$solution[8], bins=4, line.weight=100, simulate.data=TRUE, simreps=1000)
.simcounts <- round(10*log(table(factor(q4b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])


q5a <- plotTernaryGradient(tree1, bestmodel1$solution[1:3], bestmodel1$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=1000)
.simcounts <- round(10*log(table(factor(q5a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q5b <- plotTernaryGradient(tree1, bestmodel1$solution[5:7], bestmodel1$solution[8], bins=4, line.weight=100, simulate.data=TRUE, simreps=1000)
.simcounts <- round(10*log(table(factor(q5b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
dev.off()

fits$repmode_complete <- bestmodel0
fits$repmode_all <- bestmodel1



##### Global fits with trait, for model comparison ########
#########Global Fits###############
cl <- makeCluster(20, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster


lowerbound <- c(-30,-30,-30,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10)


tree0 <- completetd$phy
dat_aqterr0 <- completetd$dat$ddat_AqTerr #layerTrait(completetd[['ddat']], completetd[['Terr']], bins)
dat_aqadterr0 <- layerTrait(completetd[['ddat']], completetd[['AdTerr']], bins)
dat_repmode0 <- completetd$dat$ddat_repmode
lnLfx0.aqterr <- make.lnL2xIdenticalTernaryGrad(tree0, dat_aqterr0, bins=4, model="ARD")
lnLfx0.aqadterr <- make.lnL2xIdenticalTernaryGrad(tree0, dat_aqadterr0, bins=4, model="ARD")
lnLfx0.repmode <- make.lnL2xIdenticalTernaryGrad(tree0, dat_repmode0, bins=4, model="ARD")

tree1 <- terntd$phy
dat_aqterr1 <- terntd$dat$ddat_AqTerr #layerTrait(terntd[['ddat']], terntd[['Terr']], bins)
dat_aqadterr1 <- layerTrait(terntd[['ddat']], terntd[['AdTerr']], bins)
dat_repmode1 <- terntd$dat$ddat_repmode

lnLfx1.aqterr <- make.lnL2xIdenticalTernaryGrad(tree1, dat_aqterr1, bins=4, model="ARD")
lnLfx1.aqadterr <- make.lnL2xIdenticalTernaryGrad(tree1, dat_aqadterr1, bins=4, model="ARD")
lnLfx1.repmode <- make.lnL2xIdenticalTernaryGrad(tree1, dat_repmode1, bins=4, model="ARD")


fit0s.aqterr_id <- nloptr(c(fits$global_complete$solution, 0.0003229057, 0.0001208332), lnLfx0.aqterr, opts=opts, lb = lowerbound, ub=upperbound)
fit1s.aqterr_id <- nloptr(c(fits$global_all$solution, 0.0001977236, 9.359797e-05), lnLfx1.aqterr, opts=opts, lb = lowerbound, ub=upperbound)
fit0s.aqadterr_id <- nloptr(c(fits$global_complete$solution, 0.0002, 0.0003), lnLfx0.aqadterr, opts=opts, lb = lowerbound, ub=upperbound)
fit1s.aqadterr_id <- nloptr(c(fits$global_all$solution, 0.0001, 0.0001), lnLfx1.aqadterr, opts=opts, lb = lowerbound, ub=upperbound)
fit0s.repmode_id <- nloptr(c(fits$global_complete$solution, 0.0004, 0.0002), lnLfx0.repmode, opts=opts, lb = lowerbound, ub=upperbound)
fit1s.repmode_id <- nloptr(c(fits$global_all$solution, 0.0002, 0.003), lnLfx1.repmode, opts=opts, lb = lowerbound, ub=upperbound)

#fit0s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
#                            lower=lowerbound, upper=upperbound, 
#                            parallel=list(loginfo=TRUE))

#fit1s[[i]] <- optimParallel(par=fit0s[[i]]$par, fn=lnLfx1, 
#                            lower=lowerbound,
#                            upper=upperbound, parallel=list(loginfo=TRUE))
#  print(i)
#}

#o <-order(sapply(fit0s, function(x) x$objective))
#o2 <-order(sapply(fit1s, function(x) x$objective))
#sapply(fit0s, function(x) x$objective)[o]
#sapply(fit1s, function(x) x$objective)[o2]


#f0Qs <- do.call(cbind, lapply(fit0s[o], function(x){qq <- r2Q(x$solution[1:3], x$solution[4], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit1s[o2], function(x){qq <- r2Q(x$solution[1:3], x$solution[4], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

#bestmodel0 <- fit0s.id[[1]]
#bestmodel1 <- fit1s.id[[1]]

#pdf('../output/GlobalTernary_121224.pdf')
#q0 <- plotTernaryGradient(tree0, bestmodel0$solution[1:3], bestmodel0$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
#.simcounts <- round(10*log(table(factor(q0$simdat$tip_states, 1:16))),0)
#.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
#TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])##

#q1 <- plotTernaryGradient(tree1, bestmodel1$solution[1:3], bestmodel1$solution[4], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
#.simcounts <- round(10*log(table(factor(q1$simdat$tip_states, 1:16))),0)
#.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
#TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
#dev.off()

fits$aqterr_complete_independent <- fit0s.aqterr_id
fits$aqterr_all_independent <- fit1s.aqterr_id
fits$aqadterr_complete_independent <- fit0s.aqadterr_id
fits$aqadterr_all_independent <- fit1s.aqadterr_id
fits$repmode_complete_independent <- fit0s.repmode_id
fits$repmode_all_indpendent <- fit1s.repmode_id

saveRDS(fits, "../output/fits_12122024.rds")
