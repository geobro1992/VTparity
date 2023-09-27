source("./phylotern/phyloternFunctions.R")
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
cl <- makeCluster(30, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

pps <- list(c(-0.6, 0.7, -0.2, -5, 0.01))

lowerbound <- c(-10,-10,-10,-10,1e-8)
upperbound <-c(10,10,10,10,100)
upperbound.draw <- c(10,10,10,10, 0.01)

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
  fit0s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
                              lower=lowerbound, upper=upperbound, 
                              parallel=list(loginfo=TRUE))
  
  fit1s[[i]] <- optimParallel(par=fit0s[[i]]$par, fn=lnLfx1, 
                              lower=lowerbound,
                              upper=upperbound, parallel=list(loginfo=TRUE))
  print(i)
}

o <-order(sapply(fit0s, function(x) x$value))
o2 <-order(sapply(fit1s, function(x) x$value))


f0Qs <- do.call(cbind, lapply(fit0s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
cor(f0Qs) 
f1Qs <- do.call(cbind, lapply(fit1s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
cor(f1Qs)

bestmodel0 <- fit0s[o][[1]]
bestmodel1 <- fit1s[o2][[1]]

pdf('../output/GlobalTernary.pdf')
q0 <- plotTernaryGradient(tree0, bestmodel0$par[1:4], bestmodel0$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
.simcounts <- round(10*log(table(factor(q0$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q1 <- plotTernaryGradient(tree1, bestmodel1$par[1:4], bestmodel1$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=25)
.simcounts <- round(10*log(table(factor(q1$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
dev.off()

fits$global_complete <- bestmodel0
fits$global_all <- bestmodel1

#########Aquatic & Terrestrial###############
cl <- makeCluster(30, type="FORK")     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster
bins <- 4
pps <- list(c(c(-0.6, 0.7, -0.2, -5, 0.01), c(-0.6, 0.7, -0.2, -5, 0.01), 0.05, 0.05))

lowerbound <- c(-10,-10,-10,-10,  1e-8,-10,-10,-10,-10,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
upperbound.draw <- c( 10, 10, 10, 10, 0.1, 10, 10, 10, 10, 0.1, 0.1, 0.1)

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

for(i in 1:length(pps)){
  fit2s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
                              lower=lowerbound, upper=upperbound, 
                              parallel=list(loginfo=TRUE))
  
  fit3s[[i]] <- optimParallel(par=fit2s[[i]]$par, fn=lnLfx1, 
                              lower=lowerbound,
                              upper=upperbound, parallel=list(loginfo=TRUE))
  print(i)
}

o <-order(sapply(fit2s, function(x) x$value))
o2 <-order(sapply(fit3s, function(x) x$value))

#f0Qs <- do.call(cbind, lapply(fit2s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit3s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

bestmodel0 <- fit2s[o][[1]]
bestmodel1 <- fit3s[o2][[1]]

par(mfrow=c(1,2))
#pdf('../output/AqTerrTernary.pdf')
q2a <- plotTernaryGradient(tree0, bestmodel0$par[1:4], bestmodel0$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q2b <- plotTernaryGradient(tree0, bestmodel0$par[6:9], bestmodel0$par[10], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q2b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])


q3a <- plotTernaryGradient(tree1, bestmodel1$par[1:4], bestmodel1$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q3b <- plotTernaryGradient(tree1, bestmodel1$par[6:9], bestmodel1$par[10], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q3b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
#dev.off()

fits$aqterr_complete <- bestmodel0
fits$aqterr_all <- bestmodel1

#########Parity Mode###############
bins <- 4
pps <- list(c(c(-0.6, 0.7, -0.2, -5, 0.01), c(-0.6, 0.7, -0.2, -5, 0.01), 0.05, 0.05))

lowerbound <- c(-10,-10,-10,-10,  1e-8,-10,-10,-10,-10,  1e-8, 1e-8, 1e-8)
upperbound <- c( 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
upperbound.draw <- c( 10, 10, 10, 10, 0.1, 10, 10, 10, 10, 0.1, 0.1, 0.1)

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
  fit4s[[i]] <- optimParallel(pps[[i]], lnLfx0, 
                              lower=lowerbound, upper=upperbound, 
                              parallel=list(loginfo=TRUE))
  
  fit5s[[i]] <- optimParallel(par=fit4s[[i]]$par, fn=lnLfx1, 
                              lower=lowerbound,
                              upper=upperbound, parallel=list(loginfo=TRUE))
  print(i)
}

o <-order(sapply(fit4s, function(x) x$value))
o2 <-order(sapply(fit5s, function(x) x$value))

#f0Qs <- do.call(cbind, lapply(fit2s[o], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f0Qs) 
#f1Qs <- do.call(cbind, lapply(fit3s[o2], function(x){qq <- r2Q(x$par[1:4], x$par[5], bins=4); diag(qq) <- 0; v <- as.vector(qq); v<- v[v!=0]}))
#cor(f1Qs)

bestmodel0 <- fit4s[o][[1]]
bestmodel1 <- fit5s[o2][[1]]
bestmodel0$value
bestmodel1$value

par(mfrow=c(1,2))
#pdf('../output/RepModeTernary.pdf')
q4a <- plotTernaryGradient(tree0, bestmodel0$par[1:4], bestmodel0$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q4a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q4b <- plotTernaryGradient(tree0, bestmodel0$par[6:9], bestmodel0$par[10], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q4b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])


q5a <- plotTernaryGradient(tree1, bestmodel1$par[1:4], bestmodel1$par[5], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q5a$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])

q5b <- plotTernaryGradient(tree1, bestmodel1$par[6:9], bestmodel1$par[10], bins=4, line.weight=100, simulate.data=TRUE, simreps=100)
.simcounts <- round(10*log(table(factor(q5b$simdat$tip_states, 1:16))),0)
.simcols <- make.transparent(viridis::viridis(max(.simcounts))[.simcounts],0.2)
TernaryTiles(pi$centers[,1], pi$centers[,2], pi$centers[,3], 4, col=.simcols[(pi$labels)])
#dev.off()

fits$repmode_complete <- bestmodel0
fits$repmode_all <- bestmodel1

saveRDS(fits, "../output/fits_8132022.rds")
