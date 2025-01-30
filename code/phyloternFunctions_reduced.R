
TernaryBin <- function(xyz, bins){
  bins2 <- bins^2
  bins4 <- 2*bins2
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="") #Give the cells a state label  
  ddat <- apply(xyz, 2, cut,breaks=seq(0,1, length.out=bins+1), labels=1:bins) #Find where the data falls in the state space
  ddat <- match(apply(ddat, 1, paste, collapse=""), cellnames) #match the data's state to the state label
  ddat 
}

count_ambiguous <- function(vec, nstates){
  vals <- lapply(vec, function(x) strsplit(x, "&")[[1]])
  counts <- lapply(vals, function(x) table(factor(x, levels=1:nstates)))
  counts <- lapply(counts, function(x) x/sum(x))
  counts <- do.call(rbind, counts)
  counts <- apply(counts, 2, sum)
  return(counts)
}

r2Q <- function(r,c, bins){
  r <- c(r[1:2], 0, r[3])
  bins2 <- bins^2
  bins4 <- 2*bins2
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="") #Give the cells a state label  
  Q <- as.matrix(dist(cells))
  Q[Q!=1] <- 0
  non0rates <- which(Q != 0, arr.ind = T)
  tQ <- Q
  centers <- t(Ternary:::TriangleCentres(bins, direction=1))
  cells <-t(Ternary:::XYToTernary(centers[,1], centers[,2], direction=1))
  labels <- TernaryBin(cells, bins)
  centers <- centers[order(labels),]
  midcells <- cells[order(labels),]
  abscenter <- c(1/3, 1/3, 1/3)
  for(i in 1:nrow(non0rates)){
    d2 <- sum((midcells[non0rates[i,2], ] - abscenter)^2)
    d1 <- sum((midcells[non0rates[i,1], ] - abscenter)^2)
    tQ[non0rates[i,1], non0rates[i,2]] <- c*exp(sum(r[1:3]*(midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ])) + (d2-d1)*r[4])
  }
  diag(tQ) <- -1*apply(tQ, 1, sum)
  return(tQ)
}


simulateTernaryData <- function(tree, r, c, bins, line.weight=10, plot.points=TRUE, ...){
  r <- c(r[1:2], 0, r[3])
  bins2 <- bins^2
  bins4 <- 2*bins2
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="") #Give the cells a state label  
  Q <- as.matrix(dist(cells))
  Q[Q!=1] <- 0
  non0rates <- which(Q != 0, arr.ind = T)
  tQ <- Q
  centers <- t(Ternary:::TriangleCentres(bins))
  cells <-t(Ternary:::XYToTernary(centers[,1], centers[,2]))
  labels <- TernaryBin(cells, bins)
  
  #Ternary::TernaryPlot(alab="A", blab="B", clab="C", grid.lines = 3, grid.minor.lines = 0) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
  #points(centers, pch=as.character(labels))
  centers <- centers[order(labels),]
  midcells <- cells[order(labels),]
  abscenter <- c(1/3, 1/3, 1/3)
  Ternary::TernaryPlot(alab="A", blab="B", clab="C", grid.lines = bins, grid.minor.lines = 0,...) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
  for(i in 1:nrow(non0rates)){
    #Q[non0rates[i,1], non0rates[i,2]] <- (
    #  r[4]^(sum((midcells[non0rates[i,2], ] - abscenter)^2) - sum((midcells[non0rates[i,1], ] - abscenter)^2))*
    #  sum(r[1:3]^(bins*(midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ]))*abs((midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ])))*
    #  c
    #  )
    d2 <- sum((midcells[non0rates[i,2], ] - abscenter)^2)
    d1 <- sum((midcells[non0rates[i,1], ] - abscenter)^2)
    tQ[non0rates[i,1], non0rates[i,2]] <- c*exp(sum(r[1:3]*(midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ])) + (d2-d1)*r[4])
    arrows(centers[non0rates[i,1], 1],centers[non0rates[i,1], 2],centers[non0rates[i,2],1], centers[non0rates[i,2],2], lwd =  tQ[non0rates[i,1], non0rates[i,2]]*line.weight, length = 0.1)
  }
  
  diag(tQ) <- -1*apply(tQ, 1, sum)
  simdat <- castor::simulate_mk_model(tree, Q=tQ)
  if(plot.points){
    TernaryPoints(midcells[simdat$tip_states, ]+rnorm(3*length(simdat$tip_states),0, 0.05), pch=21, bg=bayou:::makeTransparent("red", 10), col=bayou:::makeTransparent("red", 10)) 
  }
  
  points(centers, pch=as.character(1:bins2), col="red")
  return(list(simdat=simdat, Q=tQ, indQ=Q))
}

plotTernaryIndex <- function(bins, levels=1,...){
  bins2 <- bins^2
  bins4 <- 2*bins2
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="")
  centers <- t(Ternary:::TriangleCentres(bins, direction=1))
  cells <-t(Ternary:::XYToTernary(centers[,1], centers[,2], direction=1))
  labels <- TernaryBin(cells, bins)
  for(i in 1:levels){
    Ternary::TernaryPlot(alab="Adult Mortality", blab="Juvenile Mortality", clab="Body Size", grid.lines = bins, grid.minor.lines = 0) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
    text(centers, labels=c(labels, bins2+labels)[(1+(i-1)*bins2):((i-1)*bins2+bins2)], ...)
  }
  return(list(centers=centers, cells=cells, labels=labels))
}

get_simmap_tips <- function(simmap){
  ntips <- length(simmap$tip.label)
  tip_index <- match(1:ntips, simmap$edge[,2])
  state <- sapply(simmap$maps[tip_index], function(x) names(x)[length(x)])
  state <- setNames(state, simmap$tip.label)
  return(state)
}

plotTernaryCounts <- function(.dat, bins, levels=1, offset=0, add=FALSE,...){
  bins2 <- bins^2
  bins4 <- 2*bins2
  if(class(.dat)[1]=="simmap"){
    tiprecon <- get_simmap_tips(.dat)
    counts <- table(factor(tiprecon, levels=1:bins4))
  } else {
    counts <- .dat
  }
  
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="")
  centers <- t(Ternary:::TriangleCentres(bins, direction=1))
  cells <-t(Ternary:::XYToTernary(centers[,1], centers[,2], direction=1))
  labels <- TernaryBin(cells, bins)
  print(counts)
  for(i in 1:levels){
    if(add==FALSE) Ternary::TernaryPlot(alab="Adult Mortality", blab="Juvenile Mortality", clab="Body Size", grid.lines = bins, grid.minor.lines = 0) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
    text(centers[order(labels),]+offset, labels=counts[(1+(i-1)*bins2):((i-1)*bins2+bins2)], ...)
  }
  
}

plotTernaryGradient <- function(tree, r, c, bins, line.weight=10,simulate.data=FALSE, simreps=1, ...){
  r <- c(r[1:2], 0, r[3])
  bins2 <- bins^2
  bins4 <- 2*bins2
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="") #Give the cells a state label  
  Q <- as.matrix(dist(cells))
  Q[Q!=1] <- 0
  non0rates <- which(Q != 0, arr.ind = T)
  tQ <- Q
  centers <- t(Ternary:::TriangleCentres(bins, direction=1))
  cells <-t(Ternary:::XYToTernary(centers[,1], centers[,2], direction=1))
  labels <- TernaryBin(cells, bins)
  centers <- centers[order(labels),]
  midcells <- cells[order(labels),]
  abscenter <- c(1/3, 1/3, 1/3)
  Ternary::TernaryPlot(alab="Adult Mortality", blab="Juvenile Mortality", clab="Body Size", grid.lines = bins, grid.minor.lines = 0,...) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
  
  #Ternary::TernaryPlot(alab="A", blab="B", clab="C", grid.lines = 3, grid.minor.lines = 0) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
  #points(centers, pch=as.character(labels))
  #centers <- centers[order(labels),]
  #midcells <- cells[order(labels),]
  #abscenter <- c(1/3, 1/3, 1/3)
  #Ternary::TernaryPlot(alab="Adult Mortality", blab="Juvenile Mortality", clab="Body Size", grid.lines = bins, grid.minor.lines = 0,...) #Plot with points labeled by their numeric state value (values with more than two digits won't show up right)
  for(i in 1:nrow(non0rates)){
    #Q[non0rates[i,1], non0rates[i,2]] <- (
    #  r[4]^(sum((midcells[non0rates[i,2], ] - abscenter)^2) - sum((midcells[non0rates[i,1], ] - abscenter)^2))*
    #  sum(r[1:3]^(bins*(midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ]))*abs((midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ])))*
    #  c
    #  )
    d2 <- sum((midcells[non0rates[i,2], ] - abscenter)^2)
    d1 <- sum((midcells[non0rates[i,1], ] - abscenter)^2)
    tQ[non0rates[i,1], non0rates[i,2]] <- c*exp(sum(r[1:3]*(midcells[non0rates[i,2], ] - midcells[non0rates[i,1], ])) + (d2-d1)*r[4])
    arrows(centers[non0rates[i,1], 1],centers[non0rates[i,1], 2],centers[non0rates[i,2],1], centers[non0rates[i,2],2], lwd =  tQ[non0rates[i,1], non0rates[i,2]]*line.weight, length = 0.1)
  }
  
  diag(tQ) <- -1*apply(tQ, 1, sum)
  if(simulate.data){
    simdat <- castor::simulate_mk_model(tree, Q=tQ, Nsimulations = simreps)
    #points(centers, pch=as.character(1:bins2), col="red")
    return(list(simdat=simdat, Q=tQ, indQ=Q))
  } else {
    return(list(Q=tQ, indQ=Q))
  }

  
}

make.lnLTernaryGrad <- function(tree, dat, bins){
  bins <- bins
  bins2 <- bins^2
  lnLTernaryGrad <- function(pp){
    Q <- r2Q(pp[1:3], pp[4], bins=bins)
    fitdat <- dat
    Ntips <- length(tree$tip.label)
    Nnodes <- tree$Nnode
    Nedges <- length(tree$edge.length)
    Nstates <- bins2
    tip_priors <- matrix(1e-08/(Nstates - 1), nrow=length(tree$tip.label), ncol=bins2)
    for(i in 1:nrow(tip_priors)){
      if(class(fitdat)=="character"){
        .states <- as.numeric(strsplit(fitdat[i], "&")[[1]])
      } else{
        .states <- fitdat[i]
      }
      dens <- (1 - 1e-08)/length(.states)
      tip_priors[i, .states] <- dens
    }
    tree_edge <- as.vector(t(tree$edge)) - 1
    edge.length <- tree$edge.length
    transition_matrix <- as.vector(t(Q))
    prior_probabilities_per_tip <- as.vector(t(tip_priors))
    root_prior_type = "max_likelihood"
    root_prior_probabilities = numeric(0)
    oldest_age <- -1
    runtime_out_seconds <- 0
    exponentiation_accuracy <- 0.001
    max_polynomials <- 1000
    lnL0 <- castor:::Mk_loglikelihood_CPP(Ntips=Ntips, Nnodes=Nnodes, Nstates=Nstates, Nedges=Nedges, prior_probabilities_per_tip = prior_probabilities_per_tip, 
                                          root_prior_type = root_prior_type, tree_edge <- tree_edge, edge_length=edge.length, transition_matrix=transition_matrix,
                                          root_prior=root_prior_probabilities, oldest_age=oldest_age, runtime_out_seconds = runtime_out_seconds, 
                                          exponentiation_accuracy = exponentiation_accuracy, max_polynomials = max_polynomials)
    return(lnL0$loglikelihood*-1)
  }
  return(lnLTernaryGrad)
}

r2Q_2x <- function(pp, bins, model="ER"){
  bins2 <- bins^2
  bins4 <- 2*bins2
  Q1 <- r2Q(pp[1:3], pp[4], bins=bins)
  Q2 <- r2Q(pp[5:7], pp[8], bins=bins)
  DM1 <- DM2 <- matrix(0, nrow=dim(Q1)[1], ncol=dim(Q1)[2])
      if(model=="ER"){
        diag(DM1) <- pp[9]
        diag(DM2) <- pp[9]
      } 
      if(model=="ARD"){
        diag(DM1) <- pp[9]
        diag(DM2) <- pp[10]
      }
      if(model=="Dollo"){
        diag(DM1) <- pp[9]
        diag(DM2) <- 0
      }
  Q <- rbind(cbind(Q1, DM1), cbind(DM2, Q2))
  diag(Q) <- 0
  diag(Q) <- -1*rowSums(Q)
  rownames(Q) <- colnames(Q) <- 1:bins4
  return(Q)
}

get_tippriors <- function(dat,bins=4, layers=1){
  bins2 <- bins^2
  bins4 <- 2*bins2
  Nstates <- layers*bins2
  tip_priors <- matrix(1e-08/(Nstates - 1), nrow=length(dat), ncol=layers*bins2)
  for(i in 1:nrow(tip_priors)){
    if(class(dat)=="character"){
      .states <- as.numeric(strsplit(dat[i], "&")[[1]])
    } else{
      .states <- dat[i]
    }
    dens <- (1 - 1e-08)/length(.states)
    tip_priors[i, .states] <- dens
  }
  colnames(tip_priors) <- 1:Nstates
  return(tip_priors)
}

make.lnL2xTernaryGrad <- function(tree, dat, bins, model="ER"){
  bins2 <- bins^2
  bins4 <- 2*bins2
  cache <- list(tree=tree, dat=dat, Ntips=length(tree$tip.label), Nnodes=tree$Nnode, Nedges=length(tree$edge.length), Nstates=bins4, 
                tree_edge=as.vector(t(tree$edge)) - 1, edge_length=tree$edge.length, bins=bins, bins2=bins2, bins4=bins4)
  tip_priors <- matrix(1e-08/(cache$Nstates - 1), nrow=length(tree$tip.label), ncol=bins4)
  for(i in 1:nrow(tip_priors)){
    if(class(cache$dat)=="character"){
      .states <- as.numeric(strsplit(cache$dat[i], "&")[[1]])
    } else{
      .states <- cache$dat[i]
    }
    dens <- (1 - 1e-08)/length(.states)
    tip_priors[i, .states] <- dens
  }
  cache$tip_priors <- tip_priors
  cache$prior_probabilities_per_tip <- as.vector(t(tip_priors))
  
  biTernaryGrad <- function(pp){
    Q1 <- r2Q(pp[1:3], pp[4], bins=cache$bins)
    Q2 <- r2Q(pp[5:7], pp[8], bins=cache$bins)
    DM1 <- DM2 <- matrix(0, nrow=dim(Q1)[1], ncol=dim(Q1)[2])
#    if(model=="ER"){
#      diag(DM1) <- pp[11]
#      diag(DM2) <- pp[11]
#    } 
##    if(model=="ARD"){
      diag(DM1) <- pp[9]
      diag(DM2) <- pp[10]
#    }
#    if(model=="Dollo"){
#      diag(DM1) <- pp[11]
#      diag(DM2) <- 0
#    }
    Q <- rbind(cbind(Q1, DM1), cbind(DM2, Q2))
    diag(Q) <- 0
    diag(Q) <- -1*rowSums(Q)
    #Ntips <- length(tree$tip.label)
    #Nnodes <- tree$Nnode
    #Nedges <- length(tree$edge.length)
    #Nstates <- bins2*2
    #tree_edge <- as.vector(t(tree$edge)) - 1
    #edge.length <- tree$edge.length
    .transition_matrix <- as.vector(t(Q))
    #prior_probabilities_per_tip <- as.vector(t(tip_priors))
    #root_prior_type = "max_likelihood"
    #root_prior_probabilities = numeric(0)
    #oldest_age <- -1
    #runtime_out_seconds <- 0
    #exponentiation_accuracy <- 0.001
    #max_polynomials <- 1000
    lnL0 <- castor:::Mk_loglikelihood_CPP(Ntips=cache$Ntips, Nnodes=cache$Nnode, Nstates=cache$Nstates, Nedges=cache$Nedges, prior_probabilities_per_tip = cache$prior_probabilities_per_tip, 
                                          root_prior_type = "max_likelihood", tree_edge=cache$tree_edge, edge_length=cache$edge_length, transition_matrix=.transition_matrix,
                                          root_prior=numeric(0), oldest_age=-1, runtime_out_seconds = 0, 
                                          exponentiation_accuracy = 0.001, max_polynomials = 1000)
    return(lnL0$loglikelihood*-1)
  }
  return(biTernaryGrad)
}



make.lnL2xIdenticalTernaryGrad <- function(tree, dat, bins, model="ER"){
  bins2 <- bins^2
  bins4 <- 2*bins2
  cache <- list(tree=tree, dat=dat, Ntips=length(tree$tip.label), Nnodes=tree$Nnode, Nedges=length(tree$edge.length), Nstates=bins4, 
                tree_edge=as.vector(t(tree$edge)) - 1, edge_length=tree$edge.length, bins=bins, bins2=bins2, bins4=bins4)
  tip_priors <- matrix(1e-08/(cache$Nstates - 1), nrow=length(tree$tip.label), ncol=bins4)
  for(i in 1:nrow(tip_priors)){
    if(class(cache$dat)=="character"){
      .states <- as.numeric(strsplit(cache$dat[i], "&")[[1]])
    } else{
      .states <- cache$dat[i]
    }
    dens <- (1 - 1e-08)/length(.states)
    tip_priors[i, .states] <- dens
  }
  cache$tip_priors <- tip_priors
  cache$prior_probabilities_per_tip <- as.vector(t(tip_priors))
  
  biTernaryGrad <- function(pp){
    Q1 <- r2Q(pp[1:3], pp[4], bins=cache$bins)
    Q2 <- r2Q(pp[1:3], pp[4], bins=cache$bins)
    DM1 <- DM2 <- matrix(0, nrow=dim(Q1)[1], ncol=dim(Q1)[2])
    #    if(model=="ER"){
    #      diag(DM1) <- pp[11]
    #      diag(DM2) <- pp[11]
    #    } 
    ##    if(model=="ARD"){
    diag(DM1) <- pp[5]
    diag(DM2) <- pp[6]
    #    }
    #    if(model=="Dollo"){
    #      diag(DM1) <- pp[11]
    #      diag(DM2) <- 0
    #    }
    Q <- rbind(cbind(Q1, DM1), cbind(DM2, Q2))
    diag(Q) <- 0
    diag(Q) <- -1*rowSums(Q)
    #Ntips <- length(tree$tip.label)
    #Nnodes <- tree$Nnode
    #Nedges <- length(tree$edge.length)
    #Nstates <- bins2*2
    #tree_edge <- as.vector(t(tree$edge)) - 1
    #edge.length <- tree$edge.length
    .transition_matrix <- as.vector(t(Q))
    #prior_probabilities_per_tip <- as.vector(t(tip_priors))
    #root_prior_type = "max_likelihood"
    #root_prior_probabilities = numeric(0)
    #oldest_age <- -1
    #runtime_out_seconds <- 0
    #exponentiation_accuracy <- 0.001
    #max_polynomials <- 1000
    lnL0 <- castor:::Mk_loglikelihood_CPP(Ntips=cache$Ntips, Nnodes=cache$Nnode, Nstates=cache$Nstates, Nedges=cache$Nedges, prior_probabilities_per_tip = cache$prior_probabilities_per_tip, 
                                          root_prior_type = "max_likelihood", tree_edge=cache$tree_edge, edge_length=cache$edge_length, transition_matrix=.transition_matrix,
                                          root_prior=numeric(0), oldest_age=-1, runtime_out_seconds = 0, 
                                          exponentiation_accuracy = 0.001, max_polynomials = 1000)
    return(lnL0$loglikelihood*-1)
  }
  return(biTernaryGrad)
}



recodeTraits <- function(dat, columns, bins){
  
  .dat <- dat[,columns]
  
  bins2 <- bins^2
  bins4 <- 2*bins2
  
  combos <- combinat::combn(rep(1:bins,3), 3) ## Create all combinations of values 1:bins
  sums <- apply(combos, 2, sum) #Sum the 3 values for the 3 traits
  combos <- combos[,sums %in% c(bins+2, bins+1)] #Only retain those that sum to bins+2 and bins+1 (fit in the triangle)
  cells <- t(combos[,!duplicated(apply(combos, 2,paste, collapse=""))]) #remove duplicated bin combinations
  colnames(cells) <- c("A", "B", "C") #label traits
  cells <- cells[order(cells[,3],(bins+1)-cells[,1],(bins+1)-cells[,2], decreasing=TRUE),] #reorder the cells to be prettier
  cellnames <- apply(cells, 1, paste, collapse="") #Give the cells a state label  
  ddat <- apply(.dat, 2, cut,breaks=seq(0,1, length.out=bins+1), labels=1:bins) #Find where the data falls in the state space
  ddat[is.na(ddat)] <- "."
  ddat <- apply(ddat, 1, paste, collapse="")
  ddat <- lapply(1:length(ddat), function(x) grep(ddat[x], cellnames)) #match the data's state to the state label
  ddat <- sapply(ddat, paste, collapse="&")
  dat <- mutate(dat, ddat=ddat) 
  return(dat)
}

layerTrait <- function(ddat, xlayer, bins){
   bins2 <- bins^2
   bins4 <- 2*bins2
   splitdat <- strsplit(ddat, "&")
   for(i in 1:length(xlayer)){
     if(is.na(xlayer[i])){
       splitdat[[i]] <- c(as.numeric(splitdat[[i]]), as.numeric(splitdat[[i]])+bins2)
     } else {
       splitdat[[i]] <- as.numeric(splitdat[[i]])+xlayer[i]*bins2
     }
   }
   sapply(splitdat, paste, collapse="&")
}
