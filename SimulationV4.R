################################################################################
rm(list = ls())
# for tree structure
library(ape)
# for multicore processing(optional)
library(parallel)
#set.seed(1)
################################################################################
#Parameter nCell: Number of Cells
#          nMutation: Number of Mutations postions           
#Returns   mutationCellMat: Matrix of Cells and the Mutations in each cell
#          treeReal: Randomly generated tree structure
#          mutationEdgeMat: Matrix of edges and the Mutations on each edge

genReal <- function(nCell, nMutation) {
  cellID <- paste0("C", 1:nCell)
  mutationID <- paste0("M", 1:nMutation)
  
  treeReal <- rtree(nCell, tip.label = cellID)
  
  dataEdge <- data.frame(from = treeReal$edge[, 1], to = treeReal$edge[, 2], length = treeReal$edge.length)
  nEdge <- nrow(dataEdge)
  rownames(dataEdge) <- edgeID <- paste0("E", 1:nEdge)
  
  edgeCellMat <- matrix(0L, nEdge, nCell, dimnames = list(edgeID, cellID))
  for(i in 1:nEdge){
    tempNode <- dataEdge$to[i]
    if(tempNode <= nCell){
      edgeCellMat[i, treeReal$tip.label[tempNode]] <- 1L
    }else {
      tempSubtree <- extract.clade(treeReal, tempNode)
      edgeCellMat[i, tempSubtree$tip.label] <- 1L
    }
  }
  
  mutationEdgeMat <- t(rmultinom(nMutation, 1, dataEdge$length)); dimnames(mutationEdgeMat) <- list(mutationID, edgeID)
  mutationCellMat <- matrix(0L, nMutation, nCell, dimnames = list(mutationID, cellID))
  for (i in 1:nEdge)
    mutationCellMat[mutationID[mutationEdgeMat[, i] == 1], cellID[edgeCellMat[i,] == 1]] <- 1L
  
  return (list(mutationCellMat, treeReal, mutationEdgeMat))
}
################################################################################
#Parameter mutationCellMat: Matrix of Cells and the Mutations in each cell
#          percentDrop: percent of cells that have allele drop out
#Returns   dropIndex: Cells that have allels dropped

funDropAllele <- function(mutationCellMat, percentDrop){
  dropStatus <- sample(c(1,0), dim(mutationCellMat)[2],prob = c(percentDrop, 1-percentDrop), replace = TRUE)
  dropIndex <- which(dropStatus ==1)
  return(dropIndex)
}
################################################################################
#Parameter realData: Matrix of Cells and the Mutations in each cell
#          seqDepth: average sequence depth
#          e: error rate
#Returns   dataObs: 3D array of cells, mutations and base pair readings(A,T,C,G)

genObs <- function(realData, seqDepth, e, dropCells) {
  nMutation <- dim(realData)[1]
  nCell <- dim(realData)[2]
  
  cellID <- paste0("C", 1:nCell)
  mutationID <- paste0("M", 1:nMutation)

  depthMat <- matrix(rpois(nMutation * nCell, seqDepth), nMutation, nCell,dimnames = list(mutationID, cellID))
  
  dataObs <-array(NA,dim = c(nMutation, nCell, 4), dimnames = list(mutationID, cellID, c("A", "T", "C", "G")))
  
  for (i in 1:nMutation) {
    for (j in 1:nCell) {
      tempDepth <- depthMat[i, j]
      #checks if allele is dropped
      if (j %in% dropCells ) {
        #check if it is mutated
        if( realData[i, j] == 1){
          #randomly choose an allele to drop
          if (sample(c(1, 0), 1) == 1) {
            tempProb <- c(e / 3, 1 - e, e / 3, e / 3)#TT
            }
          else{
            tempProb <- c(1 - e , e / 3, e / 3, e / 3) #AA
          }
        }
        else{
          tempProb <- c(1 - e , e / 3, e / 3, e / 3) #AA
        }
        tempDepth <-tempDepth/2
      }
      
      else if (realData[i, j] == 0) {
        tempProb <- c(1 - e , e / 3, e / 3, e / 3) #AA
      }else{
        tempProb <- c(0.5 - e / 3, 0.5 - e / 3, e / 3, e / 3) #AT
      }
      dataObs[i, j, ] = rmultinom(1, tempDepth, tempProb)
    }
  }

  return(dataObs)
}


################################################################################
#Parameter dataObs: 3D array of cells, mutations and base pair readings(A,T,C,G)
#          e: error rate
#Returns  estimateMutationCellMat: Matrix of Cells and the estimated Mutations in each cell

estMutation <- function (dataObs, e){
  AAProb <- c(1 - e, e / 3, e / 3, e / 3) #AA
  ACProb <- c(0.5 - e / 3, e / 3, 0.5 - e / 3, e / 3) #AC
  AGProb <- c(0.5 - e / 3, e / 3, e / 3, 0.5 - e / 3) #AG
  ATProb <- c(0.5 - e / 3, 0.5 - e / 3, e / 3, e / 3) #AT
  CCProb <- c(e / 3, e / 3, 1 - e, e / 3) #CC
  CGProb <- c(e / 3, e / 3, 0.5 - e / 3, 0.5 - e / 3) #CG
  CTProb <- c(e / 3, 0.5 - e / 3, 0.5 - e / 3, e / 3) #CT
  GGProb <- c(e / 3, e / 3, e / 3, 1 - e) #GG
  GTProb <- c(e / 3, 0.5 - e / 3, e / 3, 0.5 - e / 3) #GT
  TTProb <- c(e / 3, 1 - e, e / 3, e / 3) #TT
  probList <- list(AA = AAProb, AC = ACProb, AG = AGProb, AT = ATProb, CC = CCProb, CG = CGProb, CT = CTProb, GG = GGProb, GT = GTProb, TT = TTProb)

  
  funGenotypeProb <- function(nVec){
    genotypeProb <- rep(NA, length(probList)); names(genotypeProb) <- names(probList)
    for(i in 1:length(genotypeProb))genotypeProb[i] <- dmultinom(nVec, prob = probList[[i]])
    genotypeProb <- genotypeProb / sum(genotypeProb)
    return(genotypeProb)
  }
  
  funPriorProb <- function(count){
    countProb <- count/sum(count)
    p <- countProb[1] #A
    q <- countProb[2] #T
    r <- countProb[3] #C
    s <- countProb[4] #G
    priorMat <- c(p^2, 2*p*r, 2*p*s, 2*p*q, r^2, 2*r*s, 2*q*r, s^2, 2*q*s, q^2) ;names(priorMat) <- names(probList)
    return (priorMat)
  }
  
  SumMat <- (t(apply(dataObs, 1, colSums)))
  priorMat <-  apply(SumMat, 1, funPriorProb)
  priorArray <- aperm(replicate(dim(dataObs)[2], priorMat),c(2,3,1))
  
  postProbGenotype <- apply(dataObs, 1:2, funGenotypeProb); postProbGenotype <- aperm(postProbGenotype, c(2, 3, 1))
  #multicore apply 
  #postProbGenotype <- parApply(clusterl, dataObs, 1:2, funGenotypeProb); postProbGenotype <- aperm(postProbGenotype, c(2, 3, 1))
  
  #calculating prob with dropout
  alleleDropPrior <-.2 #change this to correct value
  postProbGenotype[,,4]<-alleleDropPrior*(.5*(postProbGenotype[,,1]+postProbGenotype[,,10]))+(1-alleleDropPrior)*(postProbGenotype[,,4])
  
  #Prior Calculation
  # postProbGenotype <- postProbGenotype*priorArray
   scale <- function(vec){return (vec/sum(vec))}
   postProbGenotype <- apply( postProbGenotype, 1:2, scale); postProbGenotype <- aperm(postProbGenotype, c(2, 3, 1))
  
  #multicore apply 
  # #postProbGenotype <- parApply(clusterl, postProbGenotype, 1:2, scale); postProbGenotype <- aperm(postProbGenotype, c(2, 3, 1))

  estimateMutationCellMat <- (apply(postProbGenotype, c(1, 2), which.max) != 1) + 0L
  return (estimateMutationCellMat)
}

###########################################################################
#Parameter t: tree structure
#          est: true if it is estimated tree: false if it is real tree
#Returns  edgeMat: Matrix of edges and cells under each edge

genEdgeMat <- function(t, est=FALSE){
  
  
  nEdge <- length(t$edge.length)
  nCell <- length(t$tip.label)
  
  edgeMat <- matrix (0,nEdge,nCell)
  if(est==TRUE){dimnames(edgeMat) <- list(paste0("Edge", 0:(nEdge-1)), paste0("C", 0:(nCell-1)))}
  else {dimnames(edgeMat) <- list(paste0("Edge", 1:nEdge), paste0("C", 1:nCell))}
  
  for(i in 1:nEdge){
    tempNode <- t$edge[i, 2]
    if(tempNode <= nCell){
      edgeMat[i, t$tip.label[tempNode]] <- 1L
    }else {
      tempSubtree <- extract.clade(t, tempNode)
      edgeMat[i, tempSubtree$tip.label] <- 1L
    }
  }
  
  return(edgeMat)
} 
###########################################################################
#Parameter nMat: observed data from one mutation postion
#          t: tree structure
#          real: true if it is real true, false if it is inferred tree
#          edgeCellMat: Matrix of edges and cells under each edge
#Returns  edgeProb: Posterior Probability of the likelihood of the mutation being on each edge 

funEdgeProb <- function(nMat, t, real =TRUE , edgeCellMat){
  if(real==TRUE){
    edgelength <- t$edge.length/sum(t$edge.length)
  }else{edgelength <- (t$edge.length)[-1]/sum((t$edge.length)[-1])}
  
  e=.01
  AAProb <- c(1 - e, e / 3, e / 3, e / 3) #AA
  ATProb <- c(0.5 - e / 3, 0.5 - e / 3, e / 3, e / 3) #AT
  
  
  cellNormProb <- c(apply(nMat, 1, dmultinom, prob = AAProb))
  cellMutationProb <- c(apply(nMat, 1, dmultinom, prob = ATProb))
  edgeProb <- rep(NA, nrow(edgeCellMat)); names(edgeProb) <- rownames(edgeCellMat)
  for(i in 1:nrow(edgeCellMat)){
    nMat0 <- cellNormProb[which(edgeCellMat[i, ] == 0)]
    nMat1 <- cellMutationProb[which(edgeCellMat[i, ] == 1)]
    
    edgeProb[i] <- prod(nMat0) * prod(nMat1)* edgelength[i] 
    
  }
  edgeProb <- edgeProb / sum(edgeProb)
  return(edgeProb)
}

###########################################################################
#Parameter real: real matrix of cells and mutations in each cell
#          est: estimated matrix of cells and mutations in each cell
#Returns   sensitivity: true postivesdectected / total true postives
#          FDR:  false postives / total postive dectection
funStat <- function(real, est){
  sensitivity <- 0
  FDR <- 0
  realTotal <- 0
  for(i in 1:nrow(est)){
    for(j in 1:ncol(est)){
      if(real[i,j] == 1)realTotal <- realTotal+1
      if(est[i,j] == 1 & real[i,j] == 1){
        sensitivity <- sensitivity+1
      } 
      else if(est[i,j] == 1 & real[i,j] !=1 ){
        FDR <- FDR +1
        
      }
      
    }
  }
  FDR <- FDR/(FDR+sensitivity)
  sensitivity <- sensitivity/realTotal
  return (list(sensitivity, FDR))
  
}

###########################################################################
#Parameter estMutationCellMat: matrix of estimated mutation status of each cell at each location
#Returns   estTree: inferred tree
genEstTree <- function(estMutationCellMat){
  estclusterData <- cbind("C0" = 0, estMutationCellMat)
  estimatedDist <- dist(t(estclusterData), method = "manhattan")
  estTree <- root(fastme.bal(estimatedDist), "C0")
  return(estTree)
}
###########################################################################
simMethod <-function(nCell = 200, nMutation = 1000, seqDepth = 50, e = .01)  {
  
  cellID <- paste0("C", 1:nCell)
  mutationID <- paste0("M", 1:nMutation)
  
  #generate real data
  realdata <- genReal(nCell, nMutation)
  mutationCellMat <- realdata[[1]]
  realTree <- realdata[[2]]
  mutationEdgeMat <- realdata[[3]]
  dropCells <- funDropAllele(mutationCellMat,.2)
  
  #generate observed data
  dataObs <- genObs(mutationCellMat, seqDepth , e, dropCells)
  
  #standard single sample anaylsis
  estMutationCellMat <- estMutation(dataObs, e)
  stat1 <- funStat(mutationCellMat, estMutationCellMat)

  #infer tree
  estTree <- genEstTree(estMutationCellMat)
  
  #view Trees
  #plot(realTree ,main = "Real Tree",  edge.width = 2, cex = 1, label.offset = 0.1)
  #plot(estTree ,main = "Estimated Tree",edge.width = 2, cex = 1, label.offset = 0.1)
  #edgelabels(cex =.7, c(0:38));nodelabels(cex = .8);tiplabels(cex = .8)
  
  #Joint Analysis Assuming Phylogenetic Tree Topology
  EdgeCellMat <- genEdgeMat(realTree)
  RpostProbEdge <- apply(dataObs, 1, funEdgeProb, t = realTree, edgeCellMat = EdgeCellMat); RpostProbEdge <- t(RpostProbEdge)
 
  RpostEdgeMat <- matrix(0, nMutation, ncol(RpostProbEdge))
  for (i in 1:nMutation) {
    RpostEdgeMat[i, which.max(RpostProbEdge[i, ])] <- 1L
  }
  RpostMutationMat <- matrix(0L, nMutation, nCell, dimnames = list(mutationID, cellID))
  for (i in 1:ncol(RpostProbEdge)) {
    RpostMutationMat[mutationID[RpostEdgeMat[, i] == 1], cellID[EdgeCellMat[i, ] == 1]] <- 1L
  }
  stat2 <- funStat(mutationCellMat, RpostMutationMat)
  
  #Joint Analysis Using Inferred Phylogenetic Tree
  estEdgeCellMat <- genEdgeMat(estTree, est = TRUE)[-1, -1]
  postProbEdge <- apply(dataObs, 1, funEdgeProb, t = estTree, edgeCellMat = estEdgeCellMat,real =FALSE); postProbEdge <- t(postProbEdge)
  
  postEdgeMat <- matrix(0, nMutation, ncol(postProbEdge))
  for (i in 1:nMutation) {
    postEdgeMat[i, which.max(postProbEdge[i, ])] <- 1L
  }
  postMutationMat <- matrix(0L, nMutation, nCell, dimnames = list(mutationID, cellID))
  for (i in 1:ncol(postProbEdge)) {
    postMutationMat[mutationID[postEdgeMat[, i] == 1], cellID[estEdgeCellMat[i, ] == 1]] <-
      1L
  }
  stat3 <- funStat(mutationCellMat, postMutationMat)
  
  return(list(stat1, stat2, stat3))
}
########################################################################
#multicore processing
# numCores <- detectCores()
# clusterl <- makeCluster(numCores)
# clusterEvalQ(clusterl,  {e <-.01
# AAProb <- c(1 - e, e / 3, e / 3, e / 3) #AA
# ACProb <- c(0.5 - e / 3, e / 3, 0.5 - e / 3, e / 3) #AC
# AGProb <- c(0.5 - e / 3, e / 3, e / 3, 0.5 - e / 3) #AG
# ATProb <- c(0.5 - e / 3, 0.5 - e / 3, e / 3, e / 3) #AT
# CCProb <- c(e / 3, e / 3, 1 - e, e / 3) #CC
# CGProb <- c(e / 3, e / 3, 0.5 - e / 3, 0.5 - e / 3) #CG
# CTProb <- c(e / 3, 0.5 - e / 3, 0.5 - e / 3, e / 3) #CT
# GGProb <- c(e / 3, e / 3, e / 3, 1 - e) #GG
# GTProb <- c(e / 3, 0.5 - e / 3, e / 3, 0.5 - e / 3) #GT
# TTProb <- c(e / 3, 1 - e, e / 3, e / 3) #TT
# probList <- list(AA = AAProb, AC = ACProb, AG = AGProb, AT = ATProb, CC = CCProb, CG = CGProb, CT = CTProb, GG = GGProb, GT = GTProb, TT = TTProb)})
########################################################################

r <- simMethod(200, 1000, 5, .01)


