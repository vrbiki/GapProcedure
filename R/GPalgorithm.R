GPalgorithm <- function(DNAbin.data, stage=1, dist.matrix=NULL, submod=NULL,percent){

  # DNAbin.data <- hivDNAbin
  if (class(DNAbin.data)!="DNAbin")
    stop("input sequence data must be in DNAbin form")
  if (!is.matrix(DNAbin.data))
    stop("input sequence data must be in matrix form (check that all sequences have the same length)")

  if (is.null(dim(DNAbin.data)))
    return(list(dist=0, classification=1, G=1))
  N <- nrow(DNAbin.data)

  # calculate the distance matrix if not already provided
  if (is.null(dist.matrix)){
    if (submod=="aJC69"){
      dist.time <- system.time(dist.matrix <- ajc69(DNAbin.data))
    } else if (submod=="aK80"){
      dist.time <- system.time(dist.matrix <- ak80(DNAbin.data))
    } else{
      stop("The substition model must be set to \'aJC69\' or \'aK80\'")
    }
  } else {
    dist.time <- 0
    if (!is.matrix(dist.matrix))
      stop("The specified dist matrix is not a matrix")
    if (!isSymmetric(dist.matrix))
      stop("The specified dist matrix is not a symmetric")
    if (any(dist.matrix<0))
      stop("The specified dist matrix is has negative values")
  }

  sortdist <- t(apply(dist.matrix, 1, sort))                    # sorts the rows of dist.matrix
  orddist <-  t(apply(dist.matrix, 1, order))                   # gives the sorted order of above
  half <- 1:(N*percent)                                         # restricts the search for the gap
  maxdifno0s <- apply(sortdist[,half],1,diffpossitive)          # gives the location of largest gap
  if (any(is.na(maxdifno0s)))
    return(list(dist=dist.matrix, gaps=maxdifno0s, prob=TRUE))

  Nmat <- isNeighbours(orddist, maxdifno0s)                # calculates the neighbours logical matrix
  Umat <- mgcv::uniquecombs(Nmat)                       # unique patterns
  uindex <- attr(Umat,"index")                             # mapping of obs to unique pattern
  Pmat <- Umat*as.vector(table(uindex))
  Zmat <- sweep(Pmat, 2,  colSums(Pmat), "/")
  # Zmat <- Zmat[order(rowSums(Umat),decreasing=TRUE),] # prop.matrix from old cold

  cls <- map2cluster(Zmat)

  gapout <- list(NeighourMatrix = Nmat, UniqueMatrix = Umat, PropMatrix = Pmat, IntegerCutoffs = maxdifno0s, DistCalcTime=dist.time)

  return(list(dist=dist.matrix, FreqMatrix = Zmat, classification=cls, G=length(unique(cls)), extra=gapout, prob=FALSE))

}

