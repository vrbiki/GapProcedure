#' Maps a probabilty matrix (in which each column sums to 1) to classification vector.
#'
#' @param z An $M$ by $N$ probability matrix.  Each column must sum to 1.
#' @param warn A logical indicating whether or not a warning should be displayed if any column indexes do not appear in the returned classification vector.
#' @return An integer vector of length $N$.  The $i$th element provides the row index ($1,2,\dots, M$) of the maximum value of the $i$th column of \code{z}.
#' @details Ties are resolved by randomly sampling a single number from the indices which attain a maximum.
#' @examples
#' set.seed(1234)
#' (ex.matrix <- matrix(sample(100,12), ncol=4, nrow=3))
#' (prob.matrix <- sweep(ex.matrix,2,colSums(ex.matrix),"/"))
#' (mapclass <- map2cluster(prob.matrix))
#'
#' # Create a new rows corresponding to the maximum value of each column
#' (ex.matrix2 <- rbind(ex.matrix, ex.matrix[cbind(mapclass,1:ncol(ex.matrix))]))
#' (prob.matrix2 <- sweep(ex.matrix2,2,colSums(ex.matrix2),"/"))
#'
#' # in the presence of a tie, sequences are assigned to the largest cluster
#' (map2cluster(prob.matrix2, random=TRUE))
#' (map2cluster(prob.matrix2, random=TRUE))

map2cluster <- function(z,warn=TRUE, random=FALSE){
  nobs <- ncol(z)
  zout <- numeric(nobs)
  G <- nrow(z)
  J <- 1:G

  if (!all.equal(colSums(z),rep(1,nobs)))
    stop("each column of `z` must sum to 1")

  if (!random){

    #zout <- apply(z,2,which.max)
    for (i in 1:nobs) {
      zcol <- z[,i]
      maxZ <- which(zcol==max(zcol))

      if (length(maxZ)!=1){
        zout[i] <- J[maxZ[which.max(rowSums(z)[maxZ])]]
      } else {
        # map observations to the cluster which has the highest fuzzyZ value
        zout[i] <- J[maxZ]
      }

    }

  } else {

    for (i in 1:nobs) {
      zcol <- z[,i]
      maxZ <- which(zcol==max(zcol))

      if (length(maxZ)!=1){
        zout[i] <- sample(maxZ,1)
      } else {
        # map observations to the cluster which has the highest fuzzyZ value
        zout[i] <- J[maxZ]
      }

    }

  }
  # test
  return(zout)
}
