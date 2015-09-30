#' Maps singleton class labels to the class label \code{0}.
#'
#' @param x A vector of class labels.
#' @return A new (or identical) vector of class lables.
#' @details If there are no singleton (i.e., one member) clusters \code{rmsingletons} will return a vector identical to \code{x}.  Otherwise, all singleton class labels are replaced by zero.
#' @examples
#' rmsingletons(c(4,4,8,1,1,1,3,3,6))
#' rmsingletons(c(4,4,1,1,1,3,3,2,2))
rmsingletons <- function(cvec){
  tb <- table(cvec)
  if (length(which(tb==1))==0){
    # warning("there are no singletons")
    return(cvec)
  }
  clusters.with.one.member <- as.numeric(names(tb)[which(tb==1)])
  ind.in.singleton.cluster <- which(cvec%in%clusters.with.one.member)
  cnew <- cvec
  cnew[ind.in.singleton.cluster] <- 0
  return(cnew)
}
