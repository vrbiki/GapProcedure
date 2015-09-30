#' Extract clusters of a particular size
#'
#' @param classification.vector A classification vector.
#' @param small Specifies the minimum number of members the cluster contain
#' @param big The minimum number of members the cluster must contain
#' @return groups A vector of class labels that satisfy the requirements described in \code{Details}.
#' @return ind A vector of indices indicating which elements of \code{classification.vector} belong to one of the clusters stored in \code{groups}.
#' @details This function returns the cluster labels (and members) which have between \code{small} and \code{big} members
extractClusterMembers <-  function(classification.vector, small=NULL, big=NULL){
  # classification.vector=x; small = 2; big=5
  if (small> big)
    stop("small must be smaller than big")

  tfc <- table(classification.vector)
  (groupLabels <- as.numeric(names(tfc[tfc%in% small:(big-1)])))  # clusters that have between small to big members
  groupLabels <- setdiff(groupLabels, 0)                              # disclude the singleton group (group label of 0) if applicable
  seqIndex <- which(classification.vector %in% groupLabels)           # sequences belonging to on of the clusters stored in "groupLabels"
  return(list(ind=seqIndex, groups=groupLabels))
}
