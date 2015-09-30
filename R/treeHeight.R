#' Determines the height of a rooted tree
#'
#' @param rooted.tree A rooted tree in Newick format.
#' @return The height of a the inputted \code{rooted.tree}
treeHeight <- function(rooted.tree){
  if (!is.rooted(rooted.tree)) {
    stop("A rooted tree is required")
  }
  edge.from <- rooted.tree$edge[, 1]
  edge.to <- rooted.tree$edge[, 2]
  edge.id <- which(edge.from == (length(rooted.tree$tip.label) +
                                   1))
  n.from <- edge.to[edge.id[1]]
  ret <- rooted.tree$edge.length[edge.id[1]]
  repeat {
    edge.id <- which(edge.from == n.from)
    if (length(edge.id) == 2) {
      n.from <- edge.to[edge.id[1]]
      ret <- ret + rooted.tree$edge.length[edge.id[1]]
    }
    else {
      break
    }
  }
  ret
}
