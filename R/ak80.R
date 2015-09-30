#' Adjusted Kimura 1980 Pairwise Distances
#'
#' @param DNAbin.data Aligned DNA sequence data in \code{DNAbin} (matrix) form.
#' @return A pairwise distance matrix
#' @details This function computes a matrix of pairwise distances from aligned DNA sequences using an ajusted version of the Kimura 1980 model of DNA evolution.  This adjustment takes into account ambiguous IUPAC codes. For details on these calculations see the Supplmentary Information in Vrbik et al.
#' @references Irene Vrbik, David A Stephens, Michel Roger and Bluma G Brenner. The Gap Procedure: for the identification of phylogenetic clusters in HIV-1 sequence data.
#' @references Kiumra, M. \emph{A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences}. Journal of molecular evolution \strong{16}(2), 111-120 (1980).
#' @examples
#' library(ape)
#' data(woodmouse)
#' ak80(woodmouse)
ak80 <- function(DNAbin.data){
  if (class(DNAbin.data)!="DNAbin")
    stop("input sequence data must be in DNAbin form")
  if (!is.matrix(DNAbin.data))
    stop("input sequence data must be in matrix form (check that all sequences have the same length)")
  nx <- asNumeric(DNAbin.data)
  out <- AdjK80(nx, .TSmat, .TVmat)
  if (any(is.nan(out)))
    return("NaNs produced in computing H80 distance matrix")
  if (!all(is.finite(out)))
    return("Non-finite distances produced in computing H80 distance matrix")

  return(out)
}



