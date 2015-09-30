#' A four group simulation.
#'
#' A dataset containing the class of 100 simulated DNA sequences of length 800.
#'
#' @docType data
#' @keywords datasets
#' @name simulation
#' @usage data(simulation)
#' @format A data frame with 100 sequences (rows) and 801 columns.  The first column \code{X.class} corresponds to the class label and the remaining columns (\code{X1}, ... \code{X800}) give the DNA nucleotide code for the corresponding site.
NULL


#' The Gap Procedure
#'
#' @description This function implements the Gap Procedure algorithm described in Vrbik et al.
#' for aligned DNA sequence data. The partition of the data is based on gaps in pairwise distances.
#' @param x Aligned DNA sequence data in \code{DNAbin} (matrix) form.
#' @param submod Specifies the evolutionary model to be used for calculating the pairwise genetic
#' distances between the sequences stored in \code{x}.  Options are: the "aJC", "aK80".
#' @param distance.matrix A pre-specified distance matrix used to determine the clusters within \code{x}.
#' @param outlier.adj Restricts the search for the larger gap to the first (100*\code{outlier.adj})\% of sorted pairwise distances.
#' @return
#' \describe{
#'   \item{dist}{A pairwise distance matrix calculated accorinding to \code{submod}.  Returns \code{distance.matrix} if this argument was specified.}
#'   \item{classification}{A vector indicating the clusters of DNA sequences found by the Gap Procedure.}
#'   \item{AlgorithmOutputs}{Supplementary outputs (e.g. matrices, CPU times and distances between gaps)
#'   calculated by the Gap Procedure algorithm.}
#' }
#' @details This function performs the Gap Procedure algorithm as described in Vrbik et al.  The input sequences must be aligned DNA sequences in \code{DNAbin} format.   If no distance matrix is specified in \code{distance.matrix}, this function computes a matrix of pairwise distances an ajusted version of the Kimura 1980 model of DNA evolution.  This adjustment takes into account ambiguous IUPAC codes. For details on these calculations see the Supplmentary Information in Vrbik et al.
#' @references Irene Vrbik, David A Stephens, Michel Roger and Bluma G Brenner. The Gap Procedure: for the identification of phylogenetic clusters in HIV-1 sequence data.s
#' @seealso \code{\link{ajc69},\link{ak80}}
#' @examples
#' library(ape)
#' data(woodmouse)
#' test <- GapProcedure(woodmouse,submod="aK80",outlier.adj=1)
#'
#' # to see the classifcaiton vector:
#' test$classification
#'
#' # to return the pairwise distance matrix:
#' test$dist
GapProcedure <- function(x, submod="aK80", distance.matrix=NULL, outlier.adj=0.9){

  if (class(x)!="DNAbin")
    stop("input sequence data must have class \'DNAbin\'")
  if (!is.matrix(x))
    stop("input sequence data must be in matrix form (check that all sequences have the same length)")
  if (!(submod %in% c("aJC69","aK80")))
    stop("\'submod\' must be set to: \'aJC69\', or \'aK80\'")

  GPout <- GPalgorithm(x, dist.matrix=distance.matrix, submod=submod,percent=outlier.adj)

  GPdist <- GPout$dist
  if (GPout$prob){
    warning("There is insufficient genetic diversity: No gaps found in the first ", round(nrow(x)*outlier.adj), " ordered pairwise distances")
    return(list(dist=GPdist, AlgorithmOutputs=GPout))
  }


  GPclass <- afan(GPout$classification)
  N <- dim(x)[1]
  GapSizes <- GPout$maxdiffs
  DistCutoffs <- GPout$distcutoff

  return(list(dist=GPdist, classification=GPclass,  AlgorithmOutputs=GPout))
}
