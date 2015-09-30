#' Converts the matrix index to a vector index
#'
#' @param row.ind A integer indicating the row index of a matrix
#' @param col.ind An integer indiciation the column index of a matrix
#' @param mat.dim A vector of length 2 providing the dimensions of the matrix
#' @return out An interger indicating the scalar index of the matrix corresponding to row index `row.ind` and column index `col.ind`.
#' @examples
#' (mat <- matrix(1:25, nrow=5))
#' matrix2vecInd(1,3,dim(mat))
#'
#' (mat2 <- matrix(1:30, nrow=5))
#' matrix2vecInd(1,3,dim(mat2))
#'
#' (mat3 <- matrix(1:30, ncol=5))
#' matrix2vecInd(1,3,dim(mat3))
#'
#' matrix2vecInd(c(1,4,6),c(3,1,5),dim(mat3))
matrix2vecInd <- Vectorize(

  matrix2vec <- function(row.ind, col.ind, mat.dim){
  nr <- mat.dim[1]
  nc <- mat.dim[2]

  if (row.ind > nr)
    stop(c("\'row.ind\' must not exceed ",nr))
  if (col.ind > nc)
    stop(c("\'col.ind\' must not exceed ",nc))

  out <- nr*(col.ind-1) + row.ind

  return(out)}

  , vectorize.args=c("row.ind", "col.ind"))


# testmat <- mat3
# for (i in 1:length(testmat)){
#   for (rr in 1:nrow(testmat)){
#     for (cc in 1:ncol(testmat)){
#         sind <- testmat[rr,cc]
#         sind2 <- matrix2vecInd(rr,cc,dim(testmat))
#         if (sind!=sind2)
#           stop("they are not the same")
#     }
#   }
# }

