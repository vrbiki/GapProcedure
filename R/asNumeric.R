#' Maps DNAbin class to the numbers 0 through 16.
#'
#' @param DNAbin.data DNA sequence data in \code{DNAbin} (matrix) form.
#' @return A matrix (or list) of DNA sequences with code values of 0--16
#' @details The DNAbin class is mapped to the number 0 through 16.  Below is the scheme
asNumeric <- function(DNAbin.data){
  if (class(DNAbin.data)!="DNAbin")
    stop("input data must have class \'DNAbin\'")
  f <- function(xx) {
    ans <- .ns[match(as.numeric(xx), .bs)]
    if (is.matrix(xx)) {
      dim(ans) <- dim(xx)
      dimnames(ans) <- dimnames(xx)
    }
    ans
  }
  if (is.list(DNAbin.data)){
    lapply(DNAbin.data, f)
  } else {
    f(DNAbin.data)
  }
}
