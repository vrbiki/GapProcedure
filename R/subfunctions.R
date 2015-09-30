#' Converts class labels to the lowest set of possible integers.
#'
#' @param x A vector of class labels.
#' @return A new (or identical) vector of class lables.
#' @details If the largest class label in \code{x} is equal to the number of groups and the smallest group label is equal to 1, the class vector returned by \code{afan} will be identical.  If however, potential lower class labels are skiped, a new vector of class labels is returned.  This vector will ensure that no potential class labels have been skipped.
#' @examples
#' afan(c(4,4,8,1,1,1,3,3,6))
#' afan(c(4,4,1,1,1,3,3,2,2))
afan <- function(x){
  if (!is.vector(x))
    stop("x must be a vector")
  out <- as.numeric(as.factor(x))
  return(out)
}

#' Returns the length of the unique elements in a vector x
#'
#' @param x A vector
#' @return An interger indicating the number of unique elements in x
lu <- function(x){length(unique(x))}

#' Replaces singletons (i.e., elements with value 0) with a unqiue positive valued class labels
#'
#' @param class.vec A vector of class labels
#' @return A new (or identical) class label with stictly positive class labels.
renameClusters <- function(class.vec){  # class.vec <- classlab
  out <- rep(NA, length(class.vec))
  zeros <- which(class.vec==0)
  nonzeros <- which(class.vec!=0)
  out[zeros] <- 0
  out[nonzeros] <- afan(class.vec[nonzeros])
  return(out)
}



#' Sorts a vector and computers the difference between sucessive numbers
#'
#' @param x A vector of non-negative values.
#' @return An integer indicating where largest difference in the sorted \code{x} vector is obsvered.
#' @details Note that (1) comparisons with 0 are not considered for the gap search (2) in the presence of a tie, the first large gap will be returned (see examples).
#' @examples
#' # Consider this vector of sorted values
#' x <- c(0,3,6.1,7,8)
#'
#' # The largest gap in x occurs between the second and third value.
#' diffpossitive(x)
#'
#' # This value should not change after we've shuffled x
#' sx <- sample(x)
#' diffpossitive(x)
#' cat("The gap occurs between sorted position", diffpossitive(x),"and", diffpossitive(x)+1)
#'
#' # Note that ties will return the location of the first gap
#' diffpossitive(c(1,2,4,8,9,10,14))
#' cat("The gap occurs between sorted position", diffpossitive(x),"and", diffpossitive(x)+1)
diffpossitive <- function(x, needtosort=FALSE,value=FALSE){ #
  if (needtosort) x <- sort(x)
  if (any(x<0))
    stop("the values of x must be non-negative")
  zeros<- which(x==0)
  if (length(zeros)==0){
    wm <- which.max(diff(x))
  } else {
    if (length(zeros)>=(length(x) - 1))
      return(NA)
    wm <- which.max(diff(x[-zeros]))
    wm <- wm + length(zeros)
  }
  if (value==TRUE){
    return(diff(x)[wm])
  } else {
    return(wm)
  }
}


#' Counts the number of bits in the A G C T position in the \code{DNAbin} bit-level coding scheme
#'
#' @param xL A DNA sequence in integer form
#' @param shift A logical indicating whether an initial shift needs to be made
#' @return A vector of integer values indicating the number of bits in position A G C T at each site of \code{xL}
#' @examples
#' # IUPAC code "A" has Bit-level code 1000 1000 which corresponds to the number 136
#' # there is 1 bit in the first four positions (corresponding to the binary coding system for A G C T  in the DNAbin class)
#' bit4Count(136)
#'
#' # IUPAC code "B" (c, g or t) has Bit-level code 0111 0000 which corresponds to the number 112
#' # there are 3 bits in the first four positions (corresponding to the binary coding system for A G C T  in the DNAbin class)
#' bit4Count(112)
bit4Count <- function(xL, shift=TRUE) {
  # xL <- 136
  if (shift) xL <- bitwShiftR(xL, 4)      # shift bits, to the right 4 times
  count = 0L;                             # stores the counts as integers
  while (xL > 0) {                        # until all bits are zero
    if (bitwAnd(xL, 1) == 1)              # check lower bit
      count <- count + 1L
    xL <- bitwShiftR(xL, 1)              # shift bits (1 Right) removes lower bit
  }
  return(count)
}

#' Counts the number of posible pairwise matches between \code{xL} and \code{yL},
#'
#' @param xL A DNA sequence in integer form.
#' @param yL A DNA sequence in integer form.
#' @return A vector of integer values indicating the number of matching character comparisons between the IUPAC symbols \code{xL} and \code{yL}.
#' @examples
#' # IUPAC code "V" (a, c or g) has Bit-level code 1110 0000 which corresponds to the number 224
#' # IUPAC code "H" (a, c or t) has Bit-level code 1011 0000 which corresponds to the number 176
#' # eg. V (a, c or g) and H (a, c or t) has 9 pairwise comparisons aa, ac, at, ca, cc, ct, ga, gc, gt
#' # and has 2 matching comparison (aa, cc).
#' bitMatch(224,176)
bitMatch <- function(xL,yL){
  xL <- bitwShiftR(xL, 4)     # shift bits, to the right 4 times
  yL <- bitwShiftR(yL, 4)     # shift bits, to the right 4 times
  bit4Count(bitwAnd(xL, yL),shift=FALSE)
}

vecbitMatch <- Vectorize(bitMatch)


#' Prints a vector separataed by commas
#'
#' @param x A vector
#' @return A character strings without quotes (class \code{noquote})
#' @examples
#' pcsv(c(3,12,4,5))
pcsv <- function(x){
  if (!is.vector(x))
    stop("x must be a vector")
  return(noquote(paste(x, sep="' '", collapse=",")))
}





