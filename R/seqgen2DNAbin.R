#' Converts the outputed sequences from \code{seq.gen} to `DNAbin` format
#'
#' @param read.seqgen.output The output from running seq.gen
#' @return The input sequences in DNAbin format
seqgen2DNAbin <- function(read.seqgen.output){
# read.seqgen.output: "code.type" "info"      "nseq"      "seqlen"    "seqname"   "org.code"  "org"       "byrow"
x <- read.seqgen(read.seqgen.output)
if (!(x$code.type == "NUCLEOTIDE"))
  stop("read.seqgen.output must be DNA nucleotide type")
N <- x$nseq
L <- x$seqlen
return(as.DNAbin(as.alignment(x$org.code)))
}
