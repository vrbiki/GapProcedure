#' Plots the sorted pairwise distances and indicates the position of the largest gap.
#'
#' @param GPout The output object created by the \code{GapProcedure} function.
#' @param classvec A vector indicating the class memberships of the data used to produce \code{GPout}.  The default (NULL) takes the classification vector stored in \code{GPout}.
#' @param Nindex A vector of sequence indices. If set to NULL (default), then the the first sequence in each class labels is considered.
#' @param show.singletons A logical indicating if singleton clusters should be plotted (default=\code{TRUE})
#' @param hcutoff A scalar indicating where a horizontal line should be drawn
#' @param savefile A character vector.  If set to \code{NULL} (default), the figures are not saved and ploted in R, otherwise the figure will be saved to \code{<savefile>Obs<s>Cluster<g>.pdf}, where \code{s} is the sequence index and \code{g} denotes the cluster number.
#' @param ... Options for \code{pdf()}
#' @details For each class label ($g = 1, \dots, G$), a side-by-side boxplot is created.
#' @examples
#' data(simulation)
#' cls <- simulation[,1]
#' dat <- as.DNAbin(as.alignment(as.matrix(simulation[,-1])))
#'
#' GP <- GapProcedure(dat)
#' GP$classification
#'
#' boxplotWB(GP)
#' GapPlot(GP)
GapPlot <- function(GPout, classvec=NULL, Nindex=NULL,  hcutoff = NULL,show.singletons=TRUE, savefile=NULL,...){

  if (is.null(classvec))
    classvec <- GPout$classification
  if (!show.singletons)
    classvec <- rmsingletons(classvec)
  distmat <- GPout$dist
  N <- ncol(distmat)
  if (is.null(Nindex)){
    first.ind.in.cluster <- match(setdiff(sort(unique(classvec)),0),classvec)
  } else {
    first.ind.in.cluster <- Nindex
  }


  distmat <- GPout$dist
  sortdist <- t(apply(distmat, 1, sort))                     # sorts the rows of dist.matrix
  orddist <-  t(apply(distmat, 1, order))                    # gives the sorted order of above
  half <- 1:(N*0.9)                                           # restricts the search for the gap
  ci.ind <- apply(sortdist[,half],1,diffpossitive)            # gives the location of largest gap
  ci <- apply(sortdist[,half],1, function(x) diffpossitive(x, value=TRUE)) # value of index above
  di <- sortdist[cbind(1:nrow(sortdist), ci.ind)]



  for (g in seq_along(first.ind.in.cluster)){ #  g=2
    ii <- first.ind.in.cluster[g]
    if (!is.null(savefile))
      pdf(paste(savefile,"Obs",ii, "Cluster", classvec[ii],".pdf",sep=""),...)

    #---------------------------------------------------
    #plot(sortdist[ii,], ylab= expression(d(X[i], X["[j]"])), xlab="j", cex.lab=1, main=paste("Cluster",g, "Sequence", ii))
    plot(sortdist[ii,], ylab= expression(d(X[i], X["[j]"])), xlab="j", col=0, main=paste("Cluster",g, "Sequence", ii))
    nn <- seq(ci.ind[ii])     # neighbours (nij = 1)
    nnn <- (ci.ind[ii]+1):N   # not neighbours (nij = 0)
    text(x=nn, sortdist[ii,][nn],  labels="N", cex=0.75) # labels=as.character(orddist[ii,][nn]),
    points(x=nnn, y=sortdist[ii,][nnn])

    abline(v=ci.ind[ii], col="gray", lty=2)
    abline(h=di[ii], col="blue", lty=3)
    segments(x0=ci.ind[ii], y0=sortdist[ii,][ci.ind[ii]], x1=ci.ind[ii], y1=sortdist[ii,][ci.ind[ii]+1], col=2)
    text(x = ci.ind[ii]+8, y = mean(sortdist[ii,][c(ci.ind[ii],(ci.ind[ii]+1))]),
         expression(c[i]), cex = 1.3, col=2)
    points(x = par("usr")[1], y=di[ii], pch = 16, cex=1.3, col = "blue", las = 1, xpd = TRUE)
    text(x = par("usr")[1], y=di[ii], #labels= expression(d[i]),
         labels= expression(paste(d[i], "*")),
         cex=1.1, col = "blue", las = 1, xpd = TRUE, pos=4)
    points(x=ci.ind[ii], y = par("usr")[3], pch = 16, cex=1.3, col = "gray", las = 1, xpd = TRUE)
    text(x=ci.ind[ii], y = par("usr")[3], labels= expression(paste(italic(k), "*", sep="")),
         cex=1.1, col = "gray", las = 1, xpd = TRUE, pos=3)
    #---------------------------------------------------

    if (!is.null(savefile)) dev.off()
    if (is.null(savefile))  readline(prompt="Press [enter] to continue")
  }

  return(first.ind.in.cluster)
}
