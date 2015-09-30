#' Plots the side-by-side boxplot of the within-cluster and between-cluster pairwise distances.
#'
#' @param GPout The output object created by the \code{GapProcedure} function.
#' @param classvec A vector indicating the class memberships of the data used to produce \code{GPout}.  The default (NULL) takes the classification vector stored in \code{GPout}.
#' @param Gindex A vector of cluster indices. If set to NULL (default), then the entire set of class labels is considered.
#' @param savefile A character vector.  If set to \code{NULL} (default), the figures are not saved and ploted in R, otherwise the figure will be saved to \code{<savefile>Cluster<g>.pdf}, where \code{g} denotes the cluster number.
#' @param ... options for \code{pdf()}
#' @details For each class label (\eqn{g = 1, \dots, G}), a side-by-side boxplot is created.
#' @seealso \code{\link[grDevices]{pdf}}
boxplotWB <- function(GPout, classvec=NULL, Gindex=NULL,  hcutoff = 0.05,show.singletons=FALSE, savefile=NULL, ...){


  if (is.null(classvec))
    classvec <- GPout$classification
  if (!show.singletons)
    classvec <- rmsingletons(classvec)
  distmat <- GPout$dist
  N <- ncol(distmat)
  if (is.null(Gindex)){
    Gindex <- seq(lu(classvec))
  } else {
    if (!all(Gindex %in% unique(classvec)))
      stop(paste(c("Problem with the \'Gindex\': \'classvec\' does not contain the label(s)",setdiff(Gindex,unique(classvec))), collapse = ", "))
  }

  mat <- array(1:length(distmat), dim=dim(distmat))
  rupper.tri <- Vectorize(function(x,p) {
    if (x==p)
      return(NULL)
    this <- (x+1):p
    return(cbind(x,this))
    }, vectorize.args="x")


  SSb = SSw = Nw = Nb = NULL
  for (g in Gindex){ #  g=1


    #---------------------------------------------------
    if (!is.null(savefile))
      pdf(paste(savefile,"Cluster", g,".pdf",sep=""),...)

    obs.in.gth.cluster <- which(classvec==g)
    obs.not.in.gth.cluster <- setdiff(1:N, obs.in.gth.cluster)

    upper.tri.g <- do.call(rbind, rupper.tri(obs.in.gth.cluster,N))
    upper.tri.ind <- mat[upper.tri.g]

    within.pairs <- combn(obs.in.gth.cluster, 2)
    wi.ind <- matrix2vecInd(within.pairs[1,],within.pairs[2,], dim(distmat))
    bw.ind <- upper.tri.ind[!(upper.tri.ind%in%wi.ind)]
    Nw[g] <- length(wi.ind)
    Nb[g] <- length(bw.ind)

    wi <- distmat[wi.ind]
    bt <- distmat[bw.ind]
    SSw[g] <- sum(wi)
    SSb[g] <- sum(bt)

    boxplot(wi, bt, names=c(expression(S[w]),expression(S[b])),main=paste("Cluster", g))
    if (!is.null(savefile)) dev.off()
    #---------------------------------------------------
    if (is.null(savefile)) readline(prompt="Press [enter] to continue")
  }

  # this should be true : ((sum(Nw) + sum(Nb)) == sum(upper.tri(distmat))) #[1] TRUE

  return(list(SSw=SSw, Nw=Nw, SSb=SSb, Nb=Nb))
}

# this <- boxplotWB(GP)
# library(clusterCrit)
# allind <- intCriteria(GP$dist, part=as.integer(GP$classification), crit="all")
# intCriteria(GP$dist, part=as.integer(GP$classification), crit="McClain_Rao")
#
# (sum(this$SSw)/sum(this$Nw))/(sum(this$SSb)/sum(this$Nb))
# sum(this$Nb)*sum(this$SSw)/(sum(this$Nw)*sum(this$SSb))
#
#
# all.equal((sum(this$SSw) + sum(this$SSb)), sum(GP$dist[upper.tri(GP$dist)]))
# ((sum(this$Nw) + sum(this$Nb)) == sum(upper.tri(distmat))) #[1] TRUE
#
#
#
#
