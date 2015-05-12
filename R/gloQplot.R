#' gloQplot
#'
#' quick and dirty plots from gloDatamix data
#'
#' @param data.object gloDatamix input data.frame
#' @param ylim plotting limit
#' @param mean plot mean trace?
#' @param median plot median trace?
#' @param legend plot a legend?
#' @param ... pass additional parameters to the plot() function
#'
#' @author Daniel Münch <daniel@@muench.bio>
#'
#' @return plot

gloQplot <- function(data.object, ylim = NULL, mean = F, median = F, legend = F, ...) {
  #f1 <- which(names(data.object) == "data0")
  #fl <- length(data.object)
  #frames <- length(names(data.object)) - f1  + 1
  f1     <- min(grep("data",names(data.object)))
  fl     <- max(grep("data",names(data.object)))
  frames <- length(grep("data",names(data.object)))
  n      <- dim(data.object)[1]
  col    <- rainbow(n)

  if(is.null(ylim)) ylim <- c(min(data.object[,f1:fl]),max(data.object[,f1:fl]))
  plot(as.numeric(data.object[1,f1:fl]), type="l", col=col[1], ylim=ylim, ...)
  for(i in 2:n) lines(as.numeric(data.object[i,f1:fl]), col=col[i])
  if (mean == T)   lines(as.numeric(apply(data.object[,f1:fl],2,mean)),col="#00000088", lty=2, lwd=5)
  if (median == T) lines(as.numeric(apply(data.object[,f1:fl],2,median)),col="#00000088", lty=3, lwd=5)
  if (legend == T) legend("topleft",legend=paste(data.object$TOdour, data.object$NOConc), fill=col)
}
