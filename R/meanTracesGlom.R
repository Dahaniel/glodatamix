#' meanTraces
#'
#' calculate and plot mean traces from gloDatamix data
#'
#' @param data.object gloDatamix data.frame
#' @param rec.frames number of frames measured
#' @param st.shift number of frames to shift the stimulus bars
#' @param ylim.lower limits of the y axis, if nothing is specified the max/min of the dataset will be used; it is also okay to only specify one limit, then only the other one will be calculated
#' @param ylim.upper limits of the y axis, if nothing is specified the max/min of the dataset will be used; it is also okay to only specify one limit, then only the other one will be calculated
#' @param add also plot single recordings ("single"), standard deviation ("sd") or standard error of mean ("sem"); default is "single"
#' @param method use "median" or "mean", default is "mean"
#' @param pdf plot into a PDF?
#' @param mfrow nr of rows/columns to plot
#' @param receptor adds a column containing the receptor name provided here
#' @param return return the results in the end or jus produce a pdf?
#'
#' @author Daniel Münch <daniel@@muench.bio>
#'
#' @return list and PDF

meanTraces <- function( data.object, rec.frames = NA, st.shift = NA, ylim.lower = NA, ylim.upper = NA, add = "single", method = "mean", pdf = T,  mfrow = c(4,3), receptor = NA, return = TRUE) {

# if rec.frame is not defined get it from first row of dataset
  if (is.na(rec.frames)) {
   rec.frames <- data.object$NNoFrames[1]
   print(paste("number of frames was set to",rec.frames))
  }

# get stimulus times:
  if (!is.null(data.object$NStim_on) & is.numeric(data.object$NStim_on)) {
    st1.on <- data.object$NStim_on[1] + 1
    st1.off <- data.object$NStim_off[1] + 1
    cat(paste("stimulus 1 was set to frame",st1.on,"and",st1.off,"\n"))
  }

  if (!is.null(data.object$Nstim2ON & is.numeric(data.object$Nstim2ON))) {
    st2.on <- data.object$Nstim2ON[1] + 1
    st2.off <- data.object$Nstim2OFF[1] + 1
    cat(paste("stimulus 2 was set to frame",st2.on,"and",st2.off,"\n"))
  }

### get odor names
  odors <- sort(levels(data.object$TOdour))

### get colors per animal ###
  animalcolors <- data.frame("animal"=levels(data.object$Tanimal),"color"=rainbow(length(levels(data.object$Tanimal))))

### define position of datapoints     # maybe change to use column names instead of
  firstframe <- length(data.object) - rec.frames + 1
  lastframe  <- length(data.object)

### define stimulus times
  if (!is.na(st.shift) & exists("st1.on")) {
    st1.on <- st1.on + st.shift
    st1.off <- st1.off + st.shift
  }

  if (!is.na(st.shift) & exists("st2.on")) {
    st2.on <- st2.on + st.shift
    st2.off <- st2.off + st.shift
  }

### define ylims
  if (is.na(ylim.lower)) ylim.lower <-  min(data.object[,firstframe:lastframe])
  if (is.na(ylim.upper)) ylim.upper <-  max(data.object[,firstframe:lastframe])


### get glomeruli tested
  gloms <- sort(levels(as.factor(data.object$NGloTag)), decreasing=F)

### get data.object.name for filenames
  data.object.name <- deparse(substitute(data.object))

### select glomwise
  calculated.median <- data.frame()
  calculated.mean <- data.frame()
  calculated.sd <- data.frame()
  calculated.sem <- data.frame()
  calculated.quantile <- data.frame()

  for (n in 1:length(gloms)) {
    glomx       <- subset(data.object, NGloTag==gloms[n])
    glomx.odors <- sort(levels(as.factor(as.character(glomx$TOdour))))
      if (pdf==T) {
      ### plot single and median/medians
      if (method=="median") pdf(file=paste("medianTraces_",data.object.name,"_",add,"_",gloms[n],".pdf",sep=""), paper="a4", height=11.69 , width=8.27)
      if (method=="mean") pdf(file=paste("meanTraces_",data.object.name,"_",add,"_",gloms[n],".pdf",sep=""), paper="a4", height=11.69 , width=8.27)
      par(mfrow=mfrow)
    }

    pb <- txtProgressBar(min = 0, max = length(glomx.odors), style = 3)
    for (i in 1:length(glomx.odors)) {
      glomx.odorx <- subset(glomx, TOdour==glomx.odors[i])
      glomx.odorx.concentrations <- sort(levels(as.factor(as.character(glomx.odorx$NOConc))))
      for (j in 1:length(glomx.odorx.concentrations)) {
        odorx          <- subset(glomx.odorx, NOConc==glomx.odorx.concentrations[j])
        odorx.median   <- apply(odorx[,firstframe:lastframe],2,median)
        odorx.mean     <- apply(odorx[,firstframe:lastframe],2,mean)
        odorx.sd       <- apply(odorx[,firstframe:lastframe],2,sd)
        odorx.quantile <- apply(odorx[,firstframe:lastframe],2,quantile)
        odorx.n        <- dim(odorx)[1]
        odorx.sem     <- odorx.sd/sqrt(odorx.n)

        if (pdf==T) {
          plot(1:(rec.frames), odorx[1,firstframe:lastframe], type="n", ylim=c(ylim.lower,ylim.upper), main=paste(glomx.odors[i], "(",glomx.odorx.concentrations[j],")"," g",gloms[n], "( n=", odorx.n, ")",sep=""), xlab="frames", ylab=expression(Delta*F/F))
          if(exists("st1.on")) polygon(y=c(ylim.lower,ylim.lower,ylim.upper,ylim.upper), x=c(st1.on,st1.off,st1.off,st1.on), density=-1, col="#efefef", border=NA)
          if(exists("st2.on")) polygon(y=c(ylim.lower,ylim.lower,ylim.upper,ylim.upper), x=c(st2.on,st2.off,st2.off,st2.on), density=-1, col="#efefef", border=NA)

          ### calculating colors for single measurement lines using: "grey(((2/3*x)+1):((2/3*x)+x) / (2*x))" for not getting black or white lines
          # x <- dim(odorx)[1]
          # color <- grey(((2/3*x)+1):((2/3*x)+x) / (2*x))
          color <- grey(5:(dim(odorx)[1]+4) / (dim(odorx)[1]+6))

          if (add == "single") for (k in 1:dim(odorx)[1]) lines(1:rec.frames, odorx[k,firstframe:lastframe], type="l", lty="11",col=color[k])
          if (add == "single_colourperanimal") {
            plotanimals <- match(odorx$Tanimal, animalcolors$animal)
            plotanimals <- animalcolors[plotanimals,]
            for (k in 1:dim(odorx)[1]) lines(1:rec.frames, odorx[k,firstframe:lastframe], type="l", lty="11",col=as.character(plotanimals[k,"color"]))
            legend(x="topright",as.character(plotanimals$animal), text.col=as.character(plotanimals$color), cex=0.4, bty="n")            }

          if (add == "sd" & method=="median")  polygon(c(1:rec.frames,rec.frames:1),c(odorx.median-odorx.sd,rev(odorx.median+odorx.sd)),col="#ffccff",border=NA)
          if (add == "sem" & method=="median")  polygon(c(1:rec.frames,rec.frames:1),c(odorx.median-odorx.sem,rev(odorx.median+odorx.sem)),col="#ffccff",border=NA)
          if (add == "sd" & method=="mean")  polygon(c(1:rec.frames,rec.frames:1),c(odorx.mean-odorx.sd,rev(odorx.mean+odorx.sd)),col="#ffccff",border=NA)
          if (add == "sem" & method=="mean")  polygon(c(1:rec.frames,rec.frames:1),c(odorx.mean-odorx.sem,rev(odorx.mean+odorx.sem)),col="#ffccff",border=NA)

          if (add == "quant2575")  polygon(c(1:rec.frames,rec.frames:1),c(odorx.quantile[2,],rev(odorx.quantile[4,])),col="#ffccff",border=NA)

          if (method == "median") lines(1:rec.frames, odorx.median, type="l",col="#ff00ff99", lwd=1.8)
          if (method == "mean") lines(1:rec.frames, odorx.mean, type="l",col="#ff00ff99", lwd=1.8)
        }

        calculated.odorx.median <- data.frame("receptor"=receptor,"odor"=glomx.odors[i], "concentration"=glomx.odorx.concentrations[j],"glomerulus"=gloms[n],"n"=odorx.n, t(odorx.median))
        calculated.odorx.mean <- data.frame("receptor"=receptor,"odor"=glomx.odors[i], "concentration"=glomx.odorx.concentrations[j],"glomerulus"=gloms[n],"n"=odorx.n,t(odorx.mean))
        calculated.odorx.sd <- data.frame("receptor"=receptor,"odor"=glomx.odors[i], "concentration"=glomx.odorx.concentrations[j],"glomerulus"=gloms[n],"n"=odorx.n,t(odorx.sd))
        calculated.odorx.sem <- data.frame("receptor"=receptor,"odor"=glomx.odors[i], "concentration"=glomx.odorx.concentrations[j],"glomerulus"=gloms[n],"n"=odorx.n,t(odorx.sem))
        calculated.odorx.quantile <- data.frame("receptor"=receptor,"odor"=glomx.odors[i], "concentration"=glomx.odorx.concentrations[j],"glomerulus"=gloms[n],"n"=odorx.n,"quantile"=rownames(odorx.quantile),odorx.quantile)

        calculated.median <- rbind(calculated.median,calculated.odorx.median)
        calculated.mean <- rbind(calculated.mean,calculated.odorx.mean)
        calculated.sd <- rbind(calculated.sd,calculated.odorx.sd)
        calculated.sem. <- rbind(calculated.sem,calculated.odorx.sem)
        calculated.quantile <- rbind(calculated.quantile,calculated.odorx.quantile)

       }
    setTxtProgressBar(pb, i)
    }
    close(pb)
    if (pdf==T) dev.off()
    print(paste("glomerulus ",gloms[n], " [",n,"/",length(gloms),"] ..... done .....", sep=""))
  }



  if (return == TRUE) {
    list(median = calculated.median,
    mean        = calculated.mean,
    sd          = calculated.sd,
    sem         = calculated.sem,
    quantile    = calculated.quantile)
  }
}
