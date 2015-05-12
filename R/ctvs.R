#' ctvs
#'
#' collection of 'CurveToValue' functions
#'
#' @param data.object input data (gloDatamix format)
#' @param ctv.type type of CTV to calculate
#' @param rec.frames number of frames recorded, will be autodetected if empty
#' @param bg.firstframe for CTVs with background subtraction this is the position of the first bg frame
#' @param bg.lastframe 2nd frame for background (end of window)
#' @param sig.firstframe signal onset
#' @param sig.lastframe signal offset
#' @param wa calc "width @@" xx*extremum, needs 2 arguments, e.g. width when 0.2 of peakmax is reached to when peak falls to .9 of peakmax
#' @param plot perform a barplot in the end?
#' @param report write CTV report file (gives single traces with signal, bg ... marked and values in legend)
#' @param ylim ylim for report plots, given as 2 value vector, e.g. c(-1,5) plots from -1 to
#' @param mfrow numer of columns and rows to plot per page
#' @param suffix optionally, ad describing suffix to filenames, name of data.object is used when empty
#' @param meanTraces does data come from meanTraces.R?
#' @param fp numerical vector of length 4 giving positions for ctv.w for putative peak on and offset
#'
#' @author Daniel Münch <daniel@@muench.bio>

ctvs <- function(data.object, ctv.type, rec.frames=NA, bg.firstframe=NA, bg.lastframe=NA, sig.firstframe=NA, sig.lastframe=NA, wa = c(.2,.7), plot=F, report=F, ylim=NA, mfrow=c(4,3), suffix="", meanTraces=F, fp=c(29,37,41,49)) {

  #get data.object name
  data.object.name <- deparse(substitute(data.object))
  data.object.name <- gsub("[^a-zA-Z0-9 :._-]", "", data.object.name)
  if (suffix == "") suffix <- data.object.name    #use data.object.name if suffix is empty

  ### get number of frames recorded
  if (is.na(rec.frames)){
    rec.frames <- length(names(data.object)) - which(names(data.object) == "data0") + 1
    cat(paste("number of frames was set to",rec.frames,"\n"))
  }

  ### divide headers from data
  ctv.data <- subset(data.object, select=(data0:get(paste("data",rec.frames-1,sep=""))))
  ctv.headers <- subset(data.object, select=(-(data0:get(paste("data",rec.frames-1,sep="")))))


  ### function to find extremum
  extremum <- function(x) {
    min <- min(x)
    #	minpos <- which.min(x)
    min.abs <- abs(min)
    max <- max(x)
    #	maxpos <- which.max(x)
    ifelse (min.abs > max, return(min), return(max))
  }

  which.extremum <- function(x) {
    min <- min(x)
    minpos <- which.min(x)
    min.abs <- abs(min)
    max <- max(x)
    maxpos <- which.max(x)
    ifelse (min.abs > max, return(minpos), return(maxpos))
  }


    {
    if (is.na(ylim[1]) == T) ylim.lower <- min(ctv.data)
      else ylim.lower <- ylim[1]
    }
    {
    if (is.na(ylim[2]) == T) ylim.upper <- max(ctv.data)
      else ylim.upper <- ylim[2]
    }


######################### START OF CTV SECTION ##############

# CTV 1 -------------------------------------------------------------------
  # looks for max between sig.firstframe and sig.lastframe, no background substraction
  if (ctv.type == 1) {
    ctv <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, max)
	  # data for plots
	  sig.pos <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, which.max) + sig.firstframe -1	# position of signal
	  sig <- ctv
  }

# CTV 22 ------------------------------------------------------------------
  ### takes sig.firstframe and substracts mean of 3 frames around bg.firstframe
  if (ctv.type == 22) {
    sig <- apply(ctv.data[,(sig.firstframe-1):(sig.firstframe+1)],1,mean)
    bg <- apply(ctv.data[,(bg.firstframe-1):(bg.firstframe+1)],1,mean)
    ctv <- sig - bg

  	#for plots
	  sig.pos <- rep(sig.firstframe, length(sig))
	  bg.pos <- rep(bg.firstframe, length(bg))
  }

# CTV 23 ------------------------------------------------------------------


  ### takes mean of 3 frames around sig.firstframe
  if (ctv.type == 23) {
    ctv <- apply(ctv.data[,(sig.firstframe-1):(sig.firstframe+1)],1,mean)

    #for plots
    sig.pos <- rep(sig.firstframe, length(ctv))
  }

# CTV 25 ------------------------------------------------------------------
### Anjas CTV, "an average of 16 frames before first stimulus onset were subtracted from an average of 21 frames starting 1 1/2 sec after first stimulus onset"

  if (ctv.type == 25) {
    sig <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, mean)
    bg <- apply(ctv.data[bg.firstframe:bg.lastframe], 1, mean)
    ctv <- sig - bg
  }

# CTV 26 ------------------------------------------------------------------
  ### like 25 but taking the sum under the curve
  if (ctv.type == 26) {
    sig <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, sum)
    bg <- apply(ctv.data[bg.firstframe:bg.lastframe], 1, mean)
    ctv <- sig - bg
  }

# CTV 99 ------------------------------------------------------------------
    ### looks for extremum between sig.firstframe and sig.lastframe and substracts mean of background between bg.firstframe and bg.lastframe
    ### needs function "extremum"
  if (ctv.type == 99) {
    ctv.a <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, extremum)  # find extremum
    ctv.b <- apply(ctv.data[bg.firstframe:bg.lastframe], 1, mean)       # background to substract
    ctv <- ctv.a - ctv.b
    ### write data for plots, not needed for anything else
    bg <- ctv.b 	# substracted background
	  sig <- ctv.a	# calculated signal
	  sig.pos <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, which.extremum) + sig.firstframe -1	# position of signal
  }

# CTV 100 -----------------------------------------------------------------
    ### takes average of 3 frames around extremum between sig.firstframe and sig.lastframe and substracts mean of background between bg.firstframe and bg.lastframe
    ### needs function "extremum" (top of this function)
  if (ctv.type == 100) {
    ctv <- numeric()
    bg <- numeric()		# only needed for plots in report
    sig <- numeric()	# only needed for plots in report
    for (i in 1:dim(ctv.data)[1]) {
      signal.data <- ctv.data[i,sig.firstframe:sig.lastframe]
      ext.pos <- which.extremum(signal.data)+sig.firstframe-1
      ctv.a <- mean(as.numeric(ctv.data[i,(ext.pos-1):(ext.pos+1)]))
      ctv.b <- mean(as.numeric(ctv.data[i,bg.firstframe:bg.lastframe]))
      ctvx <- ctv.a - ctv.b
      ctv <- c(ctv,ctvx)
       ### write data for plots, not needed for anything else
      bg <- c(bg,ctv.b)	# background
      sig <- c(sig,ctv.a)	# signal
    }
    sig.pos <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, which.extremum) + sig.firstframe -1 	# position of signal
  }


# CTV 101 -----------------------------------------------------------------
  ### takes average of 3 frames around maximum between sig.firstframe and sig.lastframe and substracts mean of background between bg.firstframe and bg.lastframe
  if (ctv.type == 101) {
    ctv <- numeric()
    bg  <- numeric()		# only needed for plots in report
    sig <- numeric()	# only needed for plots in report
    for (i in 1:dim(ctv.data)[1]) {
      signal.data <- ctv.data[i,sig.firstframe:sig.lastframe]
      ext.pos     <- which.max(signal.data)+sig.firstframe-1
      ctv.a       <- mean(as.numeric(ctv.data[i,(ext.pos-1):(ext.pos+1)]))
      ctv.b       <- mean(as.numeric(ctv.data[i,bg.firstframe:bg.lastframe]))
      ctvx        <- ctv.a - ctv.b
      ctv         <- c(ctv,ctvx)
      ### write data for plots, not needed for anything else
      bg  <- c(bg,ctv.b)	# background
      sig <- c(sig,ctv.a)	# signal
    }
    sig.pos <- apply(ctv.data[sig.firstframe:sig.lastframe], 1, which.max) + sig.firstframe -1 	# position of signal
  }

# CTV 102 -----------------------------------------------------------------
### copy of CTV -102 in IDL
### takes average of 3 frames around maximum between fram 23 and 40
  if (ctv.type == 102) {
    ctv <- numeric()
    sig <- numeric()
    for (i in 1:dim(ctv.data)[1]) {
      signal.data <- ctv.data[i,23:40]
      max.pos <- which.max(signal.data)+23-1
      ctvx <- mean(as.numeric(ctv.data[i,(max.pos-1):(max.pos+1)]))
      ctv <- c(ctv,ctvx)
   		### write data for plots, not needed for anything else
      sig <- c(sig,ctvx)	# signal
    }
    sig.pos <- apply(ctv.data[23:40], 1, which.max) + 23 - 1 	# position of signal for plotting
  }

# CTV 102a ----------------------------------------------------------------
  ### LIKE 102 BUT ONLY TILL FRAME 38
  ### takes average of 3 frames around maximum between fram 23 and 38
  if (ctv.type == "102a") {
    ctv <- numeric()
    sig <- numeric()
    for (i in 1:dim(ctv.data)[1]) {
      signal.data <- ctv.data[i,23:38]
      max.pos <- which.max(signal.data)+23-1
      ctvx <- mean(as.numeric(ctv.data[i,(max.pos-1):(max.pos+1)]))
      ctv <- c(ctv,ctvx)
      ### write data for plots, not needed for anything else
      sig <- c(sig,ctvx)	# signal
    }
    sig.pos <- apply(ctv.data[23:38], 1, which.max) + 23 - 1 	# position of signal for plotting
  }

# CTV missanga ------------------------------------------------------------
  if (ctv.type == "missanga1") {
    ctv <- numeric()
    bg  <- numeric()
    sig <- numeric()
    for (i in 1:dim(ctv.data)[1]) {
      ctv.bg  <- mean(as.numeric(ctv.data[i,15:20]))
      ctv.sig <- mean(as.numeric(ctv.data[i,c(25:35,45:55,65:70,75:77,80:102)]))
      ctv <- c(ctv, ctv.sig - ctv.bg)
      bg  <- c(bg, ctv.bg)		# bg values for plots
      sig <- c(sig, ctv.sig)	# sig values for plots
    }
    # polygone for the report
    pol.bg.x  <- c(15, 15, 20, 20)
    pol.bg.y  <- c(ylim.lower, ylim.upper, ylim.upper, ylim.lower)
    pol.sig.x <- c(25,25,35,35,45,45,55,55,65,65,70,70,75,75,77,77,80,80,102,102)
    pol.sig.y <- c(ylim.lower, ylim.upper, ylim.upper, ylim.lower, ylim.lower, ylim.upper, ylim.upper, ylim.lower, ylim.lower, ylim.upper, ylim.upper, ylim.lower, ylim.lower, ylim.upper, ylim.upper, ylim.lower, ylim.lower, ylim.upper, ylim.upper, ylim.lower)
  }

# ctv.w -------------------------------------------------------------------
  ###
  # sig.pos == extr.pos
  # pol.sig.x == pol for extr. seach
  # width.x .y 00 coordinates for line
  # peaks.x ??as list with object names $x and $y?
  # pit.x

  if (ctv.type == "ctv.w") {
    require(pastecs)
    ctv <- data.frame()

    extr     <- apply(ctv.data[sig.firstframe:sig.lastframe],1,extremum)
    extr.pos <- apply(ctv.data[sig.firstframe:sig.lastframe],1,which.extremum) + sig.firstframe - 1

    max     <- apply(ctv.data[sig.firstframe:sig.lastframe],1,max)
    max.pos <- apply(ctv.data[sig.firstframe:sig.lastframe],1,which.max) + sig.firstframe - 1

    min     <- apply(ctv.data[sig.firstframe:sig.lastframe],1,min)
    min.pos <- apply(ctv.data[sig.firstframe:sig.lastframe],1,which.min) + sig.firstframe - 1

    fp1.extr     <- apply(ctv.data[fp[1]:fp[2]],1,extremum)
    fp1.extr.pos <- apply(ctv.data[fp[1]:fp[2]],1,which.extremum) + fp[1] - 1

    #take peaks from fixed windows
    max1   <- apply(ctv.data[fp[1]:fp[2]],1,max)
    max1.p <- apply(ctv.data[fp[1]:fp[2]],1,which.max) + fp[1] -1
    min1   <- apply(ctv.data[fp[1]:fp[2]],1,min)
    min1.p <- apply(ctv.data[fp[1]:fp[2]],1,which.min) + fp[1] -1

    max2   <- apply(ctv.data[(fp[2]+1):fp[3]],1,max)
    max2.p <- apply(ctv.data[(fp[2]+1):fp[3]],1,which.max) + fp[2]
    min2   <- apply(ctv.data[(fp[2]+1):fp[3]],1,min)
    min2.p <- apply(ctv.data[(fp[2]+1):fp[3]],1,which.min) + fp[2]

    max3   <- apply(ctv.data[(fp[3]+1):fp[4]],1,max)
    max3.p <- apply(ctv.data[(fp[3]+1):fp[4]],1,which.max) + fp[3]
    min3   <- apply(ctv.data[(fp[3]+1):fp[4]],1,min)
    min3.p <- apply(ctv.data[(fp[3]+1):fp[4]],1,which.min) + fp[3]

    biphas <- apply(ctv.data[51:60],1,mean)

    fp.pvdiff <- numeric()
    for(i in 1:dim(ctv.data)[1]) {
      if(extr[i]>=0) {
        fp.pvdiff <- c(fp.pvdiff, max1[i]-min2[i])
        } else {
          fp.pvdiff <- c(fp.pvdiff, min1[i]-max2[i])
          }
    }

    # biphasic <- extr < 0 & max1 > abs(.5*extr)
    # biphasic <- extr < 0 & max1 > abs(min1)
    biphasic <- biphas < -0.2 & max1 > abs(min1)
    # biphasic <- biphas

    cat("calculating ctv.w...\n")
    pb <- txtProgressBar(min = 0, max = dim(ctv.data)[1], style = 3)
    for (i in 1:dim(ctv.data)[1]) {
      # reset variables to prevent them appearing in next measurement if those can not be calculated
      ctv.temp <- NA; peak.1.x <- NA; peak.2.x <- NA; peak.1.y <- NA; peak.2.y <- NA; pit.x <- NA; pit.y <- NA; peak.2.rel <- NA; peakvalleydiff <- NA; x.width <- NA; p.low <- NA; p.high <- NA; tp <- NA; two.peaks <- NA; pits.between <- NA; one.peak <- NA

      x <- as.numeric(ctv.data[i,]); names(x) <- 1:rec.frames

      if (extr[i] >= 0 | biphasic[i] == T) {
        p.low  <- which(x[sig.firstframe:max.pos[i]]>= max[i]*wa[1])[1] + sig.firstframe - 1
        p.high <- rev(which(x[max.pos[i]:rec.frames] >= max[i]*wa[2]))[1] + max.pos[i] -1
        if (p.low == max.pos[i])  p.low  <- p.low - 1                              # make sure p.low is not equal to max.pos
        if (p.high == max.pos[i]) p.high <- p.high + 1                           # make sure p.high is not equal to max.pos
        #lines(x = c(p.low,p.high), y = rep(wa*max[i],2), col=2)
        x.width <- p.high - p.low
      } else {
        p.low  <- which(x[sig.firstframe:extr.pos[i]] <= extr[i]*wa[1])[1] + sig.firstframe - 1
        p.high <- rev(which(x[extr.pos[i]:rec.frames] <= extr[i]*wa[2]))[1] + extr.pos[i] -1
        if (p.low  == extr.pos[i]) p.low  <- p.low - 1                              # make sure p.low is not equal to extr.pos
        if (p.high == extr.pos[i]) p.high <- p.high + 1                           # make sure p.high is not equal to extr.pos
        #lines(x = c(p.low,p.high), y = rep(wa*extr[i],2), col=2)
        x.width <- p.low - p.high
        }

      tp.range <- fp[c(1,4)]
      tp       <- turnpoints(x[tp.range[1]:tp.range[2]])  # changed, should search for peaks only in area wherer peaks are common, not whole width #turnpoints(x[p.low:p.high])


      if (tp$nturns > 0) {
        if (tp$firstispeak == T & !is.na(tp$firstispeak)) {typep <- c("peak", "pit")} else {typep <- c("pit", "peak")}      # I defined it for NA but don't know why it became NA...
        typepts  <- rep(typep, length.out = tp$nturns)
        tablepts <- as.data.frame(list(point = tp$tppos, type = typepts,
                                       proba = tp$proba, info = tp$info))
        tablepts$point <- tablepts$point + tp.range[1] -1                                       # set back to original frame numbers

        peaks <- subset(tablepts, type=="peak")
        pits  <- subset(tablepts, type=="pit")
        if (extr[i] < 0 & biphasic[i] == F) {temp <- peaks; peaks <- pits; pits <- temp}                       # swap pits and peaks if extremum is negativ i.e. response is inhibitory

        if (dim(peaks)[1] > 1) {
          two.peaks <- peaks[order(peaks$proba),][1:2,]                                   # select the two most probable
          two.peaks <- two.peaks[order(two.peaks$point),]                                 # sort according to time again
          pits.between <- subset(pits, point > min(two.peaks$point) & point < max(two.peaks$point))
          pit.between <- pits.between[order(pits.between$proba),][1,]

          #points(x=two.peaks$point, y=x[two.peaks$point], col=5)
          #points(x=pit.between$point, y=x[pit.between$point], col=6)

          peak.1.x <- two.peaks$point[1]
          peak.2.x <- two.peaks$point[2]
          peak.1.y <- x[peak.1.x]
          peak.2.y <- x[peak.2.x]
          pit.x    <- pit.between$point
          pit.y    <- x[pit.x]

          peak.2.rel <- peak.2.y/peak.1.y
          peakvalleydiff <- peak.1.y-pit.y
        } else {
          one.peak <- peaks[order(peaks$proba),][1,]
          peak.1.x <- one.peak$point[1]
          peak.1.y <- x[peak.1.x]
          #points(x=one.peak$point, y=x[one.peak$point], col=5)
          }
      }



      # calculate a "peakyness"
      p.whole.area <- extr[i] * x.width
      p.curve.area <- sum(x[p.low:p.high])
      p.diff       <- p.whole.area - p.curve.area
      peakyness    <- p.diff / p.curve.area       # area above curve / area under curve
      #peakyness <- paste(round(p.whole.area,2),"-",round(p.curve.area,2)," // ",round(p.diff,2),"/",round(p.curve.area,2), "==",round(peakyness,2))

      ctv.temp <- data.frame("ctv.w_width"     = x.width,
                             "ctv.w_peak.low"  = p.low,
                             "ctv.w_peak.high" = p.high,
                             "ctv.w_extr.x"    = extr.pos[i],
                             "ctv.w_extr.y"    = extr[i],
                             "ctv.w_peaks"     = sum(tp$peaks),
                             "ctv.w_pits"      = sum(tp$pits),
                             "ctv.w_p1.x"      = peak.1.x,
                             "ctv.w_p1.y"      = peak.1.y,
                             "ctv.w_p2.x"      = peak.2.x,
                             "ctv.w_p2.y"      = peak.2.y,
                             "ctv.w_pit.x"     = pit.x,
                             "ctv.w_pit.y"     = pit.y,
                             "ctv.w_pvdiff"    = peakvalleydiff,
                             "ctv.w_p2.rel"    = peak.2.rel,
                             "ctv.w_peakyness" = peakyness,
                             "fp1.max.y"       = max1[i],
                             "fp1.max.x"       = max1.p[i],
                             "fp1.min.y"       = min1[i],
                             "fp1.min.x"       = min1.p[i],

                             "fp2.max.y"       = max2[i],
                             "fp2.max.x"       = max2.p[i],
                             "fp2.min.y"       = min2[i],
                             "fp2.min.x"       = min2.p[i],

                             "fp3.max.y"       = max3[i],
                             "fp3.max.x"       = max3.p[i],
                             "fp3.min.y"       = min3[i],
                             "fp3.min.x"       = min3.p[i],

                             "fp.pvdiff"       = fp.pvdiff[i],

                             "biphasic"        = biphasic[i]
                             )

    ctv <- rbind(ctv, ctv.temp)

    setTxtProgressBar(pb, i)
    }
    close(pb)
    peak.width  <- ctv$ctv.w_width
    peak.low    <- ctv$ctv.w_peak.low
    peak.high   <- ctv$ctv.w_peak.high
    peak1.pos   <- ctv$ctv.w_p1.x
    peak1.value <- ctv$ctv.w_p1.y
    peak2.pos   <- ctv$ctv.w_p2.x
    peak2.value <- ctv$ctv.w_p2.y
    pit.pos     <- ctv$ctv.w_pit.x
    pit.value   <- ctv$ctv.w_pit.y
    extr.pos    <- ctv$ctv.w_extr.x
    extr.value  <- ctv$ctv.w_extr.y
    peaky.index <- ctv$ctv.w_peakyness
  }

# ctv.peaks ---------------------------------------------------------------
# get some peak features at fixed timepoints (mean across 5 frames)
  if(ctv.type == "ctv.peaks") {

  peak1  <- apply(ctv.data[,(fp[1]-2):(fp[1]+2)],1,mean)
  valley <- apply(ctv.data[,(fp[2]-2):(fp[2]+2)],1,mean)
  peak2  <- apply(ctv.data[,(fp[3]-2):(fp[3]+2)],1,mean)
  post   <- apply(ctv.data[,(fp[4]-2):(fp[4]+2)],1,mean)

  peak1.inh  <- apply(ctv.data[,(fp[1]):(fp[1]+4)],1,mean)
  peak2.inh  <- apply(ctv.data[,(fp[3]):(fp[3]+4)],1,mean)

  extr     <- apply(ctv.data[27:79],1,extremum)
  extr.pos <- apply(ctv.data[27:79],1,which.extremum) + 26

  peakyness = 1 - (valley / peak1)

  allwidth <- data.frame()
  for(i in 1:dim(ctv.data)[1]) {
    if (peak1[i] <= 0) width <- range(which(ctv.data[i,27:80] <= peak1[i]*.3) +26)
    if (peak1[i] >  0) width <- range(which(ctv.data[i,27:80] >= peak1[i]*.3) +26)
    allwidth <- rbind(allwidth,width)
  }
  names(allwidth) <- c("from","to")
  allwidth$width  <- allwidth$to - allwidth$from

  ctv <- data.frame(peak1  = peak1,
                    valley = valley,
                    peak2  = peak2,
                    post   = post,
                    peak1.inh = peak1.inh,
                    peak2.inh = peak2.inh,
                    extr = extr,
                    extr.pos = extr.pos,
                    peakyness = peakyness
                    )

  ctv <- cbind(ctv,allwidth)
  }



# CTV "ltk", lifetime kurtosis --------------------------------------------
if (ctv.type == "ltk") {
  ctv <- apply(ctv.data[sig.firstframe:sig.lastframe],1,ltk)
  }


###################### END OF CTV SECTION ###################

  # combine with headers again and return result
  ctv <- data.frame(ctv.headers, ctv)

  # if data comes from meanTraces.R, rename columns for plotting, then set back
  if(meanTraces == T) {
    colnames(ctv)[which(names(ctv) == "odor")] <- "TOdour"
    names(ctv)[which(names(ctv) == "concentration")] <- "NOConc"
    names(ctv)[which(names(ctv) == "glomerulus")] <- "NGloTag"
    }

  ### not ready yet
  ### if plot is set to true -> plot
  ### not working on mean traces as plots are performed animal-wise
  if (plot==T) {
    pdf(file=paste("ctv",ctv.type,suffix,"barplot.pdf",sep="_"), paper="a4", height=11.69 , width=8.27)
    par(mfrow=mfrow)

    for (i in 1:length(levels(factor(ctv$Tanimal)))) {
      animalx.name <- levels(factor(ctv$Tanimal))[i]
      animalx <- subset(ctv, Tanimal==animalx.name, select=c(TOdour,ctv))
      barplot(height=animalx$ctv, names.arg=animalx$TOdour, las=2, main=levels(ctv$Tanimal)[i])
    }
    dev.off()
  }

  if (report == T) {
    # for (i in 1:dim(data)[1]) {
    #  tracex <- subset(data.object[i,],select=c(data0:get(paste("data",rec.frames-1,sep=""))))
    #  ylim <- c(min(tracex),max(tracex))
    #  main <- paste(data.object[i,]$Tanimal," ",data.object[i,]$TOdour,sep="")
    #  plot(1:rec.frames, ylim=ylim, xlab="frames", ylab="response (deltaF/F)", main=main)
    # }

    #resort temporarely in order to plot ordered by ltk
    if (ctv.type == "ltk") {
      ctv.backup  <- ctv
      ctv.data    <- ctv.data[order(ctv$ctv),]
      ctv.headers <- ctv.headers[order(ctv$ctv),]
      ctv         <- ctv[order(ctv$ctv),] #has to be last, otherwise others are sorted in the wrong way
    }


  if (exists("bg")==F) bg <- NA
  if (exists("sig")==F) sig <- NA

  cat(paste("plotting",dim(data.object)[1],"measurements into "));cat(paste("ctv",ctv.type,suffix,"report.pdf...\n",sep="_"))

  pdf(file=paste("ctv",ctv.type,suffix,"report.pdf",sep="_"), paper="a4", height=11.69 , width=8.27)
  par(mfrow=mfrow)
    plot(1:rec.frames,ctv.data[1,], main="parameter settings",type="n",ylim=c(ylim.lower,ylim.upper), xlab="", ylab="")
    legend("center",legend=c(paste("ctv = ",ctv.type,sep=""),paste("bg.firstframe = ",bg.firstframe,sep=""),paste("bg.lastframe = ",bg.lastframe,sep=""),paste("sig.firstframe = ",sig.firstframe,sep=""),paste("sig.lastframe = ",sig.lastframe,sep=""),paste("rec.frames = ",rec.frames,sep="")), bty="n",cex=1,pch=-1)

    pb <- txtProgressBar(min = 0, max = dim(data.object)[1], style = 3)
    for (i in 1:dim(data.object)[1])
      {
      #ex.pos <- which(ctv.data[i,] == ctv$ctv[i])[1]
      plot(1:rec.frames,ctv.data[i,],
           main=paste(
             if(!is.null(ctv$Tanimal[i]))   paste("a: ",     ctv$Tanimal[i], "\n", sep=""),
             if(!is.null(ctv$receptor[i]))  paste("r: ",     ctv$receptor[i], "\n", sep=""),
             if(!is.null(ctv$TOdour[i]))    paste("o: ",     ctv$TOdour[i],  " ", sep=""),
             if(!is.null(ctv$NOConc[i]))    paste("c: ",     ctv$NOConc[i],  " ", sep=""),
             if(!is.null(ctv$NGloTag[i]) & length(levels(factor(ctv$NGloTag))) > 1) paste("g: ",     ctv$NGloTag[i], " ", sep=""),
             if(!is.null(ctv$TPharma[i]))   paste("p: ",     ctv$TPharma[i], " ", sep=""),
             if(!is.null(ctv$n[i]))         paste("n: ",     ctv$n[i]),
             sep=""),
             type="n",xlab="frames",ylab=expression(Delta*F/F),ylim=c(ylim.lower,ylim.upper))

      if (exists("pol.bg.x")==T)
        {
        polygon(y=pol.bg.y, x=pol.bg.x, density=-1, col="#F5FF99", border=NA)
        }
        else  polygon(y=c(ylim.lower,ylim.lower,ylim.upper,ylim.upper), x=c(bg.firstframe,bg.lastframe,bg.lastframe,bg.firstframe), density=-1, col="#F5FF99", border=NA)

      if (exists("pol.sig.x")==T)
        {
        polygon(y=pol.sig.y, x=pol.sig.x, density=-1, col="#C2FF99", border=NA)
        }
        else  polygon(y=c(ylim.lower,ylim.lower,ylim.upper,ylim.upper), x=c(sig.firstframe,sig.lastframe,sig.lastframe,sig.firstframe), density=-1, col="#C2FF99", border=NA)


  	  lines(1:rec.frames,ctv.data[i,], main=paste(ctv$Tanimal[i],ctv$TOdour[i]),type="l")

      if (exists("sig.pos")) points(sig.pos[i],sig[i],col="red")
      if (exists("bg.pos")) points(bg.pos[i],bg[i],col="green")

      if (exists("peak1.pos")) points(peak1.pos[i], peak1.value[i], col="blue")
      if (exists("peak1.value")) points(peak2.pos[i], peak2.value[i], col="blue")
      if (exists("pit.pos")) points(pit.pos[i], pit.value[i], col="magenta")
      if (exists("peak.width")) abline(v = c(peak.low[i], peak.high[i]), lty=2, col=2)

      if (ctv.type == "ctv.w")
        {
        legend("topleft",legend=c(paste("p1 =",peak1.pos[i]),paste("p2 =",peak2.pos[i]),paste("width =",peak.width[i]),paste("peak-valley diff. =",round(ctv$ctv.w_pvdiff[i],3)),paste("peak2 / peak1 =",round(ctv$ctv.w_p2.rel[i],3)),paste("peakiness =",round(peaky.index[i],2)), paste("biphasic:", biphasic[i])), bty="n",cex=.6,pch=-1)
        } else {
          legend("topleft",legend=c(paste("bg = ",round(bg[i],3),sep=""),paste("sig = ",round(sig[i],3),sep=""),paste("ctv = ",round(ctv$ctv[i],3),sep="")), bty="n",cex=.9,pch=-1)

        }
    setTxtProgressBar(pb,i)
    }
    close(pb)
  dev.off()
  cat("done!\n")

  }

  #restore old order if coming from ltk
  if (ctv.type == "ltk") ctv <- ctv.backup

  ### if data comes from meanTraces.R, set back column names
  if(meanTraces == T) {
    names(ctv)[which(names(ctv) == "TOdour")]  <- "odor"
    names(ctv)[which(names(ctv) == "NOConc")]  <- "concentration"
    names(ctv)[which(names(ctv) == "NGloTag")] <- "glomerulus"
  }

return(ctv)
}
