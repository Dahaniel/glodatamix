#' bleachCorrect
#'
#' correct bleaching in gloDatamix data
#'
#' @param data.object object containing the data
#' @param rec.frames number of frames recorded
#' @param bleach.begin number of the startframe for the fitting (defaults to 1)
#' @param bleach.leaveout numerical vector of frames to leave out for fitting (where the stimulus comes) (e.g.: c(27:65))
#' @param weight.BS how to weight the part before the stimulus
#' @param weight.AS how to weight the part after the stimulus
#' @param polcol color of the weight polygon
#' @param mfrow
#' @param fixedA fixed value for tparameter A
#' @param fixedB fixed value for tparameter B
#' @param fixedC fixed value for tparameter C
#' @param ref control to calculate the bleach correction on
#' @param estimate.B estimate parameter B from fits to the control?
#' @param bleach.begin.ref number of the startframe for the fitting (specific for fits to control)
#' @param bleach.leaveout.ref umerical vector of frames to leave out for fitting (where the stimulus comes) (e.g.: c(27:65))
#' @param weight.BS.ref how to weight the part before the stimulus (specific for fits to control)
#' @param weight.AS.ref how to weight the part after the stimulus (specific for fits to control)
#' @param diffParms4Ref define different fitting parameters for the reference (MOL)
#' @param plot plot?
#' @param ylim y limits for plots
#'
#' @author Daniel Münch <daniel@@muench.bio>

bleachCorrect <- function(data.object, rec.frames = NA, bleach.begin = 1, bleach.leaveout = NA, weight.BS = 8, weight.AS = 1, polcol = "#cccccc", mfrow = c(5,4), fixedA = NA, fixedB = NA, fixedC = NA, ref = "MOL", estimate.B = F, bleach.begin.ref =    NA, bleach.leaveout.ref = NA, weight.BS.ref = NA, weight.AS.ref = NA, diffParms4Ref = F, plot = T, ylim = NULL) {

  # get data.object name
  data.object.name <- deparse(substitute(data.object))
  data.object.name <- gsub("[^a-zA-Z0-9 :._-]", "", data.object.name)
  ifelse(estimate.B == T, assign("folder.estimate","_B.estimated"), assign("folder.estimate",""))
  folder <- paste("bleachCorrection.reports_",data.object.name,folder.estimate,sep="")


  if (is.na(bleach.begin.ref))
    bleach.begin.ref <- bleach.begin
  if (is.na(bleach.leaveout.ref)[1])
    bleach.leaveout.ref <- bleach.leaveout
  if (is.na(weight.BS.ref))
    weight.BS.ref <- weight.BS
  if (is.na(weight.AS.ref))
    weight.AS.ref <- weight.AS

  if(!is.null(ylim))
    ylims <- ylim

  if (is.na(fixedA) == F)
    A = fixedA
  if (is.na(fixedB) == F)
    B = fixedB
  if (is.na(fixedC) == F)
    C = fixedC

  # if rec.frame is not defined get it from first row of dataset
  if (is.na(rec.frames)) {
    rec.frames <- data.object$NNoFrames[1]
    print(paste("number of frames was automatically set to",rec.frames))
  }

  # get animals
  animals <- levels(data.object$Tanimal)

  # define weights
  weights <- rep(0, rec.frames)
  weights[bleach.begin:rec.frames] <- 1
  if (is.na(bleach.leaveout[1]) == F) {
    weights[bleach.leaveout] <- 0
    weights[bleach.begin:(bleach.leaveout[1] - 1)] <- weight.BS
    weights[(bleach.leaveout[length(bleach.leaveout)] + 1):rec.frames] <-
      weight.AS
  }

  if (diffParms4Ref == F) {
    weights.ref <- weights
  } else {
    weights.ref <- rep(0, rec.frames)
    weights.ref[bleach.begin.ref:rec.frames] <- 1
    if (is.na(bleach.leaveout.ref[1]) == F) {
      weights.ref[bleach.leaveout.ref] <- 0
      weights.ref[bleach.begin.ref:(bleach.leaveout.ref[1] - 1)] <-
        weight.BS.ref
      weights.ref[(bleach.leaveout.ref[length(bleach.leaveout.ref)] + 1):rec.frames] <-
        weight.AS.ref
    }
  }

  x <- 1:rec.frames
  bleachcorrected.data <- data.frame()
  all.corrected.measurements <- data.frame()
  not.fitted <- data.frame()

  if (plot == T)
    dir.create(path = folder)
  for (j in 1:length(animals)) {
    if (plot == T)
      pdf(file = paste(folder,"/",animals[j],"_bleachCorrectionReport.pdf",sep = ""), title = paste(animals[j],"_bleachCorrectionReport"), paper = "a4", height = 11.69 , width = 8.27)
    if (plot == T)
      par(mfrow = mfrow)

    animalx <- droplevels(subset(data.object, Tanimal == animals[j]))
    animalx.glomeruli <- levels(factor(animalx$NGloTag))
    all.corrected.measurements.animalx <- data.frame()
    for (g in 1:length(animalx.glomeruli)) {
      animalx.gx <-
        droplevels(subset(animalx, NGloTag == animalx.glomeruli[g]))
      if (estimate.B == T){
        animalx.gx.ref <- droplevels(subset(animalx.gx, TOdour == ref))
        ref.median <-apply(animalx.gx.ref[,(length(animalx.gx.ref) - rec.frames + 1):length(animalx.gx.ref)],2,median)
        min.ref <- min(ref.median)
        max.ref <- max(ref.median)
        ref.median <- as.numeric(ref.median)
        ref.nls <- nls(ref.median ~ (A * exp(-x / B) + C), trace = F, start = list(A = (max.ref - min.ref), B = 10, C = min.ref), weights = weights.ref)
        B.ref <- coef(ref.nls)["B"]
        A.ref <- coef(ref.nls)["A"]
        C.ref <- coef(ref.nls)["C"]

        fixedB <- B.ref
        B <- B.ref

        if (is.null(ylim))
          ylims <- c(min.ref,max.ref)

        if (plot == T) {
          ### plot
          plot(x, ref.median, type="n", col="darkgrey", lwd="1", ylim=ylims, xlab="frames", ylab="response (deltaF/F)", main=ref, cex.main=.7)
          # plot weights
          relative.weights.ref <- weights.ref/weight.BS.ref
          polygon.range.ref <- (max.ref-min.ref)/2
          polygon.weights.ref <- (relative.weights.ref * polygon.range.ref) +min.ref
          polygon.weights.ref[c(1,length(polygon.weights.ref))] <- min(polygon.weights.ref) # set first and last value to lowest level for plotting

          polygon(x, polygon.weights.ref, col=polcol, border=NA, density=-1)

          legend("topleft",legend=c(paste(round(A.ref,2),"*e^(-x/",round(B.ref,2),")+(",round(C.ref,2),")", sep=""), paste("glomerulus:",animalx.glomeruli[g])), bty="n",cex=.6,pch=-1)

          lines(x, ref.median, type="l", col="darkgrey", lwd="1")
          curve((A.ref * exp(-x/B.ref)+C.ref), from=1, to=rec.frames, add=T, col="#ee333399", lwd="3")
        }


        print(paste("parameter B (", round(B.ref,2),") calculated from medianTrace of",ref,"for",animals[j]))
      }



      all.corrected.measurements.animalx.gx <- data.frame()
      print(paste("calculating animal",animals[j]))
      ifelse(dim(animalx.gx)[1] > 1, PB <- T, PB <- F)
      if (PB == T)
        pb <- txtProgressBar(min = 0, max = dim(animalx.gx)[1], style = 3)
      for (i in 1:dim(animalx.gx)[1]) {
        measurement <- animalx.gx[i,]
        measurement.data <-
          as.numeric(subset(measurement, select = (data0:get(
            paste("data",rec.frames - 1,sep = "")
          ))))
        measurement.headers <-
          subset(measurement, select = (-(data0:get(
            paste("data",rec.frames - 1,sep = "")
          ))))

        min <- round(min(measurement.data))
        max <- round(max(measurement.data))

        if (is.na(fixedB) & is.na(fixedA) & is.na(fixedC))    fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(A=(max-min), B=10, C=min), weights=weights), silent=T)
        if (!is.na(fixedA) & is.na(fixedB) & is.na(fixedC))   fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(B=10, C=min), weights=weights), silent=T)
        if (is.na(fixedA) & !is.na(fixedB) & is.na(fixedC))   fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(A=(max-min), C=min), weights=weights), silent=T)
        if (is.na(fixedA) & is.na(fixedB) & !is.na(fixedC))   fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(A=(max-min), B=10), weights=weights), silent=T)
        if (!is.na(fixedB) & !is.na(fixedA) & is.na(fixedC))    fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(C=min), weights=weights), silent=T)
        if (!is.na(fixedB) & is.na(fixedA) & !is.na(fixedC))    fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(B=10), weights=weights), silent=T)
        if (is.na(fixedB) & !is.na(fixedA) & !is.na(fixedC))    fit <- try(nls(measurement.data ~ (A*exp(-x/B)+C), trace=F,start=list(A=(max-min)), weights=weights), silent=T)


        if (class(fit) != "try-error") {
          if (is.na(fixedA) == T) A <- coef(fit)["A"]
          if (is.na(fixedB) == T) B <- coef(fit)["B"]
          if (is.na(fixedC) == T) C <- coef(fit)["C"]

          corrected.measurement.data <- numeric()
          for (k in 1:length(measurement.data)) {
            y0 <- measurement.data[k]
            y <- y0 - (A*exp(-k/B)+C)
            corrected.measurement.data <- c(corrected.measurement.data, y)
          }

          corrected.measurement <- measurement
          corrected.measurement[((length(measurement)-rec.frames)+1):length(measurement)] <- corrected.measurement.data
          all.corrected.measurements.animalx.gx <- rbind(all.corrected.measurements.animalx.gx, corrected.measurement)

          if (is.null(ylim)) ylims <- c(min,max)

          if (plot == T) {
            plot(x, measurement.data, type="n", col="darkgrey", lwd="1", ylim=ylims, xlab="frames", ylab="response (deltaF/F)", main=paste(measurement$TOdour, "\n", " c:",measurement$NOConc, " g:",measurement$NGloTag, sep=""), cex.main=.7)
            # plot weights
            relative.weights <- weights/weight.BS
            polygon.range <- (max-min)/2
            polygon.weights <- (relative.weights * polygon.range) +min
            polygon.weights[c(1,length(polygon.weights))] <- min(polygon.weights) # set first and last value to lowest level for plotting

            polygon(x, polygon.weights, col=polcol, border=NA, density=-1)

            legend("topleft",legend=c(paste(round(A,2),"*e^(-x/",round(B,2),")+(",round(C,2),")", sep="")), bty="n",cex=.6,pch=-1)
            if(is.null(ylim))  ylims <- c(-1,abs(max-min))
            lines(x, measurement.data, type="l", col="darkgrey", lwd="1")
            curve((A * exp(-x/B)+C), from=1, to=rec.frames, add=T, col="#ee333399", lwd="3")
            plot(x, corrected.measurement.data, type="n", ylim=ylims, col="black", lwd="1", xlab="frames", ylab="response (deltaF/F)", main=paste(measurement$TOdour, " (corrected) \n", " c:",measurement$NOConc, " g:",measurement$NGloTag, sep=""), cex.main=.7)
            abline(h=0, col="#bbbbbbbb", lwd=2)
            lines(x, corrected.measurement.data, type="l", ylim=ylims, col="black", lwd="1")
          }
        } else {
          not.fitted <- rbind(not.fitted, measurement)
          print(paste("ERROR fitting", measurement$Tanimal, "glomerulus:", measurement$NGloTag, measurement$TOdour, measurement$NOConc))

          if (plot == T) {
            if(is.null(ylim)) ylims <- c(min,max)
            plot(x, measurement.data, type="n", col="darkgrey", lwd="1", ylim=ylims, xlab="frames", ylab="response (deltaF/F)", main=paste(measurement$TOdour, "\n", " c:",measurement$NOConc, " g:",measurement$NGloTag, sep=""), cex.main=.7)
            legend("center",legend="ERROR",text.col="#FF000055", bty="n",cex=2, pch=-1)
            plot(x, measurement.data, type="l", col="darkgrey", lwd="1", ylim=ylims, xlab="frames", ylab="response (deltaF/F)", main=paste(measurement$TOdour, "\n", " c:",measurement$NOConc, " g:",measurement$NGloTag, sep=""), cex.main=.7)
          }
        }


  if(PB == T) setTxtProgressBar(pb, i)
  }
  if(PB == T) close(pb)



      all.corrected.measurements.animalx <- rbind(all.corrected.measurements.animalx, all.corrected.measurements.animalx.gx)
    } #end glomeruli loop

    if (plot == T) dev.off()
    all.corrected.measurements <- rbind(all.corrected.measurements, all.corrected.measurements.animalx)
  } #end animal loop

  # clean up "not.fitted" (remove non exosting levels)
  not.fitted <- droplevels(not.fitted)

  #  # plot not fitted data again, NOT WORKING YET if (dim(not.fitted)[1] > 0) {
  #  pdf(file="notfitted_bleachCorrectionReport", paper="a4", height=11.69 ,
  #  width=8.27) par(mfrow=mfrow)
  #
  #  not.fitted.data <- as.numeric(subset(not.fitted,
  #  select=(data0:get(paste("data",rec.frames-1,sep=""))))) not.fitted.headers
  #  <- subset(not.fitted,
  #  select=(-(data0:get(paste("data",rec.frames-1,sep="")))))
  #
  #  for (d in 1:dim(not.fitted)[1]) { plot(x, not.fitted.data[d,], type="l",
  #  col="darkgrey", lwd="1", ylim=c(min,max), xlab="frames", ylab="response
  #  (deltaF/F)", main=paste(not.fitted.headers$TOdour[d], "\n", "
  #  c:",not.fitted.headers$NOConc[d], "
  #  g:",not.fitted.headers$NGloTag[d],"\n",not.fitted.headers$Tanimal[d],
  #  sep=""), cex.main=.7) }
  #
  #  dev.off() }

  parameters <-
    c(
      rec.frames = rec.frames, bleach.begin = bleach.begin, bleach.leaveout =
        paste(bleach.leaveout[1],"to",bleach.leaveout[length(bleach.leaveout)]), weight.BS =
        weight.BS, weight.AS = weight.AS
    )
  result <-
    list(
      bleachcorrected.measurements = all.corrected.measurements, not.fitted =
        not.fitted, weights = weights, parameters = parameters
    )
  result
}

