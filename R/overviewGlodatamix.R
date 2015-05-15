#' overview gloDatamix
#'
#' prints an pdf overview per animal with additionally all mol, air and reference measurements in one plot
#'
#' @param data.object input data.fram (gloDatamix)
#' @param rec.frames number of recorded frames, read from gloDatamix file if not specified
#' @param ref name of your reference, if specified, all references will be plottet in 1st plot
#' @param mol name of your mineral oil control, if specified, all mols will be plottet in 1st plot
#' @param air name of your room-air conrol, if specified, all airs will be plottet in 1st plot
#' @param ylim_lower lower y axis limit, if empty minimum value is taken
#' @param ylim_upper upper y axis limit, if empty maximum value is taken
#' @param mfrow number of rows and columns per page
#' @param glomwise plot glomerulus wise
#' @param odorandglomwise lot odor and glomerulus wise
#' @param stimulusframes only needed when stimulus is different from what is in the gloDatamix file, used to plot a polygon. has to be defined as one vector stim. start and stop frames ("c(20,24,30,36)")
#' @param polcol color of the stimulus polygon
#' @param suffix suffic for filenames
#'
#' @author Daniel Münch <daniel@@muench.bio>
#'
###########################################################
############ Version 1.1 - 01.08.2008, 12:48 ############
###########################################################
### CHANGELOG
###
### 04.08.2010
### - !!!! removed sorting (in plot glomwise and odorandglomwise) recordings according to NRealTime which makes problems when time switches for example from 12am to 1pm
###
### 30.06.2010
### - renamed "glomwise" to "odorandglomwise", added new ("real") glomwise
###
### 20.04.2010
### - added polcol, the color of the polygon
###
### 15.04.2010
### - changed animals <- levels(data.object$Tanimal) to animals <- levels(factor(data.object$Tanimal))
###
### 05.01.2010
###  - added option "glomwise" to plot data glomerulus and odor-wise (for pharma experiments in AL)
###
### 29.10.2009
###  - added glomeruli support
###
### 25.02.2009
###  - added subfolderfolder creation
###
### 29.01.2009
###  - changed frames to rec.frames (recorded frames) as used in other functions
###
### 01.08.2008, 12:48 - v1.0.1
###  - added option to manualy set upper and lower ylim
###     - if not set, limits are calculated out of data
###     - either one or both limits can be set
###
### v1.0   15.07.2008, 15:13
###  - initial version
###
###########################################################

overviewGlodatamix <- function(	data.object, rec.frames = NA, ref = NA, mol = NA, air = NA, ylim_lower = NA, ylim_upper = NA, mfrow = c(4,3), glomwise = F, odorandglomwise = F, stimulusframes  = NA, polcol = "#efefef", suffix = "") {
  # get data.object.name for foldernames&filenames
  data.object.name <- deparse(substitute(data.object))
  data.object.name <- gsub("[^a-zA-Z0-9 :._-]", "", data.object.name)
  if (suffix == "") suffix <- data.object.name    #use data.object.name if suffix is empty

  # create folder for overviews
  oldwd <- getwd()
  folder <- "overviews"
  if (glomwise == T) folder <- paste(folder, "glomwise", sep="_")
  if (odorandglomwise == T) folder <- paste(folder, "odorandglomwise", sep="_")
  folder <- paste(folder,suffix,sep="_")
  #ifelse (is.na(suffix)==T, folder <- paste("overviews","_",Sys.time(),sep=""), folder <- paste("overviews","_",suffix,"_",Sys.time(),sep=""))
  dir.create(path=folder)
  setwd(paste(oldwd,folder, sep="/"))

  # if rec.frame is not defined get it from first row of dataset
  if (is.na(rec.frames)) {
    rec.frames <- data.object$NNoFrames[1]
    print(paste("number of frames was set to",rec.frames))
  }

  # get animal names and count
  animals <- levels(factor(data.object$Tanimal))

  # get position of datapoints within data.frame
  firstframe <- length(data.object) - rec.frames + 1
  lastframe  <- length(data.object)

  # calculate animal overviews
  for (i in 1:length(animals)) {
    animalx <- subset(data.object, Tanimal==animals[i])
    recordings <- animalx$T_dbb1
    glomeruli <- levels(droplevels(as.factor(animalx$NGloTag)))
    color.glom <- rainbow(length(glomeruli))
    firstframex <- length(animalx) - rec.frames + 1
    lastframex <- length(animalx)

    ### calculate ylim if not set


    ## old version ## ylim=c(min(animalx[,firstframex:lastframex]),max(animalx[,firstframex:lastframex]))

    if (is.na(ylim_lower)){
    ylimlower <- min(animalx[,firstframex:lastframex])
    }	else ylimlower <- ylim_lower

    if (is.na(ylim_upper)){
    ylimupper <- max(animalx[,firstframex:lastframex])
    }	else ylimupper <- ylim_upper

    ylim=c(ylimlower,ylimupper)


    pdf(file=paste(animals[i],"_",suffix,"_overview.pdf",sep=""), title=animals[i], paper="a4", height=11.69 , width=8.27)
    par(mfrow=mfrow)

    # convert stimulusframes into vector for polygon
    if(is.na(stimulusframes[1])) stimulusframes <- c(data.object$NStim_on[1], data.object$NStim_off[1], data.object$Nstim2ON[1], data.object$Nstim2OFF[1])

    polygon.x <- rep(stimulusframes, each=2)
    polygon.y <- rep(c(ylim[1], ylim[2], ylim[2], ylim[1]), length(stimulusframes)/2)


    ### plot references
    if (!is.na(ref)) {
      refx <- subset(data.object, TOdour==ref&Tanimal==animals[i])
      firstframe.ref <- length(refx) - rec.frames + 1
      lastframe.ref <- length(refx)
      j.ref <- 1:dim(refx)[1]
      color.ref <- rainbow(length(j.ref))

      plot(1:rec.frames, refx[1,firstframe.ref:lastframe.ref],type="n",main=refx$TOdour[1], ylim=ylim, xlab="frames", ylab="response (deltaF/F)")

      polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)

      for (j in j.ref) {
        lines(1:rec.frames, refx[j,firstframe.ref:lastframe.ref],type="l", col=color.ref[j])
      }
      if (j != 0) {
        legend("topright", legend=1:j, fill=color.ref[1:j], bty="n", ncol=2)
      }
    }
    ### plot mol
    if (!is.na(mol)) {
      refx <- subset(data.object, TOdour==mol&Tanimal==animals[i])
      firstframe.ref <- length(refx) - rec.frames + 1
      lastframe.ref <- length(refx)
      j.ref <- 1:dim(refx)[1]
      color.ref <- rainbow(length(j.ref))

      plot(1:rec.frames, refx[1,firstframe.ref:lastframe.ref],type="n",main=refx$TOdour[1], ylim=ylim, xlab="frames", ylab="response (deltaF/F)")
      polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)
      for (j in j.ref) lines(1:rec.frames, refx[j,firstframe.ref:lastframe.ref],type="l", col=color.ref[j])
      if (j != 0) legend("topright", legend=1:j, fill=color.ref[1:j], bty="n", ncol=2)
    }
    ### plot air
    if (!is.na(air)) {
      refx <- subset(data.object, TOdour==air&Tanimal==animals[i])
      firstframe.ref <- length(refx) - rec.frames + 1
      lastframe.ref <- length(refx)
      j.ref <- 1:dim(refx)[1]
      color.ref <- rainbow(length(j.ref))

      plot(1:rec.frames, refx[1,firstframe.ref:lastframe.ref],type="n",main=refx$TOdour[1], ylim=ylim, xlab="frames", ylab="response (deltaF/F)")
      polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)
      for (j in j.ref) lines(1:rec.frames, refx[j,firstframe.ref:lastframe.ref],type="l", col=color.ref[j])
      if (j != 0) legend("topright", legend=1:j, fill=color.ref[1:j], bty="n", ncol=2)
    }

    ### plot measurements:
    if (odorandglomwise == F && glomwise == F)  {
      for (k in 1:length(recordings)) {
      recordingx <- subset(animalx, T_dbb1==recordings[k])
      plot(1:rec.frames,recordingx[1,firstframex:lastframex],type="n",main=paste(recordingx$TOdour[1],"(",recordingx$NOConc[1],"), ",recordingx$TPharma[1]), ylim=ylim, xlab="frames", ylab="response (deltaF/F)", col="red")
      polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)
      for (l in 1:dim(recordingx)[1]) lines(1:rec.frames,recordingx[l,firstframex:lastframex],type="l", col=color.glom[l])
      #	 	lines(x=animalx$NStim_on[k]:animalx$NStim_off[k],y=rep(ylim[1], length(animalx$NStim_on[k]:animalx$NStim_off[k])), type="l", col="black")
      #		lines(x=animalx$Nstim2ON[k]:animalx$Nstim2OFF[k],y=rep(ylim[1], length(animalx$Nstim2ON[k]:animalx$Nstim2OFF[k])), type="l", col="black")
      }
    }

    if (odorandglomwise == T) {
      odorsx <- levels(factor(animalx$TOdour))
      for (k in 1:length(glomeruli)) {
        for (l in 1:length(odorsx)) {
          glomxodorx <- subset(animalx, NGloTag==glomeruli[k] & TOdour==odorsx[l])
          conc.glomxodorx <- levels(factor(glomxodorx$NOConc))
          for (c in 1:length(conc.glomxodorx)) {
            concx <- subset(glomxodorx, NOConc == conc.glomxodorx[c])

            n.recordings <- dim(concx)[1]
            col.recordings <- rainbow(n.recordings)
            PharmaNames <- concx$TPharma

            plot(1:rec.frames,concx[1,firstframex:lastframex],type="n",main=paste("glom.:",glomeruli[k],", ",concx$TOdour[1],"(",concx$NOConc[1],")"), ylim=ylim, xlab="frames", ylab="response (deltaF/F)", col="red")
            polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)
            for (m in 1:n.recordings) lines(1:rec.frames,concx[m,firstframex:lastframex],type="l", col=col.recordings[m])
            legend("topright", legend=PharmaNames, fill=col.recordings[1:m], bty="n", ncol=1)
          }
        }
      }
    }

    if (glomwise == T) {
      odorsx <- levels(factor(animalx$TOdour))
      for (k in 1:length(glomeruli)) {
        glomx          <- subset(animalx, NGloTag==glomeruli[k])
        n.recordings   <- dim(glomx)[1]
        col.recordings <- rainbow(n.recordings)
        OdorNames      <- glomx$TOdour

        plot(1:rec.frames,glomx[1,firstframex:lastframex],type="n",main=paste("glom.:",glomeruli[k]), ylim=ylim, xlab="frames", ylab="response (deltaF/F)", col="red")
        polygon(y=polygon.y, x=polygon.x, density=-1, col=polcol, border=NA)
        for (m in 1:n.recordings)lines(1:rec.frames,glomx[m,firstframex:lastframex],type="l", col=col.recordings[m])
        legend("topright", legend=OdorNames, fill=col.recordings[1:m], bty="n", ncol=1)
      }
    }

    dev.off()

    print(paste(animals[i],"_overview.pdf......created!"))
  }
  setwd(oldwd)
  print("..........done!")
}
