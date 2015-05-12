##########
# Function that fits a curve over reference odors and normalizes to it
##########

##########
# CHANGELOG
#
# 08.03.2011
# adding possibility to subtract nullreference
# !!! stopped doing that, as nullref has to be subtracted as ctv, writing new function, this here will be for normalizing traces to ref only!
##########

##########
# Scheme:
# 1. subset of reference
# 2. ctv of reference
# 3. fit curve over ctv/time
#
# A-4. divide every timepoint of every measurement through corresponding point on fitted curve -> normalized curves
# A-5. do ctv for measurements -> normalized ctvs
#
# B-4. do ctv for measurements
# B-5. divide every ctv through corresponding point on fitted curve
##########

normalizeToReference <- function(
                                    data.object, 
                                    rec.frames=NA, 
                                    reference.odor, 
                                    ctv.type, 
                                    bg.firstframe=NA, 
                                    bg.lastframe=NA, 
                                    sig.firstframe, 
                                    sig.lastframe, 
                                    pdf=T, 
                                    mode="mean", 
                                    rescale.data = F,      #if FALSE: reference set to 1, if set to TRUE: values are multiplied with first reference value afterwards
                                    ctvreport=F, 
                                    glom=NA,
                                    suffix="", 
                                    plot = T
                                    )
    {

#get data.object name
data.object.name <- deparse(substitute(data.object))
data.object.name <- gsub("[^a-zA-Z0-9 :._-]", "", data.object.name)
if (suffix == "") suffix <- data.object.name    #use data.object.name if suffix is empty



# if rec.frame is not defined get it from first row of dataset
    if (is.na(rec.frames)){
     rec.frames <- data.object$NNoFrames[1]
     cat(paste("\nnumber of frames was automatically set to",rec.frames,"\n\r"))
    }
    # print normalization to pdf
    if ((pdf==T & mode=="linFit") | (pdf==T & mode=="linFitEachGlom")) {
    pdf(file=paste("fitsForNormalization_",suffix,".pdf",sep=""), paper="a4", height=11.69 , width=8.27)
    par(mfrow=c(4,3))
    }

    # Codesammlung:
    animals <- levels(factor(data.object$Tanimal))
    normed.data <- data.frame()
    for (i in 1:length(animals)) 
        {
        animalx <- subset(data.object, Tanimal==animals[i])
        if (is.na(glom) == F) 
            {
                    animalx.ref <- subset(animalx, TOdour==reference.odor & NGloTag==glom)
            } else  animalx.ref <- subset(animalx, TOdour==reference.odor)
         
         
         # animalx.ref.ctvdata <- subset(animalx.ref, select=c(NRealTime,data0:get(paste("data",rec.frames-1,sep="")))) # removed, instead keep all the headers
         
         # calculate ctv
         # ctv auf eigene function ausgelagert
         ctv <- ctvs(animalx.ref,ctv.type,rec.frames,bg.firstframe,bg.lastframe,sig.firstframe,sig.lastframe,report=ctvreport)  # give data to "ctvs" function, leave timepoints out
         
         # ctv <- ctvs(animalx.ref.ctvdata[,-1],ctv.number,bg.firstframe,bg.lastframe,sig.firstframe,sig.lastframe)  # give data to "ctvs" function, leave timepoints out
         # ctv <- cbind(animalx.ref.ctvdata[,1],ctv)

        ## MODE: linFit
        if (mode == "linFit") 
            {
              if (dim(ctv)[1] != 1)
              {
                fit <- lm(ctv$ctv~ctv$NRealTime)
      	        if (plot == T) {
      	          plot(ctv$NRealTime,ctv$ctv, main=levels(data.object$Tanimal)[i], xlab="time [min]", ylab=expression(Delta*F/F))
      	          abline(fit)
      	        }

      	        #to get corresponding points on line
      	        #fit$coef[2]*theXvalue+fit$coef[1]
      
      	        animalx.normed <- data.frame()
      	        for (j in 1:dim(animalx)[1]) 
      	            {
      	            normto <- fit$coef[2]*animalx[j,]$NRealTime+fit$coef[1]
      	            normed <- subset(animalx[j,], select=(data0:get(paste("data",rec.frames-1,sep=""))))/normto
                    if(rescale.data == T) normed <- normed * arrange(ctv, NRealTime)$ctv[1] #rescale to first reference, sort dataframe by NRealTime first to make sure that the first measurement is selected
      	            animalx.normed <- rbind(animalx.normed,normed)
      	            }
                  cat("normalized using linFit\n\r")
              }   else 
                    {
                     if (plot == T) 
                       plot(ctv$NRealTime,ctv$ctv, main=levels(data.object$Tanimal)[i], xlab="time [min]", ylab=expression(Delta*F/F))
                    
                     normto <- ctv$ctv
                     normed <- subset(animalx, select=(data0:get(paste("data",rec.frames-1,sep=""))))/normto
                     if(rescale.data == T) normed <- normed * arrange(ctv, NRealTime)$ctv[1] #rescale to first reference, sort dataframe by NRealTime first to make sure that the first measurement is selected
                     animalx.normed <- normed
    
                  cat(paste("animal: ",animals[i]," contained only one reference, normalized to this value (",normto,")\n\r",sep=""))
                     
                    
                    }

        
        
        
          }


        ## MODE: mean
        if (mode == "mean") 
            {
	        animalx.normed <- data.frame()
	        for (j in 1:dim(animalx)[1]) 
	            {
	            normto <- mean(ctv$ctv)
	            normed <- subset(animalx[j,], select=(data0:get(paste("data",rec.frames-1,sep=""))))/normto
	            animalx.normed <- rbind(animalx.normed,normed)
	            }
            cat("normalized using mean of",normto,"\n\r")
            }



        if (mode == "linFitEachGlom") ### NOT WORKING YET; RESULTS LOOK BAD! MUST BE SOME ERROR
             {
             cat("Normalizing using a linear fit over references within each glomerulus\n")
             glomeruli <- levels(as.factor(animalx$NGloTag))
             animalx.normed <- data.frame()
             for (k in 1:length(glomeruli))
                {
	            glomx <- subset(animalx, NGloTag==glomeruli[k])
	            glomx.ref <- subset(glomx, TOdour==reference.odor)
	            ctv <- ctvs(glomx.ref,ctv.type,rec.frames,bg.firstframe,bg.lastframe,sig.firstframe,sig.lastframe,report=ctvreport, suffix=paste(animals[i],glomeruli[k],sep="_"))  # give data to "ctvs" function, leave timepoints out
	             
	            fit <- lm(ctv$ctv~ctv$NRealTime)
	            if (plot == T) {
	              plot(ctv$NRealTime,ctv$ctv, main=paste(animals[i],glomeruli[k]), xlab="time [min]", ylab=expression(Delta*F/F))
	              abline(fit)
	            }
              #to get corresponding points on line
	            #fit$coef[2]*theXvalue+fit$coef[1]

	            glomx.normed <- data.frame()
	            for (j in 1:dim(glomx)[1]) 
	                {
	                normto <- fit$coef[2]*glomx[j,]$NRealTime+fit$coef[1]
	                normeddata <- subset(glomx[j,], select=(data0:get(paste("data",rec.frames-1,sep=""))))/normto
	                if(rescale.data == T) normeddata <- normeddata * arrange(ctv, NRealTime)$ctv[1] #rescale to first reference, sort dataframe by NRealTime first to make sure that the first measurement is selected
	                normed <- cbind(glomx[j,1:(length(glomx)-rec.frames)],normeddata)
	                glomx.normed <- rbind(glomx.normed,normed)
	                }
	
                animalx.normed <- rbind(animalx.normed, glomx.normed)  
	            }
	
            cat("normalized using linFit\n\r") 
            }


        if (mode != "linFitEachGlom")
            {
            animalx.normed <- cbind(animalx[,1:(length(animalx)-rec.frames)],animalx.normed) ### ACHTUNG FEHLER BEI LINFITEACHGLOM! Daten sind umsortiert, also werden die falschen header angebunden!!!
            }


        normed.data <- rbind(normed.data, animalx.normed)
        cat(paste(dim(animalx)[1], " measurements in ", levels(data.object$Tanimal)[i],"...normalized to ", reference.odor, "\n\r", sep=""))
        
        if (is.na(glom) == F) cat(paste("normalized to glomerulus", glom, "\n\r", sep=""))

    }

    if ((pdf==T & mode=="linFit") | (pdf==T & mode=="linFitEachGlom")) 
        {
        dev.off()
        }
    if (rescale.data) cat("rescaled data to first reference value\n\r")
    print("DONE!")
    return(normed.data)
    }
