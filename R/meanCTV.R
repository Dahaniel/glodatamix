###########################################################
############ Version 0.1 - 10.11.2010 #####################
###########################################################
### INFO
### Function to calculate mean, median... from ctvs, needs 
### output from ctvs.R as input (GloDatamix with ctv column
###
###########################################################

### 110212: added paired ttest against mol

meanCTV <- function(data.object,test=F, nullref="MOL", plev=0.05, ctvcolumnname = "ctv", shift = F)
    {
    
  # select subset of data with ctv column of interest
  data.object <- subset(data.object, select = c("NGloTag", "NOConc", "TOdour", "Tanimal", ctvcolumnname))
  colnames(data.object)[which(colnames(data.object) == ctvcolumnname)] <- "ctv"
  
    if (test == T) 
        {
        data.object <- cbind(data.object, null.mean.withinanimal=rep(NA, dim(data.object)[1]))
        animals <- levels(data.object$Tanimal)
        for (l in 1:length(animals))
            {
            animalx <- subset(data.object, Tanimal==animals[l] & TOdour == nullref)
            mean.nullref <- mean(animalx$ctv)
            animalx.pos <- which(data.object$Tanimal == animals[l])
            data.object[animalx.pos,"null.mean.withinanimal"] <- mean.nullref
            }  
        }
    
    
    
    result <- data.frame()
    
    glomeruli <- levels(factor(data.object$NGloTag))

    for (i in 1:length(glomeruli))
        {
        glomerulusx <- subset(data.object, NGloTag == glomeruli[i])
        odors <- levels(factor(glomerulusx$TOdour))
        
        for (j in 1:length(odors))
            {
            odorx <- subset(glomerulusx, TOdour==odors[j])
            concentrations <- levels(factor(odorx$NOConc))
            
            for (k in 1:length(concentrations))
                {
                concentrationx <- subset(odorx, NOConc == concentrations[k])
                ctvsx <- concentrationx$ctv
                
                pvalue.t <- NA
                sig.t <- NA
                pvalue.w <- NA
                sig.w <- NA
                if (test==T)
                    {
                    ttest <- t.test(concentrationx$ctv, concentrationx$null.mean.withinanimal, paired=T)
                    pvalue.t <- ttest$p.value
                    sig.t <- pvalue.t < plev 
                    wtest <- wilcox.test(concentrationx$ctv, concentrationx$null.mean.withinanimal, paired=T)
                    pvalue.w <- wtest$p.value
                    sig.w <- pvalue.w < plev 
                    }
                
                ctvs.n        <- length(na.omit(ctvsx))
                ctvs.mean    <- mean(ctvsx, na.rm=T)
                ctvs.median  <- median(ctvsx, na.rm=T)
                ctvs.sd      <- as.numeric(sd(ctvsx, na.rm=T))  #output is a dataframe which makes renaming of the col title impossible thus convert to numeric before...
                ctvs.sem     <- ctvs.sd/sqrt(ctvs.n)
                ctvs.var     <- var(ctvsx, na.rm=T)
                ctvs.quant   <- quantile(ctvsx,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=T)
                ctvs.quant2.5 <- as.numeric(ctvs.quant["2.5%"])
                ctvs.quant97.5 <- as.numeric(ctvs.quant["97.5%"])
                ctvs.quant25 <- as.numeric(ctvs.quant["25%"])
                ctvs.quant75 <- as.numeric(ctvs.quant["75%"])
                ctvs.min     <- as.numeric(ctvs.quant["0%"])
                ctvs.max     <- as.numeric(ctvs.quant["100%"])
                
                
                resultx <-  data.frame("glomerulus"=glomeruli[i], "odor" = odors[j], "concentration" = as.numeric(concentrationx$NOConc[1]), #take from concentrationx not from concentrations[k] otherwise it will be saved as factor
                            "n" = ctvs.n, "mean.response" = ctvs.mean, "median.response" = ctvs.median, "sd" = ctvs.sd, "sem" = ctvs.sem, "variance" = ctvs.var, "quant2.5" = ctvs.quant2.5, "quant97.5" = ctvs.quant97.5, "quant25" = ctvs.quant25, "quant75" = ctvs.quant75, "min" = ctvs.min, "max" = ctvs.max, "plevel" = plev, "pvalue.t" = pvalue.t, "significant.t" = sig.t, "pvalue.w" = pvalue.w, "significant.w" = sig.w)

                result <- rbind(result, resultx)
                }
            }
        }
  
  if (shift == T)
  {
    for(i in 1:length(glomeruli))
    {
      result[result$glomerulus==glomeruli[i],"mean.response"] <- result[result$glomerulus==glomeruli[i],"mean.response"] - result[result$glomerulus==glomeruli[i] & result$odor==nullref,"mean.response"]
      
      result[result$glomerulus==glomeruli[i],"median.response"] <- result[result$glomerulus==glomeruli[i],"median.response"] - result[result$glomerulus==glomeruli[i] & result$odor==nullref,"median.response"]
      
#       null.mean <- result[result$odor==nullref,]$mean.response
#       null.median <- result[result$odor==nullref,]$median.response
#       
#       result$mean.response <- result$mean.response- null.mean
#       result$median.response <- result$median.response- null.median
      
      cat("\nData in \"mean.response\" and \"median.response\" were shifted by corresponding", nullref, "value. Other values like quantiles WERE NOT SHIFTED!\n")
    }
    
  }
    
    return(result)
    }
