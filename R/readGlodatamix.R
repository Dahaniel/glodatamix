readGlodatamix <- function(directory,
                           header= TRUE,
                           input = "gloDatamix" # gloDatamix for gloDatamix-files, "csv" for knime output
                           )
{

###### read glodatamix inside the directory, and merge them into a new assigned data.frame in R.
##example:	a<-readGlodatamix("./80109-final_all/",header=F)

wd <- getwd()
setwd(directory)
if (input == "gloDatamix") nfile<-dir(pattern="(gloDatamix)$",ignore.case=T)
if (input == "csv")        nfile<-dir(pattern="(csv)$",ignore.case=T)

if (length(nfile) < 1 & input == "gloDatamix") cat(paste("\n\nWARNING!!!\n No *.gloDatamix files found @", getwd(),", please check folder and filenames (must end with gloDatamix)\n\n"))
if (length(nfile) < 1 & input == "csv") cat(paste("\n\nWARNING!!!\n No *.csv files found @", getwd(),", please check folder and filenames (must end with csv)\n\n"))


combdata <- numeric()
for (i in 1:length(nfile)) 
{
      if (input == "gloDatamix") pcd <- read.table(nfile[i], sep="\t", header=header)
      if (input == "csv")        pcd <- read.csv(nfile[i], header=header)
      combdata <- rbind(combdata, pcd)
      print(nfile[i])
}

setwd(wd)
return(combdata)
}
