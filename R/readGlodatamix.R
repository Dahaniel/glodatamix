#' read.gdm
#'
#' @param path character, location of the *.gloDatamix files
#' @param header logical, gloDataMix files with or without header?
#' @param input character, 'gloDatamix' or 'csv'
#'
#' @return concatenates all gloDatamix files from a location in one data.frame
#' @author Daniel Münch <daniel@@muench.bio>
#'
#' @export
#' @examples \dontrun{
#' data <- read.gdm("./80109-final_all/", header=F)
#' }
read.gdm <- function(path,
                     header = TRUE,
                     input  = "gloDatamix" # gloDatamix for gloDatamix-files, "csv" for knime output
                     ) {
  wd <- getwd()
  setwd(path)
  if (input == "gloDatamix") nfile <- dir(pattern = "(gloDatamix)$", ignore.case = T)
  if (input == "csv")        nfile <- dir(pattern = "(csv)$", ignore.case = T)

  if (length(nfile) < 1 & input == "gloDatamix") warning(paste("No *.gloDatamix files found @", getwd(),", please check folder and filenames (must end with gloDatamix)\n\n"))
  if (length(nfile) < 1 & input == "csv") warning(paste("No *.csv files found @", getwd(),", please check folder and filenames (must end with csv)\n\n"))


  combdata <- numeric()
  for (i in 1:length(nfile)) {
    if (input == "gloDatamix") tmp <- read.table(nfile[i], sep="\t", header=header)
    if (input == "csv")        tmp <- read.csv(nfile[i], header=header)
    combdata <- rbind(combdata, tmp)
    message(nfile[i])
  }

  setwd(wd)
  return(combdata)
}
