#' melt.gdm
#'
#' @param gdm the data.frame containing the gloDatamix data (from read.gdm())
#'
#' @return a "melted", long version of the gdm
#' @export
#'
#' @author Daniel Münch <daniel@@muench.bio>
#'
#' @examples \dontrun{
#' data <- read.gdm("./80109-final_all/", header=F)
#' melt.gdm(data)
#' }
melt.gdm <- function(gdm) {
  require(tidyr)
  framerange <- range(grep("data", names(gdm)))
  gdm <- gather(data = gdm, key = frame, value = value, framerange[1]:framerange[2]) %>%
    mutate(frame = as.numeric(substr(frame, start = 5, stop = 99)))
}
