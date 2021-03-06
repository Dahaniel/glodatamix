% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/melt.gdm.R
\name{melt.gdm}
\alias{melt.gdm}
\title{melt.gdm}
\usage{
melt.gdm(gdm)
}
\arguments{
\item{gdm}{the data.frame containing the gloDatamix data (from read.gdm())}
}
\value{
a "melted", long version of the gdm
}
\description{
melt.gdm
}
\examples{
\dontrun{
data <- read.gdm("./80109-final_all/", header=F)
melt.gdm(data)
}
}
\author{
Daniel Münch <daniel@muench.bio>
}

