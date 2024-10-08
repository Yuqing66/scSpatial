% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalization.R
\name{wdist}
\alias{wdist}
\title{Calculate Wasserstein distance}
\usage{
wdist(bins1, bins2, tab = "counts.smooth")
}
\arguments{
\item{bins1, bins2}{An object storing multiple tabs of matrices. Raw counts or after smoothing.}

\item{tab}{The name of a list in bins object.}
}
\description{
Calculate Wasserstein distance between two wide matrices stored in bins1 and bins2 objects.
}
\examples{
wdist()
}
