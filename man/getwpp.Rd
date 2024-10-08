% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalization.R
\name{getwpp}
\alias{getwpp}
\title{Constructor for the wpp Class}
\usage{
getwpp(bins, tab)
}
\arguments{
\item{bins}{An object storing multiple tabs of matrices. Raw counts or after smoothing.}

\item{tab}{The name of a list in bins object.}
}
\description{
Construct an object of class "wpp" from a wide matrix stored in bins object.
}
\examples{
getwpp(bins, tab="counts.smooth")
}
