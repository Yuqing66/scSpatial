% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalization.R
\name{getCoords.transcript}
\alias{getCoords.transcript}
\title{Get the coordinates of transcripts in a Seurat object}
\usage{
getCoords.transcript(object, transcript, image.name)
}
\arguments{
\item{object}{Seurat object.}

\item{transcript}{The name of the transcript.}

\item{image.name}{The name of the image fov in the Seurat object.}
}
\description{
Extract the coordinates of transcripts in a fov.
}
\examples{
coords <- getCoords.transcript(srt, transcript="IFNB1", image.name="UV109fov1")
}
