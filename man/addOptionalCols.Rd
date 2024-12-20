% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic.R
\name{addOptionalCols}
\alias{addOptionalCols}
\title{Add optional columns to a data frame.}
\usage{
addOptionalCols(df, ...)
}
\arguments{
\item{df}{A data frame.}

\item{...}{Vectors to add to the data frame. Format column_name=vector. Use vector name if column_name is not provided.}
}
\value{
A data frame with optional columns added.
}
\description{
Add optional columns to a data frame if not NULL.
}
\examples{
m1 <- data.frame(a=1:5, b=letters[1:5])
c <- 2:6
m2 <- addOptionalCols(df=m1, c, d=3:7)
}
