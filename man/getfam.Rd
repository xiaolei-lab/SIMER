% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Reproduces.R
\name{getfam}
\alias{getfam}
\title{Family index and within-family index}
\usage{
getfam(sir, dam, fam.op, mode = c("pat", "mat", "pm"))
}
\arguments{
\item{sir}{the indice of sires.}

\item{dam}{the indice of dams.}

\item{fam.op}{the initial index of family indice.}

\item{mode}{"pat": paternal mode; "mat": maternal mode; "pm": paternal and maternal mode.}
}
\value{
a matrix with family indice and within-family indice.
}
\description{
Get indice of family and within-family
}
\details{
Build date: Nov 14, 2018
Last update: Apr 30, 2022
}
\examples{
s <- c(0, 0, 0, 0, 1, 3, 3, 1, 5, 7, 5, 7, 1, 3, 5, 7)
d <- c(0, 0, 0, 0, 2, 4, 4, 2, 6, 8, 8, 6, 6, 8, 4, 8)
fam <- getfam(sir = s, dam = d, fam.op = 1, mode = "pm")
fam
}
\author{
Dong Yin
}
