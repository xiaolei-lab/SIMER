% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Genotype.R
\name{geno.cvt1}
\alias{geno.cvt1}
\title{Genotype code convertor 1}
\usage{
geno.cvt1(pop.geno, ncpus = 0)
}
\arguments{
\item{pop.geno}{genotype matrix of (0, 1).}

\item{ncpus}{the number of threads used, if NULL, (logical core number - 1) is automatically used.}
}
\value{
genotype matrix of (0, 1, 2).
}
\description{
Convert genotype matrix from (0, 1) to (0, 1, 2).
}
\details{
Build date: Nov 14, 2018
Last update: Jan 30, 2025
}
\examples{
\donttest{
library(bigmemory)
options(bigmemory.typecast.warning=FALSE)
pop.geno <- matrix(sample(c(0, 1), 16, replace = TRUE), 4, 4)
pop.geno[]
bigmat <- geno.cvt1(pop.geno)
bigmat[]
pop.geno <- as.big.matrix(pop.geno, type = 'char')
bigmat <- geno.cvt1(pop.geno)
bigmat[]
}
}
\author{
Dong Yin
}
