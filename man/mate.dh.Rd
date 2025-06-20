% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Reproduces.R
\name{mate.dh}
\alias{mate.dh}
\title{Doubled haploid}
\usage{
mate.dh(SP, ncpus = 0, verbose = TRUE)
}
\arguments{
\item{SP}{a list of all simulation parameters.}

\item{ncpus}{the number of threads used, if NULL, (logical core number - 1) is automatically used.}

\item{verbose}{whether to print detail.}
}
\value{
the function returns a list containing
\describe{
\item{$reprod$pop.gen}{the generations of simulated population.}
\item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', 'assort', 'disassort', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
\item{$reprod$sex.rate}{the sex ratio of simulated population.}
\item{$reprod$prog}{the progeny number of an individual.}
\item{$geno}{a list of genotype simulation parameters.}
\item{$pheno}{a list of phenotype simulation parameters.}
}
}
\description{
Produce individuals by doubled haploid.
}
\details{
Build date: Nov 14, 2018
Last update: Jan 28, 2025
}
\examples{
\donttest{
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = list(tr1 = 10))
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "ind")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "dh")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run doubled haploid
SP <- mate.dh(SP)
}
}
\author{
Dong Yin
}
