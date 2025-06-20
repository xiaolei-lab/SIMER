% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Genotype.R
\name{annotation}
\alias{annotation}
\title{Annotation simulation}
\usage{
annotation(SP, verbose = TRUE)
}
\arguments{
\item{SP}{a list of all simulation parameters.}

\item{verbose}{whether to print detail.}
}
\value{
the function returns a list containing
\describe{
\item{$map$pop.map}{the map data with annotation information.}
\item{$map$species}{the species of genetic map, which can be "arabidopsis", "cattle", "chicken", "dog", "horse", "human", "maize", "mice", "pig", and "rice".}
\item{$map$pop.marker}{the number of markers.}
\item{$map$num.chr}{the number of chromosomes.}
\item{$map$len.chr}{the length of chromosomes.}
\item{$map$qtn.model}{the genetic model of QTN such as "A + D".}
\item{$map$qtn.index}{the QTN index for each trait.}
\item{$map$qtn.num}{the QTN number for (each group in) each trait.}
\item{$map$qtn.dist}{the QTN distribution containing "norm", "geom", "gamma" or "beta".}
\item{$map$qtn.var}{the variances for normal distribution.}
\item{$map$qtn.prob}{the probability of success for geometric distribution.}
\item{$map$qtn.shape}{the shape parameter for gamma distribution.}
\item{$map$qtn.scale}{the scale parameter for gamma distribution.}
\item{$map$qtn.shape1}{the shape1 parameter for beta distribution.}
\item{$map$qtn.shape2}{the shape2 parameter for beta distribution.}
\item{$map$qtn.ncp}{the ncp parameter for beta distribution.}
\item{$map$qtn.spot}{the QTN distribution probability in each block.}
\item{$map$len.block}{the block length.}
\item{$map$maf}{the maf threshold, markers less than this threshold will be exclude.}
\item{$map$recom.spot}{whether to generate recombination events.}
\item{$map$range.hot}{the recombination times range in the hot spot.}
\item{$map$range.cold}{the recombination times range in the cold spot.}
}
}
\description{
Generating a map with annotation information
}
\details{
Build date: Nov 14, 2018
Last update: Jul 10, 2022
}
\examples{
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = list(tr1 = 10))

# Run annotation simulation
SP <- annotation(SP)
}
\author{
Dong Yin
}
