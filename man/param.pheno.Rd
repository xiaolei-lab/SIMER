% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Parameter.R
\name{param.pheno}
\alias{param.pheno}
\title{Phenotype parameters generator}
\usage{
param.pheno(SP = NULL, ...)
}
\arguments{
\item{SP}{a list of all simulation parameters.}

\item{...}{one or more parameter(s) for phenotype simulation.}
}
\value{
the function returns a list containing
\describe{
\item{$pheno$pop}{the population information containing environmental factors and other effects.}
\item{$pheno$pop.ind}{the number of individuals in the base population.}
\item{$pheno$pop.rep}{the repeated times of repeated records.}
\item{$pheno$pop.rep.bal}{whether repeated records are balanced.}
\item{$pheno$pop.env}{a list of environmental factors setting.}
\item{$pheno$phe.type}{a list of phenotype types.}
\item{$pheno$phe.model}{a list of genetic model of phenotype such as "T1 = A + E".}
\item{$pheno$phe.h2A}{a list of additive heritability.}
\item{$pheno$phe.h2D}{a list of dominant heritability.}
\item{$pheno$phe.h2GxG}{a list of GxG interaction heritability.}
\item{$pheno$phe.h2GxE}{a list of GxE interaction heritability.}
\item{$pheno$phe.h2PE}{a list of permanent environmental heritability.}
\item{$pheno$phe.var}{a list of phenotype variance.}
\item{$pheno$phe.corA}{the additive genetic correlation matrix.}
\item{$pheno$phe.corD}{the dominant genetic correlation matrix.}
\item{$pheno$phe.corGxG}{the GxG genetic correlation matrix.}
\item{$pheno$phe.corPE}{the permanent environmental correlation matrix.}
\item{$pheno$phe.corE}{the residual correlation matrix.}
}
}
\description{
Generate parameters for phenotype data simulation.
}
\details{
Build date: Feb 21, 2022
Last update: Jul 4, 2022
}
\examples{
SP <- param.pheno(phe.model = list(tr1 = "T1 = A + E"))
str(SP)
}
\author{
Dong Yin
}
