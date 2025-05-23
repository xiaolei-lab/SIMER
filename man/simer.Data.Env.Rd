% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Data.R
\name{simer.Data.Env}
\alias{simer.Data.Env}
\title{Environmental factor selection}
\usage{
simer.Data.Env(
  jsonList = NULL,
  hiblupPath = "",
  header = TRUE,
  sep = "\\t",
  ncpus = 10,
  verbose = TRUE
)
}
\arguments{
\item{jsonList}{the list of environmental factor selection parameters.}

\item{hiblupPath}{the path of HIBLUP software.}

\item{header}{the header of file.}

\item{sep}{the separator of file.}

\item{ncpus}{the number of threads used, if NULL, (logical core number - 1) is automatically used.}

\item{verbose}{whether to print detail.}
}
\value{
the function returns a list containing
\describe{
\item{$genotype}{the path of genotype data.}
\item{$pedigree}{the filename of pedigree data.}
\item{$selection_index}{the selection index for all traits.}
\item{$breeding_value_index}{the breeding value index for all traits.}
\item{$quality_control_plan}{a list of parameters for data quality control.}
\item{$breeding_plan}{a list of parameters for genetic evaluation.}
}
}
\description{
To find appropriate fixed effects, covariates, and random effects.
}
\details{
Build date: July 17, 2021
Last update: Apr 28, 2022
}
\examples{
# Read JSON file
jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)

\dontrun{
# It needs "hiblup" solfware
jsonList <- simer.Data.Env(jsonList = jsonList)
}
}
\author{
Dong Yin
}
