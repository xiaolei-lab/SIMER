% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simer.Data.R
\name{simer.Data.Bfile2MVP}
\alias{simer.Data.Bfile2MVP}
\title{simer.Data.Bfile2MVP: To transform plink binary data to MVP package}
\usage{
simer.Data.Bfile2MVP(
  bfile,
  out = "simer",
  maxLine = 10000,
  type.geno = "char",
  threads = 10,
  verbose = TRUE
)
}
\arguments{
\item{bfile}{Genotype in binary format (.bed, .bim, .fam).}

\item{out}{the name of output file.}

\item{maxLine}{the max number of line to write to big matrix for each loop.}

\item{type.geno}{the type of genotype elements.}

\item{threads}{number of thread for transforming.}

\item{verbose}{whether to print the reminder.}
}
\value{
number of individuals and markers.
Output files:
genotype.desc, genotype.bin: genotype file in bigmemory format
phenotype.phe: ordered phenotype file, same taxa order with genotype file
map.map: SNP information
}
\description{
transforming plink binary data to MVP package.
}
\details{
Build date: Sep 12, 2018
Last update: Dec 28, 2024
}
\examples{
\donttest{
# Get bfile path
bfilePath <- file.path(system.file("extdata", "02plinkb", package = "simer"), "demo")

# Data converting
simer.Data.Bfile2MVP(bfilePath, tempfile("outfile"))
}
}
\author{
Haohao Zhang and Dong Yin
}
