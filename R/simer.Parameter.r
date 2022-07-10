# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Annotation parameters generator
#' 
#' Generate parameters for annotation data simulation.
#' 
#' Build date: Feb 24, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for map simulation.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$map$pop.map}{the map data with annotation information.}
#' \item{$map$qtn.num}{the QTN number for (each group in) each trait.}
#' \item{$map$qtn.model}{the genetic model of QTN such as 'A + D'.}
#' \item{$map$qtn.dist}{the QTN distribution containing 'norm', 'geom', 'gamma' or 'beta'.}
#' \item{$map$qtn.sd}{the standard deviations for normal distribution.}
#' \item{$map$qtn.prob}{the probability of success for geometric distribution.}
#' \item{$map$qtn.shape}{the shape parameter for gamma distribution.}
#' \item{$map$qtn.scale}{the scale parameter for gamma distribution.}
#' \item{$map$qtn.shape1}{the shape1 parameter for beta distribution.}
#' \item{$map$qtn.shape2}{the shape2 parameter for beta distribution.}
#' \item{$map$qtn.ncp}{the ncp parameter for beta distribution.}
#' \item{$map$qtn.spot}{the QTN distribution probability in each block.}
#' \item{$map$len.block}{the block length.}
#' \item{$map$maf}{the maf threshold, markers less than this threshold will be exclude.}
#' \item{$map$recom.spot}{whether to generate recombination events.}
#' \item{$map$range.hot}{the recombination times range in the hot spot.}
#' \item{$map$range.cold}{the recombination times range in the cold spot.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.annot(qtn.num = list(tr1 = 10))
#' str(SP)
param.annot <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$map)) {
    SP.map <- list(
      pop.map = NULL,
      qtn.model = "A",
      qtn.index = NULL,
      qtn.num = list(tr1 = 10),
      qtn.dist = list(tr1 = "norm"),
      qtn.sd = list(tr1 = NA),
      qtn.prob = list(tr1 = NA),
      qtn.shape = list(tr1 = NA),
      qtn.scale = list(tr1 = NA),
      qtn.shape1 = list(tr1 = NA),
      qtn.shape2 = list(tr1 = NA),
      qtn.ncp = list(tr1 = NA),
      qtn.spot = FALSE,
      len.block = 5e7,
      maf = NULL,
      recom.spot = FALSE,
      range.hot = 4:6,
      range.cold = 1:5
    )
    
    group1 <- c("pop.map", "qtn.model", "qtn.index")
    group2 <- c("qtn.num", "qtn.dist", "qtn.sd", "qtn.prob", "qtn.shape", "qtn.scale", "qtn.shape1", "qtn.shape2", "qtn.ncp")
    group3 <- c("qtn.spot", "len.block", "maf", "recom.spot", "range.hot", "range.cold")
    
    for (x in names(SP.tmp)) {
      if (x %in% names(SP.map)) {
        SP.map[[x]] <- SP.tmp[[x]]
      }
    }
    
    if (is.null(SP.map$qtn.index)) {
      nTrait <- length(SP.map$qtn.num)
    } else {
      nTrait <- length(SP.map$qtn.index)
    }
    
    if (nTrait > 1) {
      for (x in group2) {
        SP.map[[x]] <- rep(SP.map[[x]], nTrait)
        names(SP.map[[x]]) <- paste0("tr", 1:nTrait)
      }
    }
    
    for (i in 1:nTrait) {
      nGroup <- length(SP.map$qtn.num[[i]])
      if (length(SP.map$qtn.dist[[i]]) != nGroup) {
        SP.map$qtn.dist[[i]] <- rep(SP.map$qtn.dist[[i]], nGroup)
      }
      for (j in 1:nGroup) {
        if (SP.map$qtn.dist[[i]][j] == "norm") {
          SP.map$qtn.sd[[i]][j] <- 1
        } else if (SP.map$qtn.dist[[i]][j] == "geom") {
          SP.map$qtn.prob[[i]][j] <- 0.5
        } else if (SP.map$qtn.dist[[i]][j] == "gamma") {
          SP.map$qtn.shape[[i]][j] <- 1
          SP.map$qtn.scale[[i]][j] <- 1
        } else if (SP.map$qtn.dist[[i]][j] == "beta") {
          SP.map$qtn.shape1[[i]][j] <- 1
          SP.map$qtn.shape2[[i]][j] <- 1
          SP.map$qtn.ncp[[i]][j] <- 0
        } else {
          stop("QTN effect distribution should be 'norm', 'geom', 'gamma' or 'beta'!")
        }
      }
    }
    
    for (x in names(SP.map)) {
      if (all(is.na(unlist(SP.map[[x]])))) {
        SP.map[[x]] <- NULL
      }
    }
    
    if (!SP.map$qtn.spot) {
      SP.map$maf <- NULL
    }
    
    if (!SP.map$recom.spot) {
      SP.map$range.hot <- NULL
      SP.map$range.cold <- NULL
    }
    
  } else {
    SP.map <- SP$map
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.map)) {
      SP.map[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$map <- SP.map
  return(SP)
}

#' Genotype parameters generator
#' 
#' Generate parameters for genotype data simulation.
#' 
#' Build date: Feb 21, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for genotype simulation.
#'
#' @return 
#' the function returns a list containing
#' \describe{
#' \item{$geno$pop.geno}{the genotype data.}
#' \item{$geno$incols}{'1':one-column genotype represents an individual; '2': two-column genotype represents an individual.}
#' \item{$geno$pop.marker}{the number of markers.}
#' \item{$geno$pop.ind}{the number of individuals in the base population.}
#' \item{$geno$prob}{the genotype code probability.}
#' \item{$geno$rate.mut}{the mutation rate of the genotype data.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)
#' str(SP)
param.geno <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$geno)) {
    SP.geno <- list(
      pop.geno = NULL,
      incols = 1, 
      pop.marker = 1e4,
      pop.ind = 1e2,
      prob = NULL,
      rate.mut = 1e-8
    )
    
  } else {
    SP.geno <- SP$geno
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.geno)) {
      SP.geno[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$geno <- SP.geno
  return(SP)
}

#' Phenotype parameters generator
#' 
#' Generate parameters for phenotype data simulation.
#' 
#' Build date: Feb 21, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for phenotype simulation.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$pheno$pop}{the population information containing environmental factors and other effects.}
#' \item{$pheno$pop.ind}{the number of individuals in the base population.}
#' \item{$pheno$pop.rep}{the repeated times of repeated records.}
#' \item{$pheno$pop.rep.bal}{whether repeated records are balanced.}
#' \item{$pheno$pop.env}{a list of environmental factors setting.}
#' \item{$pheno$phe.model}{a list of genetic model of phenotype such as "T1 = A + E".}
#' \item{$pheno$phe.h2A}{a list of additive heritability.}
#' \item{$pheno$phe.h2D}{a list of dominant heritability.}
#' \item{$pheno$phe.h2GxG}{a list of GxG interaction heritability.}
#' \item{$pheno$phe.h2GxE}{a list of GxE interaction heritability.}
#' \item{$pheno$phe.h2PE}{a list of permanent environmental heritability.}
#' \item{$pheno$phe.var}{a list of phenotype variance.}
#' \item{$pheno$phe.corA}{the additive genetic correlation matrix.}
#' \item{$pheno$phe.corD}{the dominant genetic correlation matrix.}
#' \item{$pheno$phe.corGxG}{the GxG genetic correlation matrix.}
#' \item{$pheno$phe.corPE}{the permanent environmental correlation matrix.}
#' \item{$pheno$phe.corE}{the residual correlation matrix.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.pheno(phe.model = list(tr1 = "T1 = A + E"))
#' str(SP)
param.pheno <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$pheno)) {
    SP.pheno <- list(
      pop = NULL,
      pop.ind = 100,
      pop.rep = 1,
      pop.rep.bal = TRUE,
      pop.env = NULL,
      phe.model = list(tr1 = "T1 = A + E"),
      phe.h2A = list(tr1 = NA),
      phe.h2D = list(tr1 = NA),
      phe.h2GxG = list(tr1 = NULL),
      phe.h2GxE = list(tr1 = NULL),
      phe.h2PE = list(tr1 = NA),
      phe.var = list(tr1 = NA),
      phe.corA = NULL,
      phe.corD = NULL,
      phe.corGxG = NULL,
      phe.corPE = NULL,
      phe.corE = NULL
    )
    
    group1 <- c("pop", "pop.ind", "pop.rep", "pop.rep.bal", "pop.env")
    group2 <- c("phe.model", "phe.h2A", "phe.h2D", "phe.h2GxG", "phe.h2GxE", "phe.h2PE", "phe.var")
    group3 <- c("phe.corA", "phe.corA", "phe.corGxG", "phe.corGxE", "phe.corPE", "phe.corE")
    
    for (x in names(SP.tmp)) {
      if (x %in% names(SP.pheno)) {
        SP.pheno[[x]] <- SP.tmp[[x]]
      }
    }
    
    nTrait <- length(SP.pheno$phe.model)
    if (nTrait == 0) { nTrait <- 1 }
    
    if (nTrait > 1) {
      for (x in group2) {
        SP.pheno[[x]] <- rep(SP.pheno[[x]], nTrait)
        names(SP.pheno[[x]]) <- paste0("tr", 1:nTrait)
      }
    }
    
    for (i in 1:nTrait) {
      model.split <- unlist(strsplit(SP.pheno$phe.model[[i]], split = "\\s*\\=\\s*"))
      eff.name <- unlist(strsplit(model.split[2], split = "\\s*\\+\\s*"))
      eff.name <- unique(eff.name)
      if (SP.pheno$pop.rep > 1) {
        SP.pheno$phe.h2PE[[i]] <- 0.1
      }
      for (j in 1:length(eff.name)) {
        if (eff.name[j] == "A") {
          SP.pheno$phe.h2A[[i]] <- 0.3
        } else if (eff.name[j] == "D") {
          SP.pheno$phe.h2D[[i]] <- 0.1
        } else if (grepl(pattern = ":", x = eff.name[j])) {
          eff.split <- unlist(strsplit(eff.name[j], split = ":"))
          if (all(eff.split %in% c("A", "D"))) {
            GxG.tmp <- list(0.1)
            names(GxG.tmp) <- eff.name[j]
            SP.pheno$phe.h2GxG[[i]] <- c(SP.pheno$phe.h2GxG[[i]], GxG.tmp)
          } else {
            GxE.tmp <- list(0.1)
            names(GxE.tmp) <- eff.name[j]
            SP.pheno$phe.h2GxE[[i]] <- c(SP.pheno$phe.h2GxE[[i]], GxE.tmp)
          } # end if (all(eff.split %in% c("A", "D"))) {
        } # end if (eff.name[j] == "A") {
      } # end for (j in 1:length(eff.name)) {
    } # end for (i in 1:nTrait) {
    
    for (x in names(SP.pheno)) {
      if (all(is.na(unlist(SP.pheno[[x]])))) {
        SP.pheno[[x]] <- NULL
      }
    }
    
    if (nTrait > 1) {
      if (!is.null(SP.pheno$phe.h2A)) {
        SP.pheno$phe.corA <- diag(nTrait)
      }
      if (!is.null(SP.pheno$phe.h2D)) {
        SP.pheno$phe.corD <- diag(nTrait)
      }
      if (!is.null(SP.pheno$phe.h2GxG)) {
        SP.pheno$phe.corGxG <- rep(list(diag(nTrait)), length(SP.pheno$phe.h2GxG[[1]]))
        names(SP.pheno$phe.corGxG) <- names(SP.pheno$phe.h2GxG[[1]])
      }
      if (SP.pheno$pop.rep > 1) {
        SP.pheno$phe.corPE <- diag(nTrait)
      }
      SP.pheno$phe.corE <- diag(nTrait)
    }
    
  } else {
    SP.pheno <- SP$pheno
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.pheno)) {
      SP.pheno[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$pheno <- SP.pheno
  return(SP)
}

#' Selection parameters generator
#' 
#' Generate parameters for selection.
#' 
#' Build date: Apr 6, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for selection.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$sel$pop.sel}{the selected males and females.}
#' \item{$sel$ps}{if ps <= 1, fraction selected in selection of males and females; if ps > 1, ps is number of selected males and females.}
#' \item{$sel$decr}{whether the sort order is decreasing.}
#' \item{$sel$sel.crit}{the selection criteria, it can be 'TBV', 'TGV', and 'pheno'.}
#' \item{$sel$sel.single}{the single-trait selection method, it can be 'ind', 'fam', 'infam', and 'comb'.}
#' \item{$sel$sel.multi}{the multiple-trait selection method, it can be 'index', 'indcul', and 'tmd'.}
#' \item{$sel$index.wt}{the weight of each trait for multiple-trait selection.}
#' \item{$sel$index.tdm}{the index of tandem selection for multiple-trait selection.}
#' \item{$sel$goal.perc}{the percentage of goal more than the mean of scores of individuals.}
#' \item{$sel$pass.perc}{the percentage of expected excellent individuals.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.sel(sel.single = "comb")
#' str(SP)
param.sel <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$sel)) {
    SP.sel <- list(
      pop.sel = NULL,
      ps = c(0.8, 0.8),
      decr = TRUE,
      sel.crit = "pheno",
      sel.single = "comb",
      sel.multi = "index",
      index.wt = c(0.5, 0.5),
      index.tdm = 1,
      goal.perc = 0.1,
      pass.perc = 0.9
    )
    
  } else {
    SP.sel <- SP$sel
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.sel)) {
      SP.sel[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$sel <- SP.sel
  return(SP)
}

#' Reproduction parameters generator
#' 
#' Generate parameters for reproduction.
#' 
#' Build date: Apr 6, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for reproduction.
#'
#' @return 
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the male rate in the population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.reprod(reprod.way = "randmate")
#' str(SP)
param.reprod <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$reprod)) {
    SP.reprod <- list(
      pop.gen = 1,
      reprod.way = "randmate",
      sex.rate = 0.5,
      prog = 2
    )
    
    if (!is.null(SP.tmp$reprod.way)) {
      if (SP.tmp$reprod.way == "userped" & is.null(SP.tmp$userped)) {
        userped <- data.frame(
          index = 1:200,
          sir = c(rep(0, 100), rep(1:50, each = 2)),
          dam = c(rep(0, 100), rep(51:100, each = 2))
        )
        SP.reprod$userped <- userped
      }
    } 
    
  } else {
    SP.reprod <- SP$reprod
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.reprod)) {
      SP.reprod[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$reprod <- SP.reprod
  return(SP)
}

#' Global parameters generator
#' 
#' Generate parameters for global options.
#' 
#' Build date: Apr 16, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for global options.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$replication}{the replication times of simulation.}
#' \item{$seed.sim}{simulation random seed.}
#' \item{$out}{the prefix of output files.}
#' \item{$outpath}{the path of output files, Simer writes files only if outpath is not 'NULL'.}
#' \item{$out.format}{'numeric' or 'plink', the data format of output files.}
#' \item{$pop.gen}{the generations of simulated population.}
#' \item{$out.geno.gen}{the output generations of genotype data.}
#' \item{$out.pheno.gen}{the output generations of phenotype data.}
#' \item{$useAllGeno}{whether to use all genotype data to simulate phenotype.}
#' \item{$ncpus}{the number of threads used, if NULL, (logical core number - 1) is automatically used.}
#' \item{$verbose}{whether to print detail.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.global(out = "simer")
#' str(SP)
param.global <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  if (is.null(SP$global)) {
    SP.global <- list(
      replication = 1,
      seed.sim = runif(1, 0, 100),
      out = "simer", 
      outpath = NULL,
      out.format = "numeric",
      pop.gen = 1,
      out.geno.gen = 1,
      out.pheno.gen = 1,
      useAllGeno = FALSE,
      ncpus = 0,
      verbose = TRUE
    )
    
    if (!is.null(SP.tmp$pop.gen)) {
      SP.global$pop.gen <- SP.tmp$pop.gen
      SP.global$out.geno.gen <- 1:SP.tmp$pop.gen
      SP.global$out.pheno.gen <- 1:SP.tmp$pop.gen
    }
    
  } else {
    SP.global <- SP$global
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.global)) {
      SP.global[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$global <- SP.global
  return(SP)
}

#' Parameter generator
#' 
#' Generate parameters for Simer.
#' 
#' Build date: Apr 17, 2022
#' Last update: Jul 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ... one or more parameter(s) for simer.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$global}{a list of global parameters.}
#' \item{$map}{a list of marker information parameters.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' \item{$sel}{a list of selection parameters.}
#' \item{$reprod}{a list of reproduction parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' SP <- param.simer(out = "simer")
#' str(SP)
param.simer <- function(SP = NULL, ...) {
  
  SP <- param.global(SP = SP, ... = ...)
  SP <- param.annot(SP = SP, ... = ...)
  SP <- param.geno(SP = SP, ... = ...)
  SP <- param.pheno(SP = SP, ... = ...)
  SP <- param.sel(SP = SP, ... = ...)
  SP <- param.pheno(SP = SP, ... = ...)
  SP <- param.reprod(SP = SP, ... = ...)
  
  return(SP)
}
