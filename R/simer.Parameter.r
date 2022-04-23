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


#' Generate parameters for annotation data simulation
#' 
#' Build date: Feb 24, 2022
#' Last update: Mar 23, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for map simulation
#'
#' @return a list of map simulation parameters 
#' @export
#'
#' @examples
#' SP <- param.annot(qtn.num = 10)
#' str(SP)
param.annot <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  SP.map <- list(
    pop.map = NULL,
    qtn.num = 10,
    qtn.model = "A",
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
  
  group1 <- c("pop.map", "qtn.num", "qtn.model")
  group2 <- c("qtn.dist", "qtn.sd", "qtn.prob", "qtn.shape", "qtn.scale", "qtn.shape1", "qtn.shape2", "qtn.ncp")
  group3 <- c("qtn.spot", "len.block", "maf", "recom.spot", "range.hot", "range.cold")
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.map)) {
      SP.map[[x]] <- SP.tmp[[x]]
    }
  }
  
  nTrait <- 1
  if (is.matrix(SP.map$qtn.num)) {
    nTrait <- nrow(SP.map$qtn.num)
  }
  
  if (nTrait > 1) {
    for (x in group2) {
      SP.map[[x]] <- rep(SP.map[[x]], nTrait)
      names(SP.map[[x]]) <- paste0("tr", 1:nTrait)
    }
  }
  
  for (i in 1:nTrait) {
    if (SP.map$qtn.dist[[i]] == "norm") {
      SP.map$qtn.sd[[i]] <- 1
    } else if (SP.map$qtn.dist[[i]] == "geom") {
      SP.map$qtn.prob[[i]] <- 0.5
    } else if (SP.map$qtn.dist[[i]] == "gamma") {
      SP.map$qtn.shape[[i]] <- 1
      SP.map$qtn.scale[[i]] <- 1
    } else if (SP.map$qtn.dist[[i]] == "beta") {
      SP.map$qtn.shape1[[i]] <- 1
      SP.map$qtn.shape2[[i]] <- 1
      SP.map$qtn.ncp[[i]] <- 0
    } else {
      stop("QTN effect distribution should be 'norm', 'geom', 'gamma' or 'beta'!")
    }
  }
  
  if (nTrait == 1 & length(SP.tmp$qtn.num) > 0) {
    for (x in group2) {
      SP.map[[x]][[1]] <- rep(SP.map[[x]][[1]], length(SP.tmp$qtn.num))
    }
  }
  
  for (x in names(SP.map)) {
    if (all(is.na(unlist(SP.map[[x]])))) {
      SP.map[[x]] <- NULL
    }
  }
  
  if (!SP.map$qtn.spot) {
    SP.map$len.block <- NULL
    SP.map$maf <- NULL
  }
  
  if (!SP.map$recom.spot) {
    SP.map$range.hot <- NULL
    SP.map$range.cold <- NULL
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.map)) {
      SP.map[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$map <- SP.map
  return(SP)
}

#' Generate parameters for genotype data simulation
#' 
#' Build date: Feb 21, 2022
#' Last update: Feb 21, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for genotype simulation
#'
#' @return a list of genotype simulation parameters 
#' @export
#'
#' @examples
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)
#' str(SP)
param.geno <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  SP.geno <- list(
    pop.geno = NULL,
    incols = 1, 
    pop.marker = 1e4,
    pop.ind = 1e2,
    prob = NULL,
    rate.mut = 1e-8
  )
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.geno)) {
      SP.geno[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$geno <- SP.geno
  return(SP)
}

#' Generate parameters for phenotype data simulation
#' 
#' Build date: Feb 21, 2022
#' Last update: Mar 23, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for phenotype simulation
#'
#' @return a list of phenotype simulation parameters 
#' @export
#'
#' @examples
#' SP <- param.pheno(phe.model = list(tr1 = "T1 = A + E"))
#' str(SP)
param.pheno <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
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
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.pheno)) {
      SP.pheno[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$pheno <- SP.pheno
  return(SP)
}

#' Generate parameters for selection
#' 
#' Build date: Apr 6, 2022
#' Last update: Apr 6, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for selection
#'
#' @return a list of selection parameters 
#' @export
#'
#' @examples
#' SP <- param.sel(sel.single = "comb")
#' str(SP)
param.sel <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
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
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.sel)) {
      SP.sel[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$sel <- SP.sel
  return(SP)
}

#' Generate parameters for reproduction
#' 
#' Build date: Apr 6, 2022
#' Last update: Apr 6, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for reproduction
#'
#' @return a list of selection parameters 
#' @export
#'
#' @examples
#' SP <- param.reprod(reprod.way = "randmate")
#' str(SP)
param.reprod <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  SP.reprod <- list(
    pop.gen = 2,
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
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.reprod)) {
      SP.reprod[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$reprod <- SP.reprod
  return(SP)
}

#' Generate parameters for global options
#' 
#' Build date: Apr 16, 2022
#' Last update: Apr 16, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for global options
#'
#' @return a list of global options
#' @export
#'
#' @examples
#' SP <- param.global(out = "simer")
#' str(SP)
param.global <- function(SP = NULL, ...) {
  
  SP.tmp <- list(...)
  
  SP.global <- list(
    replication = 1,
    seed.sim = runif(1, 0, 100),
    out = "simer", 
    outpath = NULL,
    out.format = "numeric",
    pop.gen = 2,
    out.geno.gen = 1:2,
    out.pheno.gen = 1:2,
    useAllGeno = FALSE,
    ncpus = 0,
    verbose = TRUE
  )
  
  if (!is.null(SP.tmp$pop.gen)) {
    SP.global$pop.gen <- SP.tmp$pop.gen
    SP.global$out.geno.gen <- 1:SP.tmp$pop.gen
    SP.global$out.pheno.gen <- 1:SP.tmp$pop.gen
  }
  
  for (x in names(SP.tmp)) {
    if (x %in% names(SP.global)) {
      SP.global[[x]] <- SP.tmp[[x]]
    }
  }
  
  SP$global <- SP.global
  return(SP)
}

#' Generate parameters for simer
#' 
#' Build date: Apr 17, 2022
#' Last update: Apr 17, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param ... one or more parameter(s) for simer
#'
#' @return a list of simer parameters
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
