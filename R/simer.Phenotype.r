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


#' Phenotype simulation
#' 
#' Generate single-trait or multiple-trait phenotype by mixed model.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param verbose whether to print detail.
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
#' @references Kao C and Zeng Z (2002) <https://www.genetics.org/content/160/3/1243.long>
#'
#' @examples
#' # Prepare environmental factor list
#' pop.env <- list(
#'   F1 = list( # fixed effect 1
#'     level = c("1", "2"),
#'     effect = list(tr1 = c(50, 30), tr2 = c(50, 30))
#'   ), 
#'   F2 = list( # fixed effect 2
#'     level = c("d1", "d2", "d3"),
#'     effect = list(tr1 = c(10, 20, 30), tr2 = c(10, 20, 30))
#'   ),
#'   C1 = list( # covariate 1
#'     level = c(70, 80, 90),
#'     intercept = list(tr1 = 1.5, tr2 = 1.5)
#'   ),
#'   R1 = list( # random effect 1
#'     level = c("l1", "l2", "l3"),
#'     ratio = list(tr1 = 0.1, tr2 = 0.1)
#'   )
#' )
#' 
#' # Generate genotype simulation parameters
#' SP <- param.annot(qtn.num = list(tr1 = 10, tr2 = 10),
#'                   qtn.model = "A + D + A:D")
#' # Generate annotation simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(
#'   SP = SP, 
#'   pop.ind = 100,
#'   pop.rep = 2, # 2 repeated record
#'   pop.rep.bal = TRUE, # balanced repeated record
#'   pop.env = pop.env,
#'   phe.var = list(tr1 = 100, tr2 = 100),
#'   phe.model = list(
#'     tr1 = "T1 = A + D + A:D + F1 + F2 + C1 + R1 + A:F1 + E",
#'     tr2 = "T2 = A + D + A:D + F1 + F2 + C1 + R1 + A:F1 + E"
#'   )
#' )
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
phenotype <- function(SP = NULL, verbose = TRUE) {

### Start phenotype simlulation
  
  # phenotype parameters
  useAllGeno <- SP$global$useAllGeno
  if (is.null(useAllGeno)) {
    useAllGeno <- FALSE
  }
  pop <- SP$pheno$pop
  pop.ind <- SP$pheno$pop.ind
  if (is.null(pop) & !is.null(pop.ind)) {
    pop <- generate.pop(pop.ind = pop.ind)
    SP$pheno$pop <- list(pop)
  }
  if (is.data.frame(SP$pheno$pop)) {
    pop <- SP$pheno$pop
    SP$pheno$pop <- list(pop)
  } else {
    if (useAllGeno) {
      pop <- do.call(rbind, lapply(SP$pheno$pop, function(p) {
        return(p[, 1:7])
      }))
    } else {
      pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    }
  }
  pop.rep <- SP$pheno$pop.rep
  pop.rep.bal <- SP$pheno$pop.rep.bal
  pop.map <- SP$map$pop.map
  pop.map.GxG <- SP$map$pop.map.GxG
  pop.geno <- SP$geno$pop.geno
  incols <- SP$geno$incols
  pop.env <- SP$pheno$pop.env
  phe.model <- SP$pheno$phe.model
  phe.h2A <- SP$pheno$phe.h2A
  phe.h2D <- SP$pheno$phe.h2D
  phe.h2GxG <- SP$pheno$phe.h2GxG
  phe.h2GxE <- SP$pheno$phe.h2GxE
  phe.h2PE <- SP$pheno$phe.h2PE
  phe.var <- SP$pheno$phe.var
  phe.corA <- SP$pheno$phe.corA
  phe.corD <- SP$pheno$phe.corD
  phe.corGxG <- SP$pheno$phe.corGxG
  phe.corPE <- SP$pheno$phe.corPE
  phe.corE <- SP$pheno$phe.corE
  
  if (is.null(SP)) {
    stop("'SP' should be specified!")
  }
  
  pop.ind <- nrow(pop)
  
  # add environmental factor
  if (!is.null(pop.env)) {
    env.name <- names(pop.env)
    env.list <- rep(list(NULL), length(env.name))
    for (i in 1:length(env.name)) {
      new.env <- pop.env[[env.name[i]]]
      old.env <- pop[[env.name[i]]]
      if (!(env.name[i] %in% names(pop))) {
        logging.log(" Add", env.name[i], "to population...\n", verbose = verbose)
        if (!is.null(new.env$intercept)) {
          if (!is.numeric(new.env$level)) {
            stop("The levels of covariates should be numeric!")
          }
        }
        env.order <- sample(1:length(new.env$level), pop.ind, replace = T)
        env.list[[i]] <- cbind(env.list[[i]], new.env$level[env.order])
        colnames(env.list[[i]])[ncol(env.list[[i]])] <- env.name[i]
      } else {
        if (any(sort(new.env$level) != sort(unique(old.env)))) {
          stop(paste0("The levels of ", env.name[[i]], " should be ", paste(sort(unique(old.env)), collapse = ', '), "!"))
        }
      }
    }
    pop <- cbind(pop, as.data.frame(do.call(cbind, env.list)))
  }
  
  nTrait <- length(phe.model)
  phe.eff <- rep(list(NULL), nTrait)
  phe.h2E <- NULL
  phe.name <- NULL
  qtn.pos.add <- qtn.pos.dom <- qtn.pos.GxG <- NULL
  pop.qtn.add <- pop.qtn.dom <- pop.qtn.GxG <- NULL
  pop.snp.add <- pop.snp.dom <- pop.snp.GxG <- NULL
  for (i in 1:nTrait) {
    model.split <- unlist(strsplit(phe.model[[i]], split = "\\s*\\=\\s*"))
    phe.name[[i]] <- model.split[1]
    eff.name <- unlist(strsplit(model.split[2], split = "\\s*\\+\\s*"))
    if (length(unique(eff.name)) != length(eff.name)) {
      eff.name <- unique(eff.name)
      warning("Repeated effect have been removed in the model!")
    }
    for (j in 1:length(eff.name)) {
      # additive effect
      if (eff.name[j] == "A") {
        qtn.eff <- pop.map[[paste0("QTN", i, "_A")]]
        if (is.null(qtn.eff)) {
          stop("No additive effect in map!")
        }
        qtn.pos.add[[i]] <- which(!is.na(qtn.eff))
        pop.snp.add[[i]] <- pop.map[qtn.pos.add[[i]], 1]
        qtn.eff <- qtn.eff[qtn.pos.add[[i]]]
        if (useAllGeno) {
          pop.qtn.add[[i]] <- do.call(cbind, lapply(1:length(pop.geno), function(k) {
            return(pop.geno[[k]][qtn.pos.add[[i]], ])
          }))
        } else {
          pop.qtn.add[[i]] <- pop.geno[[length(pop.geno)]][qtn.pos.add[[i]], ]
        }
        if (pop.ind != ncol(pop.qtn.add[[i]]) / incols) {
          stop("The individual number should be same in both 'pop' and 'pop.qtn.add'!")
        }
        if (ncol(pop.qtn.add[[i]]) == 2 * pop.ind) {
          pop.qtn.add[[i]] <- geno.cvt1(pop.qtn.add[[i]])
        }
        # (0, 1, 2) -> (-1, 0, 1)
        pop.qtn.add[[i]] <- pop.qtn.add[[i]] - 1
        phe.add <- crossprod(pop.qtn.add[[i]], qtn.eff)
        if (length(phe.var) < i) {
          phe.var[[i]] <- var(phe.add) / phe.h2A[[i]]
        } 
        scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2A[[i]] / var(phe.add)))
        SP$map$pop.map[[paste0("QTN", i, "_A")]] <- SP$map$pop.map[[paste0("QTN", i, "_A")]] * scale
        phe.add <- phe.add * scale
        phe.add <- scale(phe.add, scale = FALSE)
        phe.eff[[i]] <- cbind(phe.eff[[i]], phe.add)
        colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_A_eff")
      
      # dominant effect
      } else if (eff.name[j] == "D") {
        qtn.eff <- pop.map[[paste0("QTN", i, "_D")]]
        if (is.null(qtn.eff)) {
          stop("No dominant effect in map!")
        }
        qtn.pos.dom[[i]] <- which(!is.na(qtn.eff))
        pop.snp.dom[[i]] <- pop.map[qtn.pos.dom[[i]], 1]
        qtn.eff <- qtn.eff[qtn.pos.dom[[i]]]
        if (useAllGeno) {
          pop.qtn.dom[[i]] <- do.call(cbind, lapply(1:length(pop.geno), function(k) {
            return(pop.geno[[k]][qtn.pos.dom[[i]], ])
          }))
        } else {
          pop.qtn.dom[[i]] <- pop.geno[[length(pop.geno)]][qtn.pos.dom[[i]], ]
        }
        if (pop.ind != ncol(pop.qtn.dom[[i]]) / incols) {
          stop("The individual number should be same in both 'pop' and 'pop.qtn.dom'!")
        }
        if (ncol(pop.qtn.dom[[i]]) == 2 * pop.ind) {
          pop.qtn.dom[[i]] <- geno.cvt1(pop.qtn.dom[[i]])
        }
        # (0, 1, 2) -> (-0.5, 0.5, -0.5)
        pop.qtn.dom[[i]][pop.qtn.dom[[i]] == 2] <- 0
        pop.qtn.dom[[i]] <- pop.qtn.dom[[i]] - 0.5
        phe.dom <- crossprod(pop.qtn.dom[[i]], qtn.eff)
        scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2D[[i]] / var(phe.dom)))
        SP$map$pop.map[[paste0("QTN", i, "_D")]] <- SP$map$pop.map[[paste0("QTN", i, "_D")]] * scale
        phe.dom <- phe.dom * scale
        phe.dom <- scale(phe.dom, scale = FALSE)
        phe.eff[[i]] <- cbind(phe.eff[[i]], phe.dom)
        colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_D_eff")
      
      # fixed effect, covariate, and random effect
      } else if (eff.name[j] %in% names(pop)) {
        new.env <- pop.env[[eff.name[j]]]
        if (!is.null(new.env$effect)) {
          if (length(new.env$level) != length(new.env$effect[[i]])) {
            stop("The length of level and eff should be same!")
          }
        } else {
          if (!is.null(new.env$intercept)) {
            new.env$effect[[i]] <- new.env$level * new.env$intercept[[i]]
          } else {
            new.env$effect[[i]] <- rnorm(length(new.env$level))
          }
        }
        env.order <- match(pop[[eff.name[j]]], new.env$level)
        phe.env <- new.env$effect[[i]][env.order]
        eff.ratio <- pop.env[[eff.name[j]]]$ratio[[i]]
        if (!is.null(eff.ratio)) {
          scale <- as.numeric(sqrt(phe.var[[i]] * eff.ratio / var(phe.env)))
          phe.env <- phe.env * scale
          phe.env <- scale(phe.env, scale = FALSE)
        }
        phe.eff[[i]] <- cbind(phe.eff[[i]], phe.env)
        colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_", eff.name[j], "_eff")
        
      # interaction effect 
      } else if (grepl(pattern = ":", x = eff.name[j])) {
        eff.split <- unlist(strsplit(eff.name[j], split = ":"))
        # GxG effect
        if (all(eff.split %in% c("A", "D"))) {
          qtn.pos.GxG[[i]] <- grep(pattern = eff.name[j], x = pop.map.GxG[, 1])
          qtn.pos.GxG[[i]] <- qtn.pos.GxG[[i]][!is.na(pop.map.GxG[, i + 1])]
          pop.qtn.GxG[[i]] <- matrix(1, length(qtn.pos.GxG[[i]]), pop.ind)
          for (pos in qtn.pos.GxG[[i]]) {
            pop.snp.GxG <- unlist(strsplit(pop.map.GxG[pos, 1], split = "_"))[1]
            pop.snp.GxG <- unlist(strsplit(pop.snp.GxG, split = "-"))
            qtn.GxG <- rep(1, pop.ind)
            for (k in 1:length(pop.snp.GxG)) {
              if (eff.split[k] == "A") {
                pop.qtn.GxG[[i]][match(pos, qtn.pos.GxG[[i]]), ] <- pop.qtn.GxG[[i]][match(pos, qtn.pos.GxG[[i]]), ] * pop.qtn.add[[i]][match(pop.snp.GxG[k], pop.snp.add[[i]]), ]
              } else if (eff.split[k] == "D") {
                pop.qtn.GxG[[i]][match(pos, qtn.pos.GxG[[i]]), ] <- pop.qtn.GxG[[i]][match(pos, qtn.pos.GxG[[i]]), ] * pop.qtn.dom[[i]][match(pop.snp.GxG[k], pop.snp.dom[[i]]), ]
              } else {
                stop("Only 'A' or 'D' in the GxG!")
              }
            }
          }
          phe.GxG <- crossprod(pop.qtn.GxG[[i]], pop.map.GxG[qtn.pos.GxG[[i]], i + 1])
          scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2GxG[[i]][[eff.name[j]]] / var(phe.GxG)))
          SP$map$pop.map.GxG[qtn.pos.GxG[[i]], i + 1] <- SP$map$pop.map.GxG[qtn.pos.GxG[[i]], i + 1] * scale
          phe.GxG <- phe.GxG * scale
          phe.GxG <- scale(phe.GxG, scale = FALSE)
          phe.eff[[i]] <- cbind(phe.eff[[i]], phe.GxG)
          colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_", gsub(pattern = ":", replacement = "x", x = eff.name[j]), "_eff")
        # GxE effect
        } else {
          eff.GxE <- paste0(phe.name[[i]], "_", eff.split, "_eff")
          if (all(eff.GxE %in% colnames(phe.eff[[i]]))) {
            phe.GxE <- rep(1, pop.ind)
            for (k in 1:length(eff.GxE)) {
              if (eff.split[k] == "A") {
                phe.GxE <- phe.GxE * scale(crossprod(pop.qtn.add[[i]], rnorm(nrow(pop.qtn.add[[i]]))))
              } else if (eff.split[k] == "D") {
                phe.GxE <- phe.GxE * scale(crossprod(pop.qtn.dom[[i]], rnorm(nrow(pop.qtn.dom[[i]]))))
              } else {
                phe.GxE <- phe.GxE * scale(phe.eff[[i]][, match(eff.GxE[k], colnames(phe.eff[[i]]))])
              }
            }
            scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2GxE[[i]][[eff.name[j]]] / var(phe.GxE)))
            phe.GxE <- phe.GxE * scale
            phe.GxE <- scale(phe.GxE, scale = FALSE)
            phe.eff[[i]] <- cbind(phe.eff[[i]], phe.GxE)
            colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_", gsub(pattern = ":", replacement = "x", x = eff.name[j]), "_eff")
            
          } else {
            stop(paste0(setdiff(eff.GxE, colnames(phe.eff[[i]])), " should be in the population!"))
          }
        }
        
      # residual effect
      } else if (eff.name[j] == "E") {
        phe.res <- rnorm(pop.ind)
        sum.ratio <- sum(unlist(lapply(which(names(pop.env) %in% eff.name), function(j) {return(pop.env[[j]]$ratio[[i]])})))
        phe.h2E[[i]] <- 1 - sum(phe.h2A[[i]], phe.h2D[[i]], unlist(phe.h2GxG[[i]]), unlist(phe.h2GxE[[i]]), phe.h2PE[[i]], sum.ratio)
        if (phe.h2E[[i]] < 0) {
          stop("Residual variance cannot be less than 0!")
        }
        scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2E[[i]] / var(phe.res)))
        phe.res <- phe.res * scale
        phe.res <- scale(phe.res, scale = FALSE)
        phe.eff[[i]] <- cbind(phe.eff[[i]], phe.res)
        colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_E_eff")
        
      } else {
        stop("There may be error in the model!")
      }
      
    }
    # repeated record model
    if (pop.rep > 1) {
      phe.PE <- rnorm(pop.ind)
      scale <- as.numeric(sqrt(phe.var[[i]] * phe.h2PE[[i]] / var(phe.PE)))
      phe.PE <- phe.PE * scale
      phe.PE <- scale(phe.PE, scale = FALSE)
      phe.eff[[i]] <- cbind(phe.eff[[i]], phe.PE)
      colnames(phe.eff[[i]])[ncol(phe.eff[[i]])] <- paste0(phe.name[[i]], "_PE_eff")
    }
    # check data quality
    if (options("simer.show.warning") == TRUE) {
      eff.cor <- cor(phe.eff[[i]])
      diag(eff.cor) <- 0
      if (sum(eff.cor > 0.5)) {
        warning("There are hign-correlations between fixed effects or fixed effects and random effects, and it will reduce the accuracy of effects in the simulation!")
      }
    }
  }
  
  # effect combination
  phe.eff <- do.call(cbind, phe.eff)
  phe.eff <- as.data.frame(phe.eff)
  # update phenotype variance
  names(phe.var) <- paste0("tr", 1:nTrait)
  SP$pheno$phe.var <- phe.var
  
  # build correlation
  if (nTrait > 1) {
    if (!is.null(phe.corA)) {
      phe.add <- phe.eff[, grep(pattern = "_A_", x = names(phe.eff))]
      phe.varA <- unlist(phe.h2A) * unlist(phe.var)
      Sigma <- diag(sqrt(phe.varA)) %*% phe.corA %*% diag(sqrt(phe.varA))
      phe.add <- build.cov(df = phe.add, Sigma = Sigma)
      phe.eff[, grep(pattern = "_A_", x = names(phe.eff))] <- phe.add
      for (i in 1:nTrait) {
        SP$map$pop.map[[paste0("QTN", i, "_A")]][qtn.pos.add[[i]]] <- c(crossprod(ginv(pop.qtn.add[[i]]), phe.add[, i]))
      }
    }
    if (!is.null(phe.corD)) {
      phe.dom <- phe.eff[, grep(pattern = "_D_", x = names(phe.eff))]
      phe.varD <- unlist(phe.h2D) * unlist(phe.var)
      Sigma <- diag(sqrt(phe.varD)) %*% phe.corD %*% diag(sqrt(phe.varD))
      phe.dom <- build.cov(df = phe.dom, Sigma = Sigma)
      phe.eff[, grep(pattern = "_D_", x = names(phe.eff))] <- phe.dom
      for (i in 1:nTrait) {
        SP$map$pop.map[[paste0("QTN", i, "_D")]][qtn.pos.dom[[i]]] <- c(crossprod(ginv(pop.qtn.dom[[i]]), phe.dom[, i]))
      }
    }
    if (!is.null(phe.corGxG)) {
      eff.name <- names(phe.corGxG)
      eff.name <- gsub(pattern = ":", replacement = "x", x = eff.name)
      for (j in 1:length(eff.name)) {
        phe.GxG <- phe.eff[, grep(pattern = paste0("_", eff.name[j], "_"), x = names(phe.eff))]
        phe.varGxG <- unlist(phe.h2GxG) * unlist(phe.var)
        Sigma <- diag(sqrt(phe.varGxG)) %*% phe.corGxG[[j]] %*% diag(sqrt(phe.varGxG))
        phe.GxG <- build.cov(df = phe.GxG, Sigma = Sigma)
        phe.eff[, grep(pattern = paste0("_", eff.name[j], "_"), x = names(phe.eff))] <- phe.GxG
        for (i in 1:nTrait) {
          SP$map$pop.map.GxG[qtn.pos.GxG[[i]], i + 1] <- c(crossprod(ginv(pop.qtn.GxG[[i]]), phe.GxG[, i]))
        }
      }
    }
    if (!is.null(phe.corPE)) {
      phe.PE <- phe.eff[, grep(pattern = "_PE_", x = names(phe.eff))]
      phe.varPE <- unlist(phe.h2PE) * unlist(phe.var)
      Sigma <- diag(sqrt(phe.varPE)) %*% phe.corPE %*% diag(sqrt(phe.varPE))
      phe.PE <- build.cov(df = phe.PE, Sigma = Sigma)
      phe.eff[, grep(pattern = "_PE_", x = names(phe.eff))] <- phe.PE
    }
    if (!is.null(phe.corE)) {
      phe.res <- phe.eff[, grep(pattern = "_E_", x = names(phe.eff))]
      phe.varE <- unlist(phe.h2E) * unlist(phe.var)
      Sigma <- diag(sqrt(phe.varE)) %*% phe.corE %*% diag(sqrt(phe.varE))
      phe.res <- build.cov(df = phe.res, Sigma = Sigma)
      phe.eff[, grep(pattern = "_E_", x = names(phe.eff))] <- phe.res
    }
  }
  
  # if it is a repeated record trait
  if (pop.rep > 1) {
    phe.eff <- do.call(rbind, rep(list(phe.eff), pop.rep))
    phe.PE <- phe.eff[, grep(pattern = "_PE_", x = names(phe.eff)), drop = FALSE]
    if (!pop.rep.bal) {
      for (i in 1:nTrait) {
        phe.PE[sample(1:nrow(phe.PE), nrow(phe.PE) * 0.3), i] <- NA
      }
    }
    phe.eff[, grep(pattern = "_PE_", x = names(phe.eff))] <- phe.PE
    phe.res <- phe.eff[, grep(pattern = "_E_", x = names(phe.eff)), drop = FALSE]
    phe.varE <- unlist(phe.h2E) * unlist(phe.var)
    phe.res[] <- rnorm(nrow(phe.res) * ncol(phe.res))
    phe.res <- scale(phe.res)
    phe.res <- sweep(phe.res, 2, sqrt(phe.varE), "*")
    if (nTrait > 1) {
      Sigma <- diag(sqrt(phe.varE)) %*% phe.corE %*% diag(sqrt(phe.varE))
      phe.res <- build.cov(df = phe.res, Sigma = Sigma)
    }
    phe.eff[, grep(pattern = "_E_", x = names(phe.eff))] <- phe.res
  }
  
  # effect summary
  phe <- as.data.frame(do.call(cbind, lapply(1:nTrait, function(i) {
    return(rowSums(phe.eff[, grep(pattern = phe.name[[i]], x = names(phe.eff))]))
  })))
  names(phe) <- phe.name
  TBV <- as.data.frame(do.call(cbind, lapply(1:nTrait, function(i) {
    return(rowSums(phe.eff[, grep(pattern = paste0(phe.name[[i]], "_A_eff"), x = names(phe.eff)), drop = FALSE]))
  })))
  names(TBV) <- paste0(phe.name, "_TBV")
  TGV <- as.data.frame(do.call(cbind, lapply(1:nTrait, function(i) {
    return(rowSums(phe.eff[, grep(pattern = paste0(phe.name[[i]], c("_A_eff", "_D_eff"), collapse = "|"), x = names(phe.eff)), drop = FALSE]))
  })))
  names(TGV) <- paste0(phe.name, "_TGV")
  
  pop <- do.call(rbind, rep(list(pop), pop.rep))
  pop <- cbind(pop, phe, TBV, TGV, phe.eff)
  
  if (useAllGeno) {
    for (i in 1:length(SP$pheno$pop)) {
      SP$pheno$pop[[i]] <- pop[pop$gen == i, ]
    }
  } else {
    SP$pheno$pop[[length(SP$pheno$pop)]] <- pop
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  }
  
  return(SP)
}

#' Population generator
#' 
#' Generate population according to the number of individuals.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param pop.ind the number of the individuals in a population.
#' @param from initial index of the population.
#' @param ratio sex ratio of males in a population.
#' @param gen generation ID of the population.
#'
#' @return a data frame of population information.
#' 
#' @export
#'
#' @examples
#' pop <- generate.pop(pop.ind = 100)
#' head(pop)
generate.pop <- function(pop.ind = 100, from = 1, ratio = 0.5, gen = 1) {
  
  pop <- data.frame(
    index = seq(from = from, length.out = pop.ind),
    gen   = rep(gen, pop.ind),
    fam   = seq(from = from, length.out = pop.ind),
    infam = seq(from = from, length.out = pop.ind),
    sir   = rep(0, pop.ind),
    dam   = rep(0, pop.ind),
    sex   = rep(1:2, c(floor(ratio*pop.ind), (pop.ind - floor(ratio*pop.ind))))
  )
  
  return(pop)
}

#' Correlation building
#' 
#' To bulid correlation of variables.
#'
#' Build date: Oct 10, 2019
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin and R
#'
#' @param df a data frame needing building correlation.
#' @param mu means of the variables.
#' @param Sigma covariance matrix of variables.
#' @param tol tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma.
#'
#' @return a data frame with expected correlation
#' 
#' @export
#' 
#' @references B. D. Ripley (1987) Stochastic Simulation. Wiley. Page 98
#'
#' @examples
#' df <- data.frame(tr1 = rnorm(100), tr2 = rnorm(100))
#' df.cov <- build.cov(df)
#' var(df.cov)
build.cov <- function(df = NULL, mu = rep(0, nrow(Sigma)), Sigma = diag(2), tol = 1e-06) {
  if (!is.data.frame(df)) {
    df.nm <- paste0("tr", 1:ncol(df))
  } else {
    df.nm <- names(df)
  }
  
  # get zero-var index
  df.var <- apply(df, 2, var)
  idx <- which(df.var == 0)
  df.t <- df[, idx]
  df[, idx] <- rnorm(nrow(df))
  
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  df <- scale(df, center = TRUE, scale = FALSE)
  df <- df %*% svd(df, nu = 0)$v
  df <- scale(df, center = FALSE, scale = TRUE)
  
  df <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(df)
  df <- t(df)
  df <- as.data.frame(df)
  names(df) <- df.nm
  df[, idx] <- df.t
  
  return(df)
}
