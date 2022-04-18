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


#' Select individuals by combination of selection method and criterion
#'
#' Build date: Sep 8, 2018
#' Last update: Apr 4, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param verbose whether to print detail
#' 
#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' SP <- param.annot(qtn.num = 10)
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' SP <- annotation(SP)
#' SP <- genotype(SP)
#' SP <- phenotype(SP)
#' SP <- selects(SP)
#' str(SP$sel$pop.sel$gen1)
selects <- function(SP = NULL, verbose = TRUE) {

# Start selection
  
  # unfold selection parameters
  pop.sel <- SP$sel$pop.sel
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.total <- do.call(rbind, SP$pheno$pop)
  ps <- SP$sel$ps
  decr <- SP$sel$decr
  sel.crit <- SP$sel$sel.crit
  sel.single <- SP$sel$sel.single
  sel.multi <- SP$sel$sel.multi
  index.wt <- SP$sel$index.wt
  index.tdm <- SP$sel$index.tdm
  goal.perc <- SP$sel$goal.perc
  pass.perc <- SP$sel$pass.perc
  pop.ind <- length(pop$index)
  
  if (!(all(ps <= 1) | all(ps > 1))) {
    stop("Please input a correct ps!")
  }
  if (sel.crit == "pheno") {
    phe.name <- grep(pattern = "TBV", x = names(pop), value = TRUE)
    phe.name <- substr(phe.name, 1, nchar(phe.name) - 4)
  } else if (sel.crit == "TBV") {
    phe.name <- grep(pattern = "TBV", x = names(pop), value = TRUE)
  } else if (sel.crit == "TGV") {
    phe.name <- grep(pattern = "TGV", x = names(pop), value = TRUE)
  } else {
    stop("'sel.crit' should be 'pheno', 'TBV', 'TGV', 'pEBVs', 'gEBVs' or 'ssEBVs'!") 
  }
  phe.pos <- match(phe.name, names(pop))
  
  ### single trait selection ###
	if (length(phe.pos) == 1) {
	  num.infam <- tapply(rep(1, pop.ind), pop$fam, sum)
	  # calculate r by A matrix
	  if (sel.single == "comb") {
	    # calculate A matrix
	    A.cal <- function(s, d) {
	      n <- length(s)
	      A <- matrix(0, nrow = n, ncol = n)
	      for (x in 1:n) {
	        if (s[x] > 0 && d[x] > 0) {
	          A[x, x] <- 1 + 0.5 * A[s[x], d[x]];
	        } else {
	          A[x, x] <- 1;
	        }
	        if (x == n) next
	        for (y in (x+1):n) {
	          axy <- 0
	          if (s[y] > 0) axy <- axy + 0.5 * A[x, s[y]]
	          if (d[y] > 0) axy <- axy + 0.5 * A[x, d[y]]
	          A[x, y] <- axy
	          A[y, x] <- axy
	        }
	      }
	      return(A)
	    }
	    # calculate average family correlation coefficient
	    cal.r <- function(pop, pop.total) {
	      if (nrow(pop.total) == nrow(pop)) { return(0.01) }
	      sir <- pop.total$sir
	      sir <- as.numeric(sir[sir != "NA"])
	      dam <- pop.total$dam
	      dam <- as.numeric(dam[dam != "NA"])
	      A.mat <- A.cal(sir, dam)
	      cor.ani <- pop[!duplicated(pop$fam), ]$index
	      cor.r <- rep(0, length(cor.ani))
	      for (i in 1:length(cor.ani)) {
	        cor.r[i] <- A.mat[cor.ani[i], cor.ani[i]+1]
	      }
	      cor.r <- mean(cor.r)
	      return(cor.r)
	    }
	    cor.r <- cal.r(pop, pop.total)
	    # calculate combination selection method
	    cal.comb <- function(pop) {
	      pop.ind <- length(pop$index)
	      num.infam <- tapply(rep(1, pop.ind), pop$fam, sum)
	      pf <- rep(tapply(pop[, phe.pos], pop$fam, mean), num.infam)
	      cor.n <- num.infam[[1]]
	      
	      if (length(num.infam) == 1) {
	        cor.t <- 2
	      } else {
	        pop.tab <- summary(aov(pop[, phe.pos]~as.factor(pop$fam)))
	        MB <- pop.tab[[1]][1, 3]
	        MW <- pop.tab[[1]][2, 3]
	        if (is.na(MW)) MW <- 0
	        cor.t <- (MB - MW) / (MB + (cor.n - 1) * MW)
	      }
	      
	      I <- pop[, phe.pos] + ((cor.r-cor.t)*cor.n / (1-cor.r) / (1+(cor.n-1)*cor.t)) * pf
	      ind.score <- cbind(pop$index, I)
	      ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing=decr), ]
	      ind.ordered <- as.vector(ind.score.ordered[, 1])
	      return(ind.ordered)
	    }
	    ind.ordered <- cal.comb(pop)
	    
	  } else {
	    if (sel.single == "ind") {
	      bf <- rep(1, pop.ind)
	      bw <- rep(1, pop.ind)
	      
	    } else if (sel.single == "fam") {
	      bf <- rep(1, pop.ind)
	      bw <- rep(0, pop.ind)
	      
	    } else if (sel.single == "infam") {
	      bf <- rep(0, pop.ind)
	      bw <- rep(1, pop.ind)
	    }
	    bfw <- cbind(bf, bw)
	    
	    # calculate pf and pw
	    cal.pfw <- function(pop) {
	      pop.ind <- length(pop$index)
	      num.infam <- tapply(rep(1, pop.ind), pop$fam, sum)
	      pf <- rep(tapply(pop[, phe.pos], pop$fam, mean), num.infam)
	      pw <- pop[, phe.pos] - pf
	      pfw <- cbind(pf, pw)
	      return(pfw)
	    }
	    
	    pfw <- cal.pfw(pop)
	    I <- apply(bfw * pfw, 1, sum)
	    ind.score <- cbind(pop$index, I)
	    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing=decr), ]
	    ind.ordered <- ind.score.ordered[, 1]
	  }
	  
	### multiple trait selection ###
	} else {
	  pheno <- pop[, phe.pos]
	  
	  goal <- (1 + goal.perc) * apply(pheno, 2, mean)
	  
	  if (sel.multi == "index") {
	    # calculate the weight of index selection
	    cal.idx <- function(pheno, index.wt) {
	      if(ncol(pheno) != length(index.wt)) {
	        stop("Column of phenotype should equal length of weight of index")
	      }
	      # phenotype covariance matrix
	      P <- var(pheno)
	      iP <- try(solve(P), silent = TRUE)
	      if (inherits(iP, "try-error")) {
	        iP <- ginv(P)
	      }
	      TBV <- pop[, grep(pattern = "TBV", x = names(pop))]
	      # BV covariance matrix
	      A <- var(TBV)
	      b <- iP %*% A %*% index.wt
	      b <- as.vector(b)
	    } 
	    b <- cal.idx(pheno = pheno, index.wt = index.wt)
	    
	    ind.score <- cbind(pop$index, apply(pheno * b, 1, sum))
	    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing = decr), ]
	    ind.ordered <- ind.score.ordered[, 1]
	    attr(ind.ordered, "names") <- NULL
	    
	  } else if (sel.multi == "indcul") {
	    flag.prior <- apply(pheno, 1, function(v) {
	      return(all(v > goal))
	    })
	    ind.ordered <- c(pop$index[flag.prior], pop$index[!flag.prior])
	    
	  } else if (sel.multi == "tdm") {
	    ind.score <- cbind(pop$index, pheno[, index.tdm])
	    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing = decr), ]
	    ind.ordered <- ind.score.ordered[, 1]
	    if (index.tdm == ncol(pheno)) {
	      logging.log(" All phenotype have selected by tandem method.\n", verbose = verbose)
	      SP$sel$index.wt <- 1
	    }
	    if (pheno[nrow(pheno) * pass.perc, index.tdm] >= goal[index.tdm]) {
	      SP$sel$index.wt <- index.tdm + 1
	    }
	    
	  } else {
	    stop("sel.multi should be index, indcul or tdm")
	  }
	}
  
  # get selected sires and dams
  ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
  ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
  nSir <- ifelse(all(ps <= 1), round(length(ind.sir) * ps[1]), ps[1])
  nDam <- ifelse(all(ps <= 1), round(length(ind.dam) * ps[2]), ps[2])
  sir.ordered <- intersect(ind.ordered, ind.sir)
  dam.ordered <- intersect(ind.ordered, ind.dam)
  if (length(sir.ordered) == 0 | nSir == 0) {
    sir.stay <- NULL
  } else {
    if (nSir <= length(ind.sir)) {
      sir.stay <- sir.ordered[1:nSir]
    } else {
      sir.stay <- c(sir.ordered, sample(sir.ordered, nSir-length(ind.sir), replace = TRUE))
    }
  }
  if (length(dam.ordered) == 0 | nDam == 0) {
    dam.stay <- NULL
  } else {
    if (nDam <= length(ind.dam)) {
      dam.stay <- dam.ordered[1:nDam]
    } else {
      dam.stay <- c(dam.ordered, sample(dam.ordered, nDam-length(ind.dam), replace = TRUE))
    }
  }
  pop.sel <- list(sir = sir.stay, dam = dam.stay)
  
  SP$sel$pop.sel[[ifelse(is.null(SP$sel$pop.sel), 1, length(SP$sel$pop.sel) + 1)]] <- pop.sel
  names(SP$sel$pop.sel)[length(SP$sel$pop.sel)] <- paste0("gen", length(SP$sel$pop.sel))
	return(SP)
}

#' Get core population ID and genotype
#'
#' Build date: May 2, 2020
#' Last update: May 2, 2020
#'
#' @author Dong Yin
#'
#' @param pop.sel ID of the sires and the dams passing selection
#' @param core.stay ID the core sires and dams
#' @param refresh refresh ratio of core population of sires and dams, only used in ps > 1
#' @param keep.max.gen if ps <= 1, fraction selected in selection of males and females; if ps > 1, ps is number of selected males and females
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param pop.total total population information
#' @param pop.geno.curr genotype matrix of current population
#' @param pop.geno.core genotype matrix of core population
#'
#' @return core population list
#' @export
#'
#' @examples
#' # no example for now
sel.core <- function(pop.sel = NULL, core.stay = NULL, refresh = rep(0.6, 2), keep.max.gen = rep(3, 2), incols = 2, pop.total, pop.geno.curr, pop.geno.core) {
  if (any(refresh < 0.5)) 
    stop("Refresh ratio of the core population should be not less than 0.5!")
  if (length(refresh) != 2)
    stop("refresh should be set only for sires and dams!")
  if (length(keep.max.gen) != 2)
    stop("keep.max.gen should be set only for sires and dams!")
  refresh[keep.max.gen == 1] <- 1
  
  gen <- pop.total$gen[nrow(pop.total)]
  index.curr <- pop.total[pop.total$gen == gen, ]$index
  gen.core.sir <- pop.total[core.stay$sir, ]$gen
  gen.core.dam <- pop.total[core.stay$dam, ]$gen
  f.gen.sir <- gen.core.sir <= gen - keep.max.gen[1]
  f.gen.dam <- gen.core.dam <= gen - keep.max.gen[2]
  sf.gen.sir <- sum(f.gen.sir)
  sf.gen.dam <- sum(f.gen.dam)
  
  num.refresh.sir <- round(length(core.stay$sir) * refresh[1])
  num.refresh.dam <- round(length(core.stay$dam) * refresh[2])
  if (num.refresh.sir < sf.gen.sir) 
    stop("Refresh sires should be not less than old sires!")
  if (num.refresh.dam < sf.gen.dam) 
    stop("Refresh dams should be not less than old dams!")
  
  if (num.refresh.sir > sf.gen.sir) {
    id.add.sir <- sample(which(!f.gen.sir), size = num.refresh.sir-sf.gen.sir)
    f.gen.sir[id.add.sir] <- TRUE
  }
  if (num.refresh.dam > sf.gen.dam) {
    id.add.dam <- sample(which(!f.gen.dam), size = num.refresh.dam-sf.gen.dam)
    f.gen.dam[id.add.dam] <- TRUE
  }
  
  sir.stay <- pop.sel$sir[1:num.refresh.sir]
  dam.stay <- pop.sel$dam[1:num.refresh.dam]
  if (incols == 2) {
    gmt.dam.curr <- match(c(sir.stay, dam.stay), index.curr)
    gmt.dam.core <- which(c(f.gen.sir, f.gen.dam))
    
    gmt.comb.curr <- getgmt(gmt.dam.curr, incols = incols)
    gmt.comb.core <- getgmt(gmt.dam.core, incols = incols)

    pop.geno.core[, gmt.comb.core] <- pop.geno.curr[, gmt.comb.curr]
  } else {
    gmt.comb <- match(c(sir.stay, dam.stay), index.curr)
    pop.geno.core[, c(f.gen.sir, f.gen.dam)] <- pop.geno.curr[, gmt.comb]
  }
  core.stay$sir[f.gen.sir] <- sir.stay
  core.stay$dam[f.gen.dam] <- dam.stay
  
  core.list <- list(core.stay = core.stay, core.geno = pop.geno.core)
  return(core.list)
}
