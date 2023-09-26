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


#' Selection
#' 
#' Select individuals by combination of selection method and criterion.
#'
#' Build date: Sep 8, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param verbose whether to print detail.
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
#' \donttest{
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = list(tr1 = 10))
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "ind")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' }
selects <- function(SP = NULL, verbose = TRUE) {

### Start selection
  
  # selection parameters
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
  pop.gen <- SP$global$pop.gen - 1
  
  if (length(pop.gen) == 0) { pop.gen <- 0  }
  
  if (pop.gen == 0 || all(ps == 1)) {
    ind.ordered <- pop$index

  } else {
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
