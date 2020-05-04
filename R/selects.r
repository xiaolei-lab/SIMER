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


#' Select individuals by combination of secletion method and criterion
#'
#' Build date: Sep 8, 2018
#' Last update: Oct 13, 2019
#'
#' @author Dong Yin
#'
#' @param pop population information of generation, family ID, within-family ID, individual ID, paternal ID, maternal ID, sex, and phenotype
#' @param decr whether to sorting with decreasing
#' @param sel.multi selection method of multiple traits with options: "tdm", "indcul" and "index"
#' @param index.wt economic weights of selection index method, its length should equals to the number of traits
#' @param index.tdm index represents which trait is being selected
#' @param goal.perc percentage of goal more than the mean of scores of individuals
#' @param pass.perc percentage of expected excellent individuals
#' @param sel.sing selection methods of single trait with the options: "ind", "fam", "infam", and "comb"
#' @param pop.total total population infarmation
#' @param pop.pheno list of all phenotype information
#' @param verbose whether to print detail

#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' 
#' effs <-
#'     cal.effs(pop.geno = pop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' 
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               FR = NULL, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' pop <- pop.pheno$pop
#' pop.pheno$pop <- NULL 
#' 
#' ind.ordered1 <-
#'     selects(pop = pop,
#'             decr = TRUE,
#'             sel.multi = "index",
#'             index.wt = c(0.5, 0.5),
#'             index.tdm = 1,
#'             goal.perc = 0.1,
#'             pass.perc = 0.9,
#'             sel.sing = "comb",
#'             pop.total = pop,
#'             pop.pheno = pop.pheno, 
#'             verbose = TRUE)
#'
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' effs <-
#'     cal.effs(pop.geno = pop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = TRUE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#'              
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               FR = NULL, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' pop <- pop.pheno$pop
#' pop.pheno$pop <- NULL
#' 
#' ind.ordered2 <-
#'     selects(pop = pop,
#'             decr = TRUE,
#'             sel.multi = "index",
#'             index.wt = c(0.5, 0.5),
#'             index.tdm = 1,
#'             goal.perc = 0.1,
#'             pass.perc = 0.9,
#'             sel.sing = "comb",
#'             pop.total = pop,
#'             pop.pheno = pop.pheno, 
#'             verbose = TRUE)
#' str(ind.ordered1)
#' str(ind.ordered2)
selects <-
    function(pop,
             decr = TRUE,
             sel.multi = "index",
             index.wt = c(0.5, 0.5),
             index.tdm = 1,
             goal.perc = 0.1,
             pass.perc = 0.9,
             sel.sing = "comb",
             pop.total = NULL,
             pop.pheno = NULL, 
             verbose = TRUE) {

# Start selection
  
  f1 <- grep(pattern = "TBV|TGV|pheno|ebv|u1", x = names(pop), value = TRUE)
  pheno <- subset(pop, select = f1)         
	if (ncol(pheno) == 1) {
	  # calculate r by A matrix
	  if (sel.sing == "comb") {
	    cor.r <- cal.r(pop.total, pop)
	  } else {
	    cor.r <- NULL
	  }
		ind.ordered <- cal.sing(pop, decr, sel.sing, cor.r)
		
	} else {
	  # select desired individuals
	  if (sel.multi == "indcul" || sel.multi == "tdm"){
	    goal <- (1 + goal.perc) * apply(pheno, 2, mean)
	  } else {
	    goal <- NULL
	  }
	  ind.ordered <- cal.multi(pop, decr, sel.multi, index.wt, index.tdm, goal, pass.perc, pop.pheno)
	}

	if (sel.multi == "tdm") {
    if (index.tdm == ncol(pheno)) {
      logging.log(" All phenotype have selected by tandem method.\n", verbose = verbose)
      index.tdm <- 1
    }
    if (pheno[nrow(pheno) * pass.perc, index.tdm] >= goal[index.tdm]) {
      index.tdm <- index.tdm + 1
    }
  }
	return(c(index.tdm, ind.ordered))
}

#' Single trait selection method
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye
#' @param decr whether to sorting with descreasing
#' @param sel.sing selection method of single trait with "ind", "fam", "infam", "comb"
#' @param cor.r average relationship of population
#'
#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' 
#' effs <-
#'     cal.effs(pop.geno = pop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' 
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               FR = NULL, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' pop <- pop.pheno$pop
#' pop.pheno$pop <- NULL
#' 
#' cor.r <- cal.r(pop, pop)
#' ind.ordered <- cal.sing(pop = pop, decr = TRUE, sel.sing = "comb", cor.r = cor.r)
#' str(ind.ordered)
cal.sing <- function(pop, decr, sel.sing, cor.r) {
  num.ind <- length(pop$index)
  num.infam <- tapply(rep(1, num.ind), pop$fam, sum)

  if (sel.sing == "ind") {
    bf <- rep(1, num.ind)
    bw <- rep(1, num.ind)

  } else if (sel.sing == "fam") {
    bf <- rep(1, num.ind)
    bw <- rep(0, num.ind)

  } else if (sel.sing == "infam") {
    bf <- rep(0, num.ind)
    bw <- rep(1, num.ind)

  } else if (sel.sing == "comb") {
    # calculate combination selection method
    cal.comb <- function(pop) {
      num.ind <- length(pop$index)
      num.infam <- tapply(rep(1, num.ind), pop$fam, sum)
      pf <- rep(tapply(pop$pheno, pop$fam, mean), num.infam)
      cor.n <- num.infam[[1]]

      if (length(num.infam) == 1) {
        cor.t <- 2
      } else {
        pop.tab <- summary(aov(pop$pheno~as.factor(pop$fam)))
        MB <- pop.tab[[1]][1, 3]
        MW <- pop.tab[[1]][2, 3]
        if (is.na(MW)) MW <- 0
        cor.t <- (MB - MW) / (MB + (cor.n - 1) * MW)
      }
      
      I <- pop$pheno + ((cor.r-cor.t)*cor.n / (1-cor.r) / (1+(cor.n-1)*cor.t)) * pf
      ind.score <- cbind(pop$index, I)
      ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing=decr), ]
      ind.ordered <- as.vector(ind.score.ordered[, 1])
      return(ind.ordered)
    }

    ind.ordered <- cal.comb(pop)
    return(ind.ordered)
  }
  bfw <- cbind(bf, bw)

  # calculate pf and pw
  cal.pfw <- function(pop) {
    num.ind <- length(pop$index)
    num.infam <- tapply(rep(1, num.ind), pop$fam, sum)
    pf <- rep(tapply(pop$pheno, pop$fam, mean), num.infam)
    pw <- pop$pheno - pf
    pfw <- cbind(pf, pw)
    return(pfw)
  }

  pfw <- cal.pfw(pop)
  I <- apply(bfw * pfw, 1, sum)
  ind.score <- cbind(pop$index, I)
  ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing=decr), ]
  ind.ordered <- ind.score.ordered[, 1]
  return(ind.ordered)
}

#' Multiple trait selection method
#'
#' Build date: Nov 14, 2018
#' Last update: Oct 13, 2019
#'
#' @author Dong Yin
#'
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye
#' @param decr whether to sorting with descreasing
#' @param sel.multi selection method of multiple traits with options: "tdm", "indcul" and "index"
#' @param index.wt economic weights of selection index method, its length should equals to the number of traits
#' @param index.tdm index represents which trait is being selected
#' @param goal goal of the trait
#' @param pass.perc percentage of expected excellent individuals
#' @param pop.pheno list of all phenotype information
#'
#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' 
#' effs <-
#'     cal.effs(pop.geno = pop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = TRUE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#'              
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               FR = NULL, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE,
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' pop <- pop.pheno$pop
#' pop.pheno$pop <- NULL 
#' 
#' f1 <- grep(pattern = "TBV|TGV|pheno|ebv", x = names(pop), value = TRUE)
#' pheno <- subset(pop, select = f1)                  
#' goal.perc <- 0.9
#' goal <- (1 + goal.perc) * apply(pheno, 2, mean)
#' 
#' ind.ordered <- cal.multi(pop = pop, decr = TRUE, sel.multi = "index",
#'     index.wt = c(0.5, 0.5), index.tdm = 1, goal = goal, 
#'     pass.perc = 0.9, pop.pheno = pop.pheno)
#' str(ind.ordered)
cal.multi <- function(pop, decr, sel.multi, index.wt, index.tdm, goal, pass.perc, pop.pheno) {
  num.ind <- length(pop$index)
  f1 <- grep(pattern = "TBV|TGV|pheno|ebv", x = names(pop), value = TRUE)
  pheno <- subset(pop, select = f1)    

  if (sel.multi == "index") {
    # calculate the weight of index selection
    b <- cal.idx(pop.pheno = pop.pheno, index.wt = index.wt)
    
    ind.score <- cbind(pop$index, apply(pheno * b, 1, sum))
    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing = decr), ]
    ind.ordered <- ind.score.ordered[, 1]
    attr(ind.ordered, "names") <- NULL

  } else if (sel.multi == "indcul") {
  	# goal <- (1 + goal.perc) * apply(pheno, 2, mean)
  	flag.prior <- apply(pheno, 1, function(v) {
  	  return(all(v > goal))
  	})
  	ind.ordered <- c(pop$index[flag.prior], pop$index[!flag.prior])

  } else if (sel.multi == "tdm") {
    ind.score <- cbind(pop$index, pheno[, index.tdm])
    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing = decr), ]
    ind.ordered <- ind.score.ordered[, 1]

  } else {
  	stop("sel.multi should be index, indcul or tdm")
  }
  return(ind.ordered)
}

#' Calculate A matrix
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param s indice of sires
#' @param d indice of dams
#'
#' @return additive kinship matrix
#' @export
#'
#' @examples
#' s <- c(0, 0, 0, 0, 1, 3, 3, 1, 5, 7, 5, 7, 1, 3, 5, 7)
#' d <- c(0, 0, 0, 0, 2, 4, 4, 2, 6, 8, 8, 6, 6, 8, 4, 8)
#' A.mat <- A.cal(s, d)
#' A.mat
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

#' Calculate average family correlation coefficient
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param pop total population information
#' @param pop.curr population information of current generation
#'
#' @return the average correlation of current population
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' pop.list <- reproduces(pop1 = basepop,
#'                        pop1.geno.id = basepop$index,
#'                        pop1.geno = basepop.geno,
#'                        ind.stay = list(sir=1:10, dam=11:100),
#'                        mtd.reprod = "randmate",
#'                        num.prog = 4,
#'                        ratio = 0.5)
#' pop.curr <- pop.list$pop
#' pop <- rbind(basepop, pop.curr)
#' cor.r <- cal.r(pop = pop, pop.curr = pop.curr)
#' cor.r
cal.r <- function(pop, pop.curr) {
  if (nrow(pop) == nrow(pop.curr)) { return(0.01) }
  sir <- pop$sir
  sir <- as.numeric(sir[sir != "NA"])
  dam <- pop$dam
  dam <- as.numeric(dam[dam != "NA"])
  A.mat <- A.cal(sir, dam)
  cor.ani <- pop.curr[!duplicated(pop.curr$fam), ]$index
  cor.r <- rep(0, length(cor.ani))
  for (i in 1:length(cor.ani)) {
    cor.r[i] <- A.mat[cor.ani[i], cor.ani[i]+1]
  }
  cor.r <- mean(cor.r)
  return(cor.r)
}

#' Calculate the weight of index selection
#'
#' Build date: Feb 8, 2020
#' Last update: Feb 8, 2020
#'
#' @author Dong Yin
#'
#' @param pop.pheno list of all phenotype information
#' @param index.wt economic weights of selection index method, its length should equals to the number of traits
#'
#' @return the weight of index selection
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' 
#' effs <-
#'     cal.effs(pop.geno = pop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = TRUE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#'              
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               FR = NULL, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE,
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' pop <- pop.pheno$pop
#' pop.pheno$pop <- NULL 
#' b <- cal.idx(pop.pheno = pop.pheno, index.wt = c(0.5, 0.5))
cal.idx <- function(pop.pheno, index.wt) {
  pn <- grep(pattern = "pheno", x = names(pop.pheno$info.pheno), value = TRUE)
  pheno <- do.call(data.frame, pop.pheno$info.pheno[names(pop.pheno$info.pheno) %in% pn])
  
  if(ncol(pheno) != length(index.wt)) {
    stop("Column of phenotype should equal length of weight of index")
  }
  
  # phenotype covariance matrix
  P <- var(pheno)
  iP <- try(solve(P), silent = TRUE)
  if (inherits(iP, "try-error")) {
    iP <- ginv(P)
  }
  
  pn <- grep(pattern = "TBV", x = names(pop.pheno$info.pheno), value = TRUE)
  TBV <- do.call(data.frame, pop.pheno$info.pheno[names(pop.pheno$info.pheno) %in% pn])
  # BV covariance matrix
  A <- var(TBV)
  b <- iP %*% A %*% index.wt
  b <- as.vector(b)
}

#' Get core population ID and genotype
#'
#' Build date: May 2, 2020
#' Last update: May 2, 2020
#'
#' @author Dong Yin
#'
#' @param ind.stay ID of the sires and the dams passing selection
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
#' pop1 <- getpop(nind = 100, from =   1, ratio = 0.5)
#' pop2 <- getpop(nind = 100, from = 101, ratio = 0.5, gen = 2)
#' pop.total <- rbind(pop1, pop2)
#' pop.geno.core <- matrix(1, 500, 80 )
#' pop.geno.curr <- matrix(2, 500, 100)
#' core.stay <- list(sir =   1: 30, dam =  51:100)
#'  ind.stay <- list(sir = 101:130, dam = 151:200)
#' info.core <- sel.core(ind.stay = ind.stay, core.stay = core.stay, incols = 1, 
#'     refresh = rep(0.6, 2), keep.max.gen = rep(1, 2), pop.total = pop.total, 
#'     pop.geno.core = pop.geno.core, pop.geno.curr = pop.geno.curr)
#' str(info.core)
#' 
#' # change the individual 101 to the individual 1
#' # the individual must be removed in the test
#' core.stay <- info.core$core.stay
#' pop.geno.core <- info.core$core.geno
#' core.stay$sir[1] <- 1
#' pop.geno.core[, 1] <- 1
#' 
#' pop3 <- getpop(nind = 100, from = 201, ratio = 0.5, gen = 3)
#' pop.total <- rbind(pop.total, pop3)
#' ind.stay <- list(sir = 201:230, dam = 251:300)
#' pop.geno.curr <- matrix(3, 500, 100)
#' info.core <- sel.core(ind.stay = ind.stay, core.stay = core.stay, incols = 1, 
#'     refresh = rep(0.6, 2), keep.max.gen = rep(2, 2), pop.total = pop.total, 
#'     pop.geno.curr = pop.geno.curr, pop.geno.core = pop.geno.core)
#' str(info.core)
sel.core <- function(ind.stay = NULL, core.stay = NULL, refresh = rep(0.6, 2), keep.max.gen = rep(3, 2), incols = 2, pop.total, pop.geno.curr, pop.geno.core) {
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
    id.add.sir <- sample(which(f.gen.sir == FALSE), size = num.refresh.sir-sf.gen.sir)
    f.gen.sir[id.add.sir] <- TRUE
  }
  if (num.refresh.dam > sf.gen.dam) {
    id.add.dam <- sample(which(f.gen.dam == FALSE), size = num.refresh.dam-sf.gen.dam)
    f.gen.dam[id.add.dam] <- TRUE
  }
  
  sir.stay <- ind.stay$sir[1:num.refresh.sir]
  dam.stay <- ind.stay$dam[1:num.refresh.dam]
  core.stay$sir[f.gen.sir] <- sir.stay
  core.stay$dam[f.gen.dam] <- dam.stay
  if (incols == 2) {
    gmt.dam.curr <- cal.genoloc(c(sir.stay, dam.stay), index.curr)
    gmt.dam.curr <- gmt.dam.curr * 2
    gmt.sir.curr <- gmt.dam.curr - 1
    gmt.comb.curr <- c(gmt.sir.curr, gmt.dam.curr)
    gmt.comb.curr[seq(1, length(gmt.comb.curr), 2)] <- gmt.sir.curr
    gmt.comb.curr[seq(2, length(gmt.comb.curr), 2)] <- gmt.dam.curr

    core.id <- c(core.stay$sir, core.stay$dam)
    core.id.gen <- c(core.stay$sir[which(f.gen.sir == TRUE)], core.stay$dam[which(f.gen.dam == TRUE)])
    gmt.dam.core <- cal.genoloc(core.id.gen, core.id)
    gmt.dam.core <- gmt.dam.core * 2
    gmt.sir.core <- gmt.dam.core - 1
    gmt.comb.core <- c(gmt.sir.core, gmt.dam.core)
    gmt.comb.core[seq(1, length(gmt.comb.core), 2)] <- gmt.sir.core
    gmt.comb.core[seq(2, length(gmt.comb.core), 2)] <- gmt.dam.core
  
    pop.geno.core[, gmt.comb.core] <- pop.geno.curr[, gmt.comb.curr]
    
  } else {
    gmt.comb <- cal.genoloc(c(sir.stay, dam.stay), index.curr)
    pop.geno.core[c(f.gen.sir, f.gen.dam)] <- pop.geno.curr[, gmt.comb]
  }
  
  core.list <- list(core.stay = core.stay, core.geno = pop.geno.core)
  return(core.list)
}
