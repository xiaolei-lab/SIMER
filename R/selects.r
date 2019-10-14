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
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye
#' @param decr whether to sorting with descreasing
#' @param sel.multi selection method of multi-trait with options: "tdm", "indcul" and "index"
#' @param index.wt economic weights of selection index method
#' @param index.tdm index represents which trait is being selected. NOT CONTROL BY USER
#' @param goal.perc percentage of goal more than mean of scores of individuals
#' @param pass.perc percentage of expected excellent individuals
#' @param sel.sing selection method of single trait with options: "ind", "fam", "infam" and "comb"
#' @param pop.total total population infarmation
#' @param pop.pheno list of all phenotype information
#' @param verbose whether to print detail

#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE, # single trait
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' basepop <- set.pheno(basepop, pop.pheno, sel.crit = "pheno")
#' ind.ordered1 <-
#'     selects(pop = basepop,
#'             decr = TRUE,
#'             sel.multi = "index",
#'             index.wt = c(0.5, 0.5),
#'             index.tdm = 1,
#'             goal.perc = 0.1,
#'             pass.perc = 0.9,
#'             sel.sing = "comb",
#'             pop.total = basepop,
#'             pop.pheno = pop.pheno, 
#'             verbose = TRUE)
#'
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = TRUE, # multiple traits
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' basepop <- set.pheno(basepop, pop.pheno, sel.crit = "pheno")
#' ind.ordered2 <-
#'     selects(pop = basepop,
#'             decr = TRUE,
#'             sel.multi = "index",
#'             index.wt = c(0.5, 0.5),
#'             index.tdm = 1,
#'             goal.perc = 0.1,
#'             pass.perc = 0.9,
#'             sel.sing = "comb",
#'             pop.total = basepop,
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
      
	if (ncol(pop$pheno) == 1) {
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
	    goal <- (1 + goal.perc) * apply(pop$pheno, 2, mean)
	  }
	  ind.ordered <- cal.multi(pop, decr, sel.multi, index.wt, index.tdm, goal, pass.perc, pop.pheno)
	}

	if (sel.multi == "tdm") {
    if (index.tdm == ncol(pop$pheno)) {
      logging.log("all phenotype have selected by tandem method~\n", verbose = verbose)
      index.tdm <- 1
    }
    if (pop$pheno[nrow(pop$pheno) * pass.perc, index.tdm] >= goal[index.tdm]) {
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
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE, # single trait
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' basepop <- set.pheno(basepop, pop.pheno, sel.crit = "pheno")
#' cor.r <- cal.r(basepop, basepop)
#' ind.ordered <- cal.sing(pop = basepop, decr = TRUE, sel.sing = "comb",
#'                         cor.r = cor.r)
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

      pop.tab <- summary(aov(pop$pheno~as.factor(pop$fam)))
      MB <- pop.tab[[1]][1, 3]
      MW <- pop.tab[[1]][2, 3]
      cor.n <- num.infam[[1]]
      if (is.na(MW)) {
        MW <- 0
        cor.t <- 2
      } else {
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
#' @param sel.multi selection method of multi-trait with "tdm", "indcul", "index"
#' @param index.wt economic weights of selection index method
#' @param index.tdm ndex represents which trait is being selected. NOT CONTROL BY USER
#' @param goal goal of the trait
#' @param pass.perc percentage of expected excellent individuals
#' @param pop.pheno list of all phenotype information
#'
#' @return individual indice sorted by scores
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = TRUE, # multiple traits
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' basepop <- set.pheno(basepop, pop.pheno, sel.crit = "pheno")
#' goal.perc <- 0.9
#' goal <- (1 + goal.perc) * apply(basepop$pheno, 2, mean)
#' ind.ordered <- cal.multi(pop = basepop, decr = TRUE, sel.multi = "index",
#'     index.wt = c(0.5, 0.5), index.tdm = 1, goal = goal, 
#'     pass.perc = 0.9, pop.pheno = pop.pheno)
#' str(ind.ordered)
cal.multi <- function(pop, decr, sel.multi, index.wt, index.tdm, goal, pass.perc, pop.pheno) {
  num.ind <- length(pop$index)

  if (sel.multi == "index") {
    if(ncol(pop$pheno) != length(index.wt)) {
      stop("Column of phenotype should equal length of weight of index")
    }
    
    # phenotype covariance matrix
    P <- var(pop$pheno)
    iP <- try(solve(P), silent = TRUE)
    if (inherits(iP, "try-error")) {
      iP <- ginv(P)
    }
    pn <- grep(pattern = "TBV", x = names(pop.pheno$info.pheno), value = TRUE)
    TBV <- do.call(cbind, pop.pheno$info.pheno[names(pop.pheno$info.pheno) %in% pn])
    # BV covariance matrix
    A <- var(TBV)
    b <- iP %*% A %*% index.wt
    b <- as.vector(b)
    
    ind.score <- cbind(pop$index, apply(pop$pheno * b, 1, sum))
    ind.score.ordered <- ind.score[order(ind.score[, 2], decreasing = decr), ]
    ind.ordered <- ind.score.ordered[, 1]

  } else if (sel.multi == "indcul") {
  	# goal <- (1 + goal.perc) * apply(pop$pheno, 2, mean)
  	flag.prior <- apply(pop$pheno, 1, function(v) {
  	  return(all(v > goal))
  	})
  	ind.ordered <- c(pop$index[flag.prior], pop$index[!flag.prior])

  } else if (sel.multi == "tdm") {
    ind.score <- cbind(pop$index, pop$pheno[, index.tdm])
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
#'                        pop1.geno = basepop.geno,
#'                        ind.stay = basepop$index,
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

#' Read sel_ind_geno to get selection condition
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param pop total population information
#' @param selPath the path of select_scheme
#' @param out path of output files
#'
#' @return None
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE, # single trait
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' str(basepop)
#' basepop <- set.pheno(basepop, pop.pheno, sel.crit = "pheno")
#' str(basepop)
#' selPath <- system.file("extdata", "01select_scheme", package = "simer")
#' out.index <- read.selgeno(pop = basepop, selPath = selPath, out = NULL)
read.selgeno <- function(pop=NULL, selPath=NULL, out=NULL) {
  if (!dir.exists(selPath)) stop("Please input a right selection path!")
  if (!is.null(out)) {
    if (!dir.exists(out)) {
      stop("Please input a right output path!")
    }
  }
  
  filenames <- dir(path = selPath)
  for (filename in filenames) {
    fn <- unlist(strsplit(filename, ".", fixed = TRUE))
	  if (length(fn) != 2) {
	    filenames <- filenames[!(filenames == filename)]
	    next
    }
	  if (substr(fn[1], 1, 13) != "select_scheme") {
	    filenames <- filenames[!(filenames == filename)]
	  }
  }
  
  out.index <- lapply(filenames, function(filename) {
    filename.old <- filename
    filename <- file.path(selPath, filename)
    scheme <- read.delim(filename, header = TRUE, stringsAsFactors = FALSE)
    num.row <- dim(pop)[1]
    out.index <- rep(0, num.row)
    for(i in 1:num.row) { # for2
      gen <- pop$gen[i]
      if (gen %in% scheme$generation) { # if1
        si <- which(scheme$generation == gen)
        fam.index <- scheme$family_index[si]
        if (fam.index == "all") {
          flag1 <- TRUE
        } else {
          fam.index <- num2num(fam.index) + pop$fam[pop$gen==gen][1] - 1
          flag1 <- pop$fam[i] %in% fam.index
        }
        infam.index <- scheme$within_family_index[si]
        if (infam.index == "all") {
          flag2 <- TRUE
        } else {
          infam.index <- num2num(infam.index)
          flag2 <- pop$infam[i] %in% infam.index
        }
        if (scheme$sex[si] == "all") {
          flag3 <- TRUE
        } else {
          flag3 <- pop$sex[i] == as.numeric(scheme$sex[si])
        }
        if (flag1 && flag2 && flag3) { # if2
          out.index[i] <- pop$index[i]
        } # end if2
      } # end if1
    } # end for2
    out.index <- out.index[out.index != 0]
    
    if (!is.null(out))
      write.table(out.index, file.path(out, paste0("index_", filename.old)), row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    return(out.index)
  })
  
  return(out.index)
}

#' Convert string number to numeric number
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param nums string number
#'
#' @return numeric number
#' @export
#'
#' @examples
#' nums <- "1,2,4:7"
#' num <- num2num(nums)
#' num
num2num <- function(nums) {
  nums <- unlist(strsplit(nums, split=","))
  t1 <- nums[nchar(nums) == 1]
  t2 <- nums[nchar(nums) != 1]
  if (length(t2) > 0) { # if2
    for (ti in 1:length(t2)) { # for3
      t3 <- unlist(strsplit(t2[ti], split=":"))
      t1 <- c(t1, seq(t3[1], t3[2]))
    } # end for3
  } # end if2
  nums <- sort(as.numeric(t1))
  return(nums)
}
