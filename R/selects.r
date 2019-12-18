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
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
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
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
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
  
  f1 <- grep(pattern = "TBV|TGV|pheno|ebv", x = names(pop), value = TRUE)
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
      logging.log("all phenotype have selected by tandem method~\n", verbose = verbose)
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
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
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
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
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
    TBV <- do.call(cbind, pop.pheno$info.pheno[names(pop.pheno$info.pheno) %in% pn])
    # BV covariance matrix
    A <- var(TBV)
    b <- iP %*% A %*% index.wt
    b <- as.vector(b)
    
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
#' @param selPath the path of breeding_plan
#' @param outpath the path of output files
#'
#' @return None
#' @export
#'
#' @examples
#' \donttest{
#' # get map file, map is neccessary
#' data(simdata)
#'
#' # run simer
#' simer.list <-
#'      simer(num.gen = 10,
#'            replication = 1,
#'            verbose = TRUE, 
#'            mrk.dense = FALSE,
#'            out = "simer", 
#'            outpath = NULL,
#'            out.format = "numeric",
#'            seed.geno = runif(1, 0, 100),
#'            seed.map = 12345,
#'            out.geno.gen = 3:5,
#'            out.pheno.gen = 1:5,
#'            rawgeno1 = rawgeno,
#'            rawgeno2 = NULL,
#'            rawgeno3 = NULL,
#'            rawgeno4 = NULL,
#'            num.ind = NULL,
#'            prob = c(0.5, 0.5),
#'            input.map = input.map,
#'            len.block = 5e7,
#'            range.hot = 4:6,
#'            range.cold = 1:5,
#'            rate.mut = 1e-8,
#'            cal.model = "AD",
#'            FR = NULL, 
#'            h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'            num.qtn.tr1 = 18,
#'            sd.tr1 = c(2, 1, 0.5, 0.5, 0.5, 0.1),
#'            dist.qtn.tr1 = rep("normal", 6),
#'            eff.unit.tr1 = rep(0.5, 6),
#'            shape.tr1 = rep(1, 6),
#'            scale.tr1 = rep(1, 6),
#'            multrait = FALSE,
#'            num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'            sd.trn = matrix(c(1, 0, 0, 2), 2, 2),
#'            gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'            h2.trn = c(0.3, 0.5), 
#'            qtn.spot = rep(0.1, 10),
#'            maf = 0,
#'            sel.crit = "pheno",
#'            sel.on = TRUE, 
#'            mtd.reprod = "randmate",
#'            userped = userped,
#'            num.prog = 2,
#'            ratio = 0.5,
#'            prog.tri = 2,
#'            prog.doub = 2,
#'            prog.back = rep(2, 5),
#'            ps = 0.8,
#'            decr = TRUE,
#'            sel.multi = "index",
#'            index.wt = c(0.5, 0.5),
#'            index.tdm = 1,
#'            goal.perc = 0.1,
#'            pass.perc = 0.9, 
#'            sel.sing = "comb") 
#' 
#' pop <- simer.list$pop                         
#' selPath <- system.file("extdata", "01breeding_plan", package = "simer")
#' out.pop <- read.selgeno(pop = pop, selPath = selPath, outpath = NULL)
#' str(out.pop)
#' }
read.selgeno <- function(pop=NULL, selPath=NULL, outpath=NULL) {
  if (!dir.exists(selPath)) stop("Please input a right selection path!")
  if (!is.null(outpath)) {
    if (!dir.exists(outpath)) {
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
	  if (substr(fn[1], 1, 13) != "breeding_plan") {
	    filenames <- filenames[!(filenames == filename)]
	  }
  }
  
  out.pop <- lapply(filenames, function(filename) {
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
    out.pop <- pop[pop$index %in% out.index, ]
    
    if (!is.null(outpath))
      write.table(out.pop, file.path(outpath, paste0("out_", filename.old)), row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    return(out.pop)
  })
  names(out.pop) <- paste("plan", 1:length(filenames), sep = "")
  
  return(out.pop)
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
