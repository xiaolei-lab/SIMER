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


#' Generate phenotype and other values according to selection criterion
#'
#' Build date: Jul 14, 2019
#' Last update: Oct 13, 2019
#'
#' @author Dong Yin
#'
#' @param effs a list with number of overlap markers, selected markers, effects of markers
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye
#' @param pop.geno genotype matrix of population, a individual has two columns
#' @param pos.map marker information of population
#' @param h2.tr1 heritability vector of trait1, corresponding to a, d, aXa, aXd, dXa, dXd
#' @param gnt.cov genetic covaiance matrix among all traits
#' @param env.cov environment covaiance matrix among all traits
#' @param sel.crit selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"
#' @param pop.total total population infarmation
#' @param sel.on whether to add selection
#' @param inner.env environment of main function of simer
#' @param verbose whether to print detail
#'
#' @return phenotype of population
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
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(basepop)
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
#' str(pop.pheno)
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
#'              multrait = TRUE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(basepop)
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
#' str(pop.pheno)
phenotype <-
    function(effs,
             pop,
             pop.geno,
             pos.map,
             h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
             gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
             env.cov = matrix(c(10, 5, 5, 100), 2, 2),
             sel.crit = "pheno", 
             pop.total = NULL, 
             sel.on = TRUE, 
             inner.env =  NULL, 
             verbose = TRUE) {

# Start phenotype
      
  multrait <- length(effs) > 2
  geno <- geno.cvt(pop.geno)
  nind <- ncol(geno)
  
  if (multrait) {
    nqt <- length(effs) / 2
    
    # calculate environment effects with correlation
    mat.env <- mvrnorm(n = nind, mu = rep(0, nqt), Sigma = env.cov, empirical = TRUE)
    
    mat.ind.a <- matrix(0, ncol(geno), nqt)
    for (i in 1:nqt) {
      eff.a <- effs[[2*i]]$eff.a
      qtn.a <- geno[effs[[2*i-1]], ]
      ind.a <- as.vector(crossprod(qtn.a, eff.a))
      mat.ind.a[, i] <- ind.a
    }
    if (!sel.on) {
      logging.log("build genetic correlation for traits...\n", verbose = verbose)
      # calculate additive effects with correlation
      # mat.ind.a <- mvrnorm(n = nind, mu = rep(0, nqt), Sigma = gnt.cov, empirical = TRUE)
      mat.ind.a <- build.cov(mat.ind.a, Sigma = gnt.cov)
      
      # resolve markers effects
      effs.adj <- get("effs", envir = inner.env)
      logging.log("adjust effects of markers...\n", verbose = verbose)
      for (i in 1:nqt) {
        qtn.a <- geno[effs[[2*i-1]], ]
        ind.a <- mat.ind.a[, i]
        eff.a <- c(crossprod(ginv(qtn.a), ind.a))
        effs.adj[[2*i]]$eff.a <- eff.a
      }
      assign("effs", effs.adj, envir = inner.env)
    }
    
    Covg <- var(mat.ind.a)
    Cove <- var(mat.env)
    h2 <- diag(Covg) / diag(Covg + Cove)
    gnt.cor <- cor(mat.ind.a)
    nts <- paste("tr", 1:nqt, sep = "")
    dimnames(gnt.cor) <- list(nts, nts)
    ind.pheno <- mat.ind.a + mat.env
    logging.log("Total additive    covariance matrix of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(Covg[i, ], "\n", sep = "\t", verbose = verbose)
    }
    logging.log("Total environment covariance matrix of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(Cove[i, ], "\n", sep = "\t", verbose = verbose)
    }
    logging.log("Heritability:", h2, "\n", verbose = verbose)
    logging.log("Genetic correlation of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(gnt.cor[i, ], "\n", sep = "\t", verbose = verbose)
    }
    info.tr <- list(Covg = Covg, Cove = Cove, h2 = h2)
    info.eff <- data.frame(ind.a = mat.ind.a, ind.env = mat.env)
    info.pheno <- data.frame(TBV = mat.ind.a, TGV = mat.ind.a, pheno = ind.pheno)
    pheno <- list(info.tr = info.tr, info.eff = info.eff, info.pheno = info.pheno)
    
    if (sel.crit == "TBV" | sel.crit == "TGV" | sel.crit == "pheno") {
	    pheno <- pheno
	    
	  } else {
	    if (length(effs) > 4) 
	      stop("hiblup does not support trait more than 3 for now!")
eval(parse(text = "tryCatch({
      suppressMessages(library(hiblup))
      geno.id <- pop$index
	    pheno1 <- cbind(pop$index, pheno$info.pheno$pheno.1, pheno$info.pheno$pheno.2)
	    pheno1 <- as.data.frame(pheno1)
	    geno.id <- as.data.frame(geno.id)
	    pedigree1 <- cbind(pop.total$index, pop.total$sir, pop.total$dam)
	    pos.map <- pos.map[, 1:3]
	    cal.model <- \"A\"
	    if (sel.crit == \"pEBVs\") {
	      ebv <- hiblup(pheno = pheno1, bivar.pos = c(2, 3), geno = NULL, map = pos.map, geno.id = geno.id, file.output = FALSE, 
                      pedigree = pedigree1, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else if (sel.crit == \"gEBVs\") {
	      ebv <- hiblup(pheno = pheno1, bivar.pos = c(2, 3), geno = geno, map = pos.map, geno.id = geno.id,  file.output = FALSE, 
                      pedigree = NULL, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else if (sel.crit == \"ssEBVs\") {
	      ebv <- hiblup(pheno = pheno1, bivar.pos = c(2, 3), geno = geno, map = pos.map, geno.id = geno.id,  file.output = FALSE, 
                      pedigree = pedigree1, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else {
	      stop(\"please select correct selection criterion!\")
	    }
	    idx <- pop$index
      for (i in 1:length(pop$index)){
        idx[i] <- which(ebv$ebv[, 1] == pop$index[i])
      }
      pheno$info.pheno$ebv <- ebv$ebv[idx, 2:3]
      pheno$info.hiblup <- ebv
      }, error=function(e) { 
           stop(\"Something wrong when running HIBLUP!\") })"
))
	  }

  } else {
    mrk1 <- effs$mrk1
    qtn1 <- geno[mrk1, ]
    eff1 <- effs$eff1
    len.eff <- length(eff1)
    if (len.eff == 1) {
      cal.model <- "A"
    } else if (len.eff == 2) {
      cal.model <- "AD"
    } else if (len.eff == 6) {
      cal.model <- "ADI"
    } else {
      stop("No usable cal.model!")
    }

    if (cal.model == "A") {
      pheno <-   cal.A(qtn1, h2.tr1, eff1, verbose=verbose)
    } else if (cal.model == "AD") {
      pheno <-  cal.AD(qtn1, h2.tr1, eff1, sel.on, inner.env, verbose=verbose)
    } else if (cal.model == "ADI") {
      pheno <- cal.ADI(qtn1, h2.tr1, eff1, sel.on, inner.env, verbose=verbose)
      cal.model <- "AD"
    } else {
      stop("cal.model should be 'A' or 'AD' or 'ADI'!")
    }
    if (sel.crit == "TBV" | sel.crit == "TGV" | sel.crit == "pheno") {
	    pheno <- pheno

	  } else {
eval(parse(text = "tryCatch({
      suppressMessages(library(hiblup))
      geno.id <- pop$index
	    pheno1 <- cbind(pop$index, pheno$info.pheno$pheno)
	    pheno1 <- as.data.frame(pheno1)
	    geno.id <- as.data.frame(geno.id)
	    pedigree1 <- cbind(pop.total$index, pop.total$sir, pop.total$dam)
	    pos.map <- pos.map[, 1:3]
	    if (sel.crit == \"pEBVs\") {
	      ebv <- hiblup(pheno = pheno1, geno = NULL, map = pos.map, geno.id = geno.id, file.output = FALSE, 
                      pedigree = pedigree1, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else if (sel.crit == \"gEBVs\") {
	      ebv <- hiblup(pheno = pheno1, geno = geno, map = pos.map, geno.id = geno.id, file.output = FALSE, 
                      pedigree = NULL, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else if (sel.crit == \"ssEBVs\") {
	      ebv <- hiblup(pheno = pheno1, geno = geno, map = pos.map, geno.id = geno.id, file.output = FALSE, 
                      pedigree = pedigree1, vc.method = c(\"HI\"), mode = cal.model, CV = NULL, R = NULL, snp.solution = FALSE)
	    } else {
	      stop(\"please select correct selection criterion!\")
	    }
	    idx <- pop$index
      for (i in 1:length(pop$index)){
        idx[i] <- which(ebv$ebv[, 1] == pop$index[i])
      }
      if (cal.model == \"A\") {
         pheno$info.pheno$ebv <- ebv$ebv[idx, 2]
         pheno$info.hiblup <- ebv
      } else if (cal.model == \"AD\" | cal.model == \"ADI\") {
         pheno$info.pheno$ebv <- ebv$ebv[idx, 2] + ebv$ebv[idx, 3]
         pheno$info.hiblup <- ebv
      }}, error=function(e) { 
           stop(\"Something wrong when running HIBLUP!\") })"
))
	  } # if criterion
  } # end if (multrait)

  return(pheno)
}

#' Calculate for genetic effects list of selected markers
#'
#' Build date: Sep 11, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param pop.geno genotype of population, a individual has two columns
#' @param cal.model phenotype model with "A", "AD", "ADI"
#' @param num.qtn.tr1 integer or integer vector, the number of QTN in the trait1
#' @param var.tr1 variances of different effects, the last 5 vector elements are corrresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a
#' @param dist.qtn.tr1 distribution of QTN's effects with options: "normal", "geometry" and "gamma", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param eff.unit.tr1 unit effect of geometric distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param shape.tr1 shape of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param scale.tr1 scale of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param multrait whether applying pair traits with overlapping, TRUE represents applying, FALSE represents not
#' @param num.qtn.trn QTN distribution matrix, diagnal elements are total QTN number of the trait, non-diagnal are QTN number of overlop qtn
#' @param eff.sd a matrix with the standard deviation of QTN effects
#' @param qtn.spot QTN probability in every blocks
#' @param maf Minor Allele Frequency, selection range is from  maf to 0.5
#' @param verbose whether to print detail
#'
#' @return a list with number of overlap markers, selected markers, effect of markers
#' @export
#'
#' @examples
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
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(effs)
cal.effs <-
    function(pop.geno,
             cal.model = "A",
             num.qtn.tr1 = c(2, 6, 10),
             var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
             dist.qtn.tr1 = rep("normal", 6),
             eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
             shape.tr1 = c(1, 1, 1, 1, 1, 1),
             scale.tr1 = c(1, 1, 1, 1, 1, 1),
             multrait = FALSE,
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = matrix(c(1, 0, 0, 0.5), 2, 2),
             qtn.spot = rep(0.1, 10),
             maf = 0, 
             verbose = TRUE) {

# Start calculation

  # combine odd and even columns genotype matrix
  geno <- geno.cvt(pop.geno)
  num.marker <- nrow(geno)
  num.ind <- ncol(geno)

  # calculate weight of every marker
  num.block <- length(qtn.spot)
  len.block <- num.marker %/% num.block
  tail.block <- num.marker %% num.block + len.block
  num.inblock <- c(rep(len.block, (num.block-1)), tail.block)
  wt.marker <- rep(qtn.spot / num.inblock, num.inblock)

  pop.maf <- rep(0, nrow(geno))
  for (i in 1:nrow(geno)) {
    v <- geno[i, ]
    pop.maf[i] <- min(c(sum(v == 0)+sum(v == 1)/2, sum(v == 2)+sum(v == 1)/2) / length(v))
  }
  wt.marker[pop.maf < maf] <- 0

  if (multrait) {
    if (nrow(num.qtn.trn) != ncol(num.qtn.trn) || any(num.qtn.trn != t(num.qtn.trn)))
      stop("num.qtn.trn should be symmetric matrix!")
    if (any(dim(num.qtn.trn) != dim(eff.sd))) 
      stop("non-conformable arrays!")
    
    num.qtn <- sum(num.qtn.trn[lower.tri(num.qtn.trn)]) + sum(diag(num.qtn.trn))
    sel.marker <- sample(1:num.marker, num.qtn, replace = FALSE, prob = wt.marker)
    k <- 1
    nqt <- nrow(num.qtn.trn)
    effs <- lapply(1:(2*nqt), function(i) { return(NULL) })
    names(effs) <- paste(rep(c("mrk", "eff"), nqt), rep(1:nqt, each=2), sep = "")
    for (i in 1:nqt) {
      for (j in i:nqt) {
        num.t <- num.qtn.trn[i, j]
        mrk.t <- sel.marker[k:(k+num.t-1)]
        effs[[2*i-1]] <- c(effs[[2*i-1]], mrk.t)
        if (i != j) {
          effs[[2*j-1]] <- c(effs[[2*j-1]], mrk.t)
        }
        k <- k + num.t
      }
      effs[[2*i]] <- list(eff.a = rnorm(length(effs[[2*i-1]]), 0, eff.sd[i, i]))
      logging.log("number of selected markers of trait", i, ":", length(effs[[2*i-1]]), "\n", verbose = verbose)
    }

  } else {
    num.qtn <- sum(num.qtn.tr1)
    sel.marker <- sort(sample(1:num.marker, num.qtn, replace = FALSE, prob = wt.marker))
    len.qtn <- length(num.qtn.tr1)
    logging.log("number of selected markers of trait 1:", num.qtn.tr1, "\n", verbose = verbose)
    
    if (cal.model == "A") {
      logging.log("Apply A model...\n", verbose = verbose)
      if (length(var.tr1) < length(num.qtn.tr1))
        stop("The length of var.tr1 should be no less than length of num.qtn.tr1!")
      eff.a <- cal.eff(num.qtn.tr1, var.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff1 <- list(eff.a=eff.a)
      
    } else if (cal.model == "AD") {
      logging.log("Apply AD model...\n", verbose = verbose)
      eff.a <- cal.eff(num.qtn.tr1, var.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff.d <- cal.eff(sum(num.qtn.tr1), var.tr1[len.qtn+1], dist.qtn.tr1[2], eff.unit.tr1[2], shape.tr1[2], scale.tr1[2])
      eff1 <- list(eff.a=eff.a, eff.d=eff.d)
      
    } else if (cal.model == "ADI") {
      logging.log("Apply ADI model...\n", verbose = verbose)
      if (num.qtn %% 2 != 0) stop("the number of qtn should be even in the ADI model!")
      # the first part of qtn
      ophalf <- num.qtn %/% 2
      eff.a  <- cal.eff(num.qtn.tr1,      var.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff.d  <- cal.eff(sum(num.qtn.tr1), var.tr1[len.qtn+1], dist.qtn.tr1[2], eff.unit.tr1[2], shape.tr1[2], scale.tr1[2])
      eff.aa <- cal.eff(ophalf,           var.tr1[len.qtn+2], dist.qtn.tr1[3], eff.unit.tr1[3], shape.tr1[3], scale.tr1[3])
      eff.ad <- cal.eff(ophalf,           var.tr1[len.qtn+3], dist.qtn.tr1[4], eff.unit.tr1[4], shape.tr1[4], scale.tr1[4])
      eff.da <- cal.eff(ophalf,           var.tr1[len.qtn+4], dist.qtn.tr1[5], eff.unit.tr1[5], shape.tr1[5], scale.tr1[5])
      eff.dd <- cal.eff(ophalf,           var.tr1[len.qtn+5], dist.qtn.tr1[6], eff.unit.tr1[6], shape.tr1[6], scale.tr1[6])
      eff1 <- list(eff.a=eff.a, eff.d=eff.d, eff.aa=eff.aa, eff.ad=eff.ad, eff.da=eff.da, eff.dd=eff.dd)
    }

    effs <- list(mrk1=sel.marker, eff1=eff1)
  }

  return(effs)
}

#' Calculate for genetic effects vector of selected markers
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param num.qtn number of QTN
#' @param eff.var variances of different effects
#' @param dist.qtn distribution of QTN's effects with options: "normal", "geometry" and "gamma"
#' @param eff.unit unit effect of geometric distribution
#' @param shape shape of gamma distribution
#' @param scale scale of gamma distribution
#'
#' @return genetic effects vector of selected markers
#' @export
#'
#' @examples
#' num.qtn <- c(2, 6, 10) # three qtn groups
#' eff.var <- c(0.4, 0.2, 0.02) # three variance of qtn group
#' eff <- cal.eff(num.qtn = num.qtn, eff.var = eff.var, dist.qtn = "normal",
#'         eff.unit = 0.5, shape = 1, scale = 1)
#' str(eff)
#'
#' num.qtn <- sum(num.qtn)
#' eff.var <- sum(eff.var)
#' eff <- cal.eff(num.qtn = num.qtn, eff.var = eff.var, dist.qtn = "normal",
#'         eff.unit = 0.5, shape = 1, scale = 1)
#' str(eff)
cal.eff <- function(num.qtn, eff.var, dist.qtn, eff.unit, shape, scale) {
  if (sum(num.qtn) == 0) return(0)
  # Judge which kind of distribution of QTN
  eff.qtn <- NULL
  if(dist.qtn == "normal") {
    for (nq in 1:length(num.qtn)) {
    	eff.qtn <- c(eff.qtn, rnorm(num.qtn[nq], 0, eff.var[nq]))
    }

  } else if(dist.qtn == "geometry") {
    for (nq in 1:length(num.qtn)) {
    	eff.qtn <- c(eff.qtn, eff.unit^(1:num.qtn[nq]))
    }

  } else if(dist.qtn == "gamma") {
    for (nq in 1:length(num.qtn)) {
    	eff.qtn <- c(eff.qtn, rgamma(num.qtn[nq], shape, scale))
    }

  } else {
    stop("please input a right QTN effect!")
  }

  return(eff.qtn)
}

#' Calculate for phenotype of A model
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param qtn QTN matrix of population
#' @param h2 heritability of trait
#' @param eff1 effects list of QTN
#' @param ind.a additive effects of QTN
#' @param verbose whether to show all information
#'
#' @return phenotype of A model
#' @export
#'
#' @examples
#' num.qtn <- c(2, 6, 10) # three qtn groups
#' nqtn <- sum(num.qtn) # total number of qtn
#' nind <- 100
#' qtn <- matrix(sample(c(0, 1, 2), nqtn*nind, replace=TRUE), nqtn, nind)
#' eff.a <- cal.eff(num.qtn = num.qtn, eff.var = c(0.4, 0.2, 0.02),
#'                  dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff1 <- list(eff.a = eff.a)
#' h2 <- c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01)
#' pheno.a1 <- cal.A(qtn = qtn, h2 = h2, eff1 = eff1, ind.a = NULL, verbose = TRUE)
#' str(pheno.a1)
#'
#' ind.a <- crossprod(qtn, eff.a)
#' pheno.a2 <- cal.A(qtn = NULL, h2 = h2, eff1 = NULL, ind.a = ind.a, verbose = TRUE)
#' str(pheno.a2)
cal.A <- function(qtn, h2, eff1, ind.a=NULL, verbose=FALSE) {
  # change code from (0, 1, 2) to (-1, 0, 1)
  # qtn.a <- qtn
  qtn.a <- qtn - 1
  eff.a <- eff1$eff.a
  if (is.null(ind.a))
    ind.a <- as.vector(crossprod(qtn.a, eff.a))

  num.ind <- length(ind.a)
  var.add <- var(ind.a)
  if (verbose) logging.log("Total additive    variance:", var.add, "\n", verbose = verbose)
  if (h2[1] > 0 & h2[1] <=1) {
    var.env <- (var.add - h2[1]*var.add) / h2[1]
  } else if (h2[1] == 0){
    var.env <- 1
    ind.a <- ind.a * 0
  } else {
    stop("Heritability of additive should be no more than 1 and no less than 0!")
  }
  ind.env <- rnorm(num.ind, 0, sqrt(var.env))
  ind.pheno <- ind.a + ind.env
  Vg <- var.add
  Ve <- var(ind.env)
  h2.new <- Vg / (Vg + Ve)
  logging.log("Total environment variance:", Ve, "\n", verbose = verbose)
  logging.log("Heritability:", h2.new, "\n", verbose = verbose)
  info.tr <- list(Vg = Vg, Ve = Ve, h2 = h2.new)
  info.eff <- data.frame(ind.a = ind.a, ind.env = ind.env)
  info.pheno <- data.frame(TBV = ind.a, TGV = ind.a, pheno = ind.pheno)
  A.list <- list(info.tr = info.tr, info.eff = info.eff, info.pheno = info.pheno, qtn.a = qtn.a)
  return(A.list)
}

#' Calculate for phenotype of AD model
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param qtn QTN matrix of population
#' @param h2 heritability of trait
#' @param eff1 effects list of QTN
#' @param sel.on whether to add selection
#' @param inner.env environment of main function of simer
#' @param verbose whether to show all information
#'
#' @return phenotype of AD model
#' @export
#'
#' @examples
#' num.qtn <- c(2, 6, 10) # three qtn groups
#' len.qtn <- length(num.qtn)
#' nqtn <- sum(num.qtn) # total number of qtn
#' nind <- 100
#' qtn <- matrix(sample(c(0, 1, 2), nqtn*nind, replace=TRUE), nqtn, nind)
#' var.tr1 <- c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001)
#' eff.a <- cal.eff(num.qtn = c(2, 6, 10), eff.var = var.tr1[1:len.qtn],
#'                  dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.d <- cal.eff(num.qtn = nqtn,        eff.var = var.tr1[len.qtn+1],
#'                  dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff1 <- list(eff.a = eff.a, eff.d = eff.d)
#' h2 <- c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01)
#' pheno.ad <- cal.AD(qtn = qtn, h2 = h2, eff1 = eff1, sel.on = TRUE, inner.env = NULL, verbose = TRUE)
#' str(pheno.ad)
cal.AD <- function(qtn, h2, eff1, sel.on = TRUE, inner.env = NULL, verbose=FALSE) {
  A.list <- cal.A(qtn, h2, eff1, verbose=FALSE)
  ind.a <- A.list$info.eff$ind.a

  # change code from (0, 1, 2) to (-0.5, 0.5, -0.5)
  qtn.d <- qtn
  qtn.d[qtn.d == 2] <- 0
  qtn.d <- qtn.d - 0.5
  eff.d <- eff1$eff.d
  ind.d <- as.vector(crossprod(qtn.d, eff.d))

  num.ind <- length(ind.a)
  var.add <- var(ind.a)
  
  if (!sel.on) {
    var.dom <- var(ind.d)
    if (var.dom != 0) {
      # adjust domimance effect according to ratio of additive variance and dominance variance
      ind.d <- ind.d * sqrt(h2[2] / var.dom * var.add / h2[1])
      eff.d <- c(crossprod(ginv(qtn.d), ind.d))
    }
    effs.adj <- get("effs", envir = inner.env)
    logging.log("adjust effects of markers...\n", verbose = verbose)
    effs.adj$eff1$eff.d <- eff.d
    assign("effs", effs.adj, envir = inner.env)
  }
  
  ind.gv <- ind.a + ind.d
  var.gv <- var(ind.gv)

  if (verbose) {
    logging.log("Total additive    variance:", var.add, "\n", verbose = verbose)
    logging.log("Total dominance   variance:", var(ind.d), "\n", verbose = verbose)
    logging.log("Total genetic     variance:", var.gv, "\n", verbose = verbose)
  }

  H2 <- h2[1] + h2[2]
  if (H2 > 0 & H2 <=1) {
    var.env <- (var.gv - H2 * var.gv) / H2
  } else if (H2 == 0){
    var.env <- 1
    ind.gv <- ind.gv * 0
  } else {
    stop("Heritability of additive and dominance should be no more than 1 and no less than 0!")
  }
  ind.env <- rnorm(num.ind, 0, sqrt(var.env))
  ind.pheno <- ind.gv + ind.env
  Vg <- c(var.add, var(ind.d))
  Ve <- var(ind.env)
  h2.new <- Vg / sum(Vg, Ve)
  logging.log("Total environment variance:", Ve, "\n", verbose = verbose)
  logging.log("Heritability:", h2.new, "\n", verbose = verbose)
  info.tr <- list(Vg = Vg, Ve = Ve, h2 = h2.new)
  info.eff <- data.frame(ind.a = ind.a, ind.d = ind.d, ind.env = ind.env)
  info.pheno <- data.frame(TBV = ind.a, TGV = ind.gv, pheno = ind.pheno)
  AD.list <- list(info.tr = info.tr, info.eff = info.eff, info.pheno = info.pheno, qtn.a = A.list$qtn.a, qtn.d = qtn.d)
  return(AD.list)
}

#' Calculate for phenotype of ADI model
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param qtn QTN matrix of population
#' @param h2 heritability of trait
#' @param eff1 effects list of QTN
#' @param sel.on whether to add selection
#' @param inner.env environment of main function of simer
#' @param verbose whether to show all information
#'
#' @return phenotype of ADI model
#' @export
#' @references Kao C and Zeng Z (2002) <https://www.genetics.org/content/160/3/1243.long>
#'
#' @examples
#' num.qtn <- c(2, 6, 10) # three qtn groups
#' nqtn <- sum(num.qtn)
#' len.qtn <- length(num.qtn)
#' ophalf <- nqtn %/% 2
#' nind <- 100
#' qtn <- matrix(sample(c(0, 1, 2), nqtn*nind, replace=TRUE), nqtn, nind)
#' var.tr1 <- c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001)
#' eff.a  <- cal.eff(num.qtn = num.qtn, eff.var = var.tr1[1:len.qtn],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.d  <- cal.eff(num.qtn =    nqtn, eff.var = var.tr1[len.qtn+1],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.aa <- cal.eff(num.qtn =  ophalf, eff.var = var.tr1[len.qtn+2],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.ad <- cal.eff(num.qtn =  ophalf, eff.var = var.tr1[len.qtn+3],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.da <- cal.eff(num.qtn =  ophalf, eff.var = var.tr1[len.qtn+4],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff.aa <- cal.eff(num.qtn =  ophalf, eff.var = var.tr1[len.qtn+5],
#'                   dist.qtn = "normal", eff.unit = 0.5, shape = 1, scale = 1)
#' eff1 <- list(eff.a = eff.a, eff.d = eff.d, eff.aa = eff.aa,
#'              eff.ad = eff.ad, eff.da = eff.da, eff.dd = eff.ad)
#' h2 <- c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01)
#' pheno.adi <- cal.ADI(qtn = qtn, h2 = h2, eff1 = eff1, sel.on = TRUE, 
#'     inner.env = NULL, verbose = TRUE)
#' str(pheno.adi)
cal.ADI <- function(qtn, h2, eff1, sel.on = TRUE, inner.env = NULL, verbose=FALSE) {
  AD.list <- cal.AD(qtn, h2, eff1, sel.on, inner.env, verbose=FALSE)
  ind.a <- AD.list$info.eff$ind.a
  ind.d <- AD.list$info.eff$ind.d
  qtn.a <- AD.list$qtn.a
  qtn.d <- AD.list$qtn.d
  
  # only for ADI model
  # qtn.a <- qtn.a - 1
  # qtn.d <- qtn.d - 0.5

  # two part of qtn
  ophalf <- 1:(nrow(qtn) %/% 2)
  edhalf <- ((nrow(qtn) %/% 2)+1):nrow(qtn)

  num.ind <- length(ind.a)
  var.add <- var(ind.a)
  qtn.aa <- qtn.a[ophalf, ] * qtn.a[edhalf, ]
  qtn.ad <- qtn.a[ophalf, ] * qtn.d[edhalf, ]
  qtn.da <- qtn.d[ophalf, ] * qtn.a[edhalf, ]
  qtn.dd <- qtn.d[ophalf, ] * qtn.d[edhalf, ]
  eff.aa <- eff1$eff.aa
  eff.ad <- eff1$eff.ad
  eff.da <- eff1$eff.da
  eff.dd <- eff1$eff.dd
  ind.aa <- as.vector(crossprod(qtn.aa, eff.aa))
  ind.ad <- as.vector(crossprod(qtn.ad, eff.ad))
  ind.da <- as.vector(crossprod(qtn.da, eff.da))
  ind.dd <- as.vector(crossprod(qtn.dd, eff.dd))
  var.aa <- var(ind.aa)
  var.ad <- var(ind.ad)
  var.da <- var(ind.da)
  var.dd <- var(ind.dd)

  if (!sel.on) {
    # adjust iteraction effect according to ratio of additive variance and iteraction variance
    if (var.aa != 0) {
      ind.aa <- ind.aa * sqrt(h2[3] / var.aa * var.add / h2[1])
      eff.aa <- c(crossprod(ginv(qtn.aa), ind.aa))
    }
    if (var.ad != 0) {
      ind.ad <- ind.ad * sqrt(h2[4] / var.ad * var.add / h2[1])
      eff.ad <- c(crossprod(ginv(qtn.ad), ind.ad))
    }
    if (var.da != 0) {
      ind.da <- ind.da * sqrt(h2[5] / var.da * var.add / h2[1])
      eff.da <- c(crossprod(ginv(qtn.da), ind.da))
    }
    if (var.dd != 0) {
      ind.dd <- ind.dd * sqrt(h2[6] / var.dd * var.add / h2[1])
      eff.dd <- c(crossprod(ginv(qtn.dd), ind.dd))
    } 
    effs.adj <- get("effs", envir = inner.env)
    logging.log("adjust effects of markers...\n", verbose = verbose)
    effs.adj$eff1$eff.aa <- eff.aa
    effs.adj$eff1$eff.ad <- eff.ad
    effs.adj$eff1$eff.da <- eff.da
    effs.adj$eff1$eff.dd <- eff.dd
    assign("effs", effs.adj, envir = inner.env)
  }
  
  ind.gv <- ind.a + ind.d + ind.aa + ind.ad + ind.da + ind.dd
  var.gv <- var(ind.gv)

  if (verbose) {
    logging.log("Total additive            variance:", var.add, "\n", verbose = verbose)
    logging.log("Total dominance           variance:", var(ind.d), "\n", verbose = verbose)
    logging.log("Total additiveXadditive   variance:", var(ind.aa), "\n", verbose = verbose)
    logging.log("Total dominanceXadditive  variance:", var(ind.ad), "\n", verbose = verbose)
    logging.log("Total dominanceXadditive  variance:", var(ind.da), "\n", verbose = verbose)
    logging.log("Total dominanceXdominance variance:", var(ind.dd), "\n", verbose = verbose)
    logging.log("Total genetic             variance:", var.gv, "\n", verbose = verbose)
  }

  H2 <- sum(h2)
  if (H2 > 0 & H2 <=1) {
    var.env <- (var.gv - H2 * var.gv) / H2
  } else if (H2 == 0){
    var.env <- 1
    ind.gv <- ind.gv * 0
  } else {
    stop("Heritability of additive and dominance and interaction of them should be no more than 1 and no less than 0!")
  }
  ind.env <- rnorm(num.ind, 0, sqrt(var.env))
  ind.pheno <- ind.gv + ind.env
  Vg <- c(var.add, var(ind.d), var(ind.aa), var(ind.ad), var(ind.da), var(ind.dd))
  Ve <- var(ind.env)
  h2.new <- Vg / sum(Vg, Ve)
  logging.log("Total environment         variance:", Ve, "\n", verbose = verbose)
  logging.log("Heritability:", h2.new, "\n", verbose = verbose)
  info.tr <- list(Vg = Vg, Ve = Ve, h2 = h2.new)
  info.eff <- data.frame(ind.a = ind.a, ind.d = ind.d, ind.aa = ind.aa, ind.ad = ind.ad, ind.da = ind.da, ind.dd = ind.dd, ind.env = ind.env)
  info.pheno <- data.frame(TBV = ind.a, TGV = ind.gv, pheno = ind.pheno)
  ADI.list <- list(info.tr = info.tr, info.eff = info.eff, info.pheno = info.pheno)
  return(ADI.list)
}

#' Convert genotype matrix from (0, 1) to (0, 1, 2)
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param pop.geno genotype matrix of (0, 1)
#'
#' @return genotype matrix of (0, 1, 2)
#' @export
#'
#' @examples
#' num.marker <- 48353
#' num.ind <- 100
#' geno1 <- genotype(num.marker = num.marker, num.ind = num.ind, verbose = TRUE)
#' geno2 <- geno.cvt(pop.geno = geno1)
#' geno1[1:5, 1:10]
#' geno2[1:5, 1:5]
geno.cvt <- function(pop.geno) {
  num.ind <- ncol(pop.geno) / 2
  v.odd <- (1:num.ind) * 2 - 1
  v.even <- (1:num.ind) * 2
  geno <- pop.geno[, v.odd] + pop.geno[, v.even]
  return(geno)
}

#' To bulid correlation of variables
#'
#' Build date: Oct 10, 2019
#' Last update: Oct 10, 2019
#'
#' @author Dong Yin and R
#'
#' @param mat matrix without correlation
#' @param mu means of the variables 
#' @param Sigma covariance matrix of variables
#' @param tol tolerance (relative to largest variance) for numerical 
#'       lack of positive-definiteness in Sigma.
#'
#' @return matrix with correlaion
#' @export
#' @references B. D. Ripley (1987) Stochastic Simulation. Wiley. Page 98
#'
#' @examples
#' Sigma <- matrix(c(14, 10, 10, 15), 2, 2)
#' Sigma
#' mat <- cbind(rnorm(100), rnorm(100))
#' mat.cov <- build.cov(mat, Sigma = Sigma)
#' var(mat.cov)
build.cov <- function(mat, mu = rep(0, nrow(Sigma)), Sigma, tol = 1e-06) {
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  
  mat <- scale(mat, center = TRUE, scale = FALSE)
  mat <- mat %*% svd(mat, nu = 0)$v
  mat <- scale(mat, center = FALSE, scale = TRUE)

  mat <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(mat)
  mat <- t(mat)
  return(mat)
}

#' Set the phenotype of the population
#'
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex
#' @param pop.pheno phenotype information
#' @param sel.crit selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"
#'
#' @return popluation information with phenotype
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
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              eff.sd = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(basepop)
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
#' basepop <- set.pheno(pop = basepop, pop.pheno = pop.pheno, sel.crit = "pheno")
#' str(basepop)
set.pheno <- function(pop, pop.pheno, sel.crit) {
  if (sel.crit == "TBV") {
    pn <- grep(pattern = "TBV", x = names(pop.pheno$info.pheno), value = TRUE)
  } else if (sel.crit == "TGV") {
    pn <- grep(pattern = "TGV", x = names(pop.pheno$info.pheno), value = TRUE)
  } else if (sel.crit == "pheno") {
    pn <- grep(pattern = "pheno", x = names(pop.pheno$info.pheno), value = TRUE)
  } else if (sel.crit == "pEBVs" | sel.crit == "gEBVs" | sel.crit == "ssEBVs") {
    pn <- grep(pattern = "ebv", x = names(pop.pheno$info.pheno), value = TRUE)
  } else {
    stop("please select correct selection criterion!")
  }
  pop$pheno <- do.call(cbind, pop.pheno$info.pheno[names(pop.pheno$info.pheno) %in% pn])
  return(pop)
}
