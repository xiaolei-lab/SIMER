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
#' Last update: Nov 5, 2019
#'
#' @author Dong Yin
#'
#' @param effs a list of selected markers and their effects
#' @param FR list of fixed effects, random effects, and their combination
#' @param pop population information of generation, family ID, within-family ID, individual ID, paternal ID, maternal ID, and sex
#' @param pop.geno genotype matrix of the population, an individual has two columns
#' @param pos.map marker information of the population
#' @param h2.tr1 heritability vector of a single trait, every element are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param gnt.cov genetic covaiance matrix among all traits
#' @param h2.trn heritability among all traits
#' @param sel.crit selection criteria with the options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", and "pheno"
#' @param pop.total total population infarmation
#' @param sel.on whether to add selection
#' @param inner.env R environment of parameter "effs"
#' @param verbose whether to print detail
#'
#' @return phenotype of population
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 49336, num.ind = 100, verbose = TRUE)
#' a <- sample(c("a1", "a2", "a3"), 100, replace = TRUE)
#' b <- sample(c("b1", "b2", "b3"), 100, replace = TRUE)
#' pop$a <- a # load your fixed  effects
#' pop$b <- b # load your random effects
#' pop.env <- environment()
#' 
#' # combination of fixed effects
#' cmb.fix <- list(tr1 = c("mu", "gen", "sex"), # trait 1
#'                 tr2 = c("mu", "diet", "season")) # trait 2
#'             
#' fix <- list( 
#'          mu = list(level = "mu", eff = 2),  				      
#'         gen = FALSE,         
#'         sex = list(level = c("1", "2"), eff = c(0.5, 0.3)), 
#'        diet = list(level = c("d1", "d2", "d3"), eff = c(0.1, 0.2, 0.3)),
#'      season = list(level = c("s1", "s2", "s3", "s4"), eff = c(0.1, 0.2, 0.3, 0.2)), 
#'           a = list(level = c("a1", "a2", "a3"), eff = c(0.1, 0.2, 0.3))) 
#' 
#' # combination and ralation of random effects
#' tr1 <- list(rn = c("sir", "PE"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))
#' tr2 <- list(rn = c("dam", "litter"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))          
#' cmb.rand <- list(tr1 = tr1, tr2 = tr2)   
#'       
#' rand <- list(
#'         sir = list(mean = 0, sd = 1),		      
#'         dam = list(mean = 0, sd = 1),	       
#'          PE = list(level = c("p1", "p2", "p3"), eff = c(0.01, 0.02, 0.03)), 
#'      litter = list(level = c("l1", "l2"), eff = c(0.01, 0.02)),    
#'           b = list(level = c("b1", "b2", "b3"), eff = c(0.01, 0.02, 0.03)))
#'   
#' FR <- list(cmb.fix = cmb.fix, fix = fix, cmb.rand = cmb.rand, rand = rand)
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
#'               FR = FR, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5),  
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = TRUE, 
#'               inner.env = pop.env, 
#'               verbose = TRUE)
#' str(pop)              
#' pop <- pop.pheno$pop
#' str(pop)
#' pop.pheno$pop <- NULL           
#' str(pop.pheno)
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
#'               FR = FR, 
#'               pop = pop,
#'               pop.geno = pop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5), 
#'               sel.crit = "pheno", 
#'               pop.total = pop, 
#'               sel.on = FALSE, 
#'               inner.env = pop.env, 
#'               verbose = TRUE)
#' str(pop)              
#' pop <- pop.pheno$pop
#' str(pop)
#' pop.pheno$pop <- NULL  
#' str(pop.pheno)
phenotype <-
    function(effs = NULL,
             FR = NULL, 
             pop = NULL,
             pop.geno = NULL,
             pos.map = NULL,
             h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
             gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
             h2.trn = c(0.3, 0.5), 
             sel.crit = "pheno", 
             pop.total = NULL, 
             sel.on = TRUE, 
             inner.env =  NULL, 
             verbose = TRUE) {

# Start phenotype
  
  pop.env <- environment()
  if (is.null(pop.geno)) {
    multrait <- effs
  } else {
    multrait <- length(effs) > 2
  }
  geno <- geno.cvt(pop.geno)
  nind <- nrow(pop)
  
  if (multrait) {
    nqt <- nrow(gnt.cov)
    
    df.ind.a <- lapply(1:nqt, function(i) { return(rep(0, nind)) })
    df.ind.a <- as.data.frame(df.ind.a)
    nts <- paste("tr", 1:nqt, sep = "")
    names(df.ind.a) <- nts
    
    if (!is.null(geno)) {
      # calculate for corelated additive effects
      for (i in 1:nqt) {
        eff.a <- effs[[2*i]]$eff.a
        qtn.a <- geno[effs[[2*i-1]], ]
        ind.a <- as.vector(crossprod(qtn.a, eff.a))
        df.ind.a[, i] <- ind.a
      }
      if (!sel.on) {
        logging.log("build genetic correlation for traits...\n", verbose = verbose)
        # calculate additive effects with correlation
        # df.ind.a <- mvrnorm(n = nind, mu = rep(0, nqt), Sigma = gnt.cov, empirical = TRUE)
        df.ind.a <- build.cov(df.ind.a, Sigma = gnt.cov)
        
        # adjust markers effects
        effs.adj <- get("effs", envir = inner.env)
        logging.log("adjust effects of markers...\n", verbose = verbose)
        for (i in 1:nqt) {
          qtn.a <- geno[effs[[2*i-1]], ]
          ind.a <- df.ind.a[, i]
          eff.a <- c(crossprod(ginv(qtn.a), ind.a))
          effs.adj[[2*i]]$eff.a <- eff.a
        }
        assign("effs", effs.adj, envir = inner.env)
      }
    } # end if (!is.null(geno))
    
    var.add <- var(df.ind.a)
    var.pheno <- diag(var.add) / h2.trn
    if (is.null(geno)) var.pheno <- sample(100, nqt)
    
    # calculate for fixed effects and random effects
    if (!is.null(FR)) {
      # calculate for fixed effects and random effects
      fr <- cal.FR(pop = pop, FR = FR, var.pheno = var.pheno, pop.env = pop.env, verbose = verbose)
      frn <- names(fr)
      var.fr <- unlist(lapply(1:length(fr), function(i) {  return(sum(apply(fr[[i]]$rand, 2, var))) }))
      fr <- lapply(1:length(fr), function(i) {return(cbind(fr[[i]]$fix, fr[[i]]$rand))})
      names(fr) <- frn
      mat.fr <- do.call(cbind, lapply(1:length(fr), function(i) { return(apply(fr[[i]], 1, sum)) }))

      # calculate for environmental effects
      var.env <- var.pheno - var.fr - diag(var.add)
      mat.env <- mvrnorm(n = nind, mu = rep(0, length(var.env)), Sigma = diag(var.env), empirical = TRUE)
      mat.env <- as.data.frame(mat.env)
      names(mat.env) <- nts
      var.env <- var(mat.env)
      h2 <- diag(var.add) / (diag(var.add) + var.fr + diag(var.env))
      
      # calculate for phenotype
      ind.pheno <- df.ind.a + mat.fr + mat.env
      
    } else {
      fr <- lapply(1:nqt, function(i) return(NULL))
      names(fr) <- paste("tr", 1:nqt, sep = "")
      
      # calculate for environmental effects
      var.env <- var.pheno - diag(var.add)
      mat.env <- mvrnorm(n = nind, mu = rep(0, length(var.env)), Sigma = diag(var.env), empirical = TRUE)
      mat.env <- as.data.frame(mat.env)
      names(mat.env) <- nts
      var.env <- var(mat.env)
      h2 <- diag(var.add) / (diag(var.add) + diag(var.env))
      
      # calculate for phenotype
      ind.pheno <- df.ind.a + mat.env
    }

    var.env <- var(mat.env)
    if(any(var(df.ind.a) == 0)) {
      gnt.cor <- matrix(0, nqt, nqt)
    } else {
      gnt.cor <- cor(df.ind.a)
    }
    logging.log("Total additive    covariance matrix of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(var.add[i, ], "\n", sep = "\t", verbose = verbose)
    }
    logging.log("Total environment covariance matrix of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(var.env[i, ], "\n", sep = "\t", verbose = verbose)
    }
    logging.log("Heritability:", h2, "\n", verbose = verbose)
    logging.log("Genetic correlation of all traits: \n", verbose = verbose)
    for (i in 1:nqt) {
      logging.log(gnt.cor[i, ], "\n", sep = "\t", verbose = verbose)
    }
    info.tr <- list(Covg = var.add, Cove = var.env, h2 = h2, gnt.cor = gnt.cor)
    for (i in 1:nqt) {
      fr[[i]]$ind.a <- df.ind.a[, i]
      fr[[i]]$ind.env <- mat.env[, i]
      if (!is.data.frame(fr[[i]])) fr[[i]] <- as.data.frame(fr[[i]])
    }
    info.pheno <- data.frame(TBV = df.ind.a, TGV = df.ind.a, pheno = ind.pheno)
    pheno <- list(info.tr = info.tr, info.eff = fr, info.pheno = info.pheno)
    
    # check data quality
    for (i in 1:nqt) {
      idx.len <- unlist(lapply(1:ncol(pheno$info.eff[[i]]), function(j) {  return(length(unique(pheno$info.eff[[i]][, j]))) }))
      info.eff.t <- pheno$info.eff[[i]][, idx.len != 1]
      info.eff.cor <- cor(info.eff.t)
      info.eff.f <- names(info.eff.t)[names(info.eff.t) %in% c("ind.d", "ind.aa", "ind.ad", "ind.da", "ind.dd")]
      info.eff.cor[info.eff.f, ] <- 0
      info.eff.cor[, info.eff.f] <- 0
      if (any(info.eff.cor[lower.tri(info.eff.cor)] > 0.5))
        warning("There are hign-correlations between fixed effects or fixed effects and random effects, and it will reduce the accuracy of effects simulation!")
    }
    
    if (sel.crit == "TBV" | sel.crit == "TGV" | sel.crit == "pheno") {
	    pheno <- pheno
	    
	  } else {
eval(parse(text = "tryCatch({
      suppressMessages(library(hiblup))
      geno.id <- as.data.frame(pop$index)
      pn <- grep(pattern = \"pheno\", x = names(pheno$info.pheno), value = TRUE)
      pheno1 <- subset(pheno$info.pheno, select = pn)
      pheno1 <- cbind(pop$index, pheno1)
      pedigree1 <- subset(pop.total, select = c(\"index\", \"sir\", \"dam\"))
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
      pheno$info.pheno$ebv <- ebv$ebv[idx, 2:ncol(ebv$ebv)]
      pheno$info.hiblup <- ebv
      }, error=function(e) { 
           stop(\"Something wrong when running HIBLUP!\") })"
))
	  }

  } else {
    # calculate for genetic effects
    info.eff <- cal.gnt(geno = geno, h2 = h2.tr1, effs = effs, sel.on = sel.on, inner.env = inner.env, verbose = verbose)
    
    if (!is.null(info.eff)) {
      # calculate for phenotype variance
      ind.a <- info.eff$ind.a
      var.pheno <- var(ind.a) / h2.tr1[1]
    } else {
      ind.a <- rep(0, nind)
      var.pheno <- sample(100, 1)
    }
    
    # calculate for fixed effects and random effects
    if (!is.null(FR)) {
      fr <- cal.FR(pop = pop, FR = FR, var.pheno = var.pheno, pop.env = pop.env, verbose = verbose)
    } else {
      fr <- NULL
    }

    # calculate for phenotype
    pheno <- cal.pheno(fr = fr, info.eff = info.eff, h2 = h2.tr1, num.ind = nind, var.pheno = var.pheno, verbose = verbose)
    
    # check data quality
    idx.len <- unlist(lapply(1:ncol(pheno$info.eff), function(i) {  return(length(unique(pheno$info.eff[, i]))) }))
    info.eff.t <- pheno$info.eff[, idx.len != 1]
    info.eff.cor <- cor(info.eff.t)
    info.eff.f <- names(info.eff.t)[names(info.eff.t) %in% c("ind.d", "ind.aa", "ind.ad", "ind.da", "ind.dd")]
    info.eff.cor[info.eff.f, ] <- 0
    info.eff.cor[, info.eff.f] <- 0
    if (any(info.eff.cor[lower.tri(info.eff.cor)] > 0.5))
      warning("There are hign-correlations between fixed effects or fixed effects and random effects, and it will reduce the accuracy of effects simulation!")
    
    if (sel.crit == "TBV" | sel.crit == "TGV" | sel.crit == "pheno") {
	    pheno <- pheno

	  } else {
eval(parse(text = "tryCatch({
      suppressMessages(library(hiblup))
      geno.id <- as.data.frame(pop$index)
      pn <- grep(pattern = \"pheno\", x = names(pheno$info.pheno), value = TRUE)
      pheno1 <- subset(pheno$info.pheno, select = pn)
      pheno1 <- cbind(pop$index, pheno1)
	    pedigree1 <- subset(pop.total, select = c(\"index\", \"sir\", \"dam\"))
	    pos.map <- pos.map[, 1:3]
	    if (\"ind.d\" %in% names(pheno$info.eff)) {
	      cal.model <- \"AD\"
	    } else {
	      cal.model <- \"A\"
	    }
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
  
  pop <- set.pheno(pop, pheno, sel.crit)
  pheno$pop <- pop
  return(pheno)
}

#' Generate genetic effects
#'
#' Build date: Nov 3, 2019
#' Last update: Nov 3, 2019
#'
#' @author Dong Yin
#'
#' @param geno genotype matrix of the population, an individual has two columns
#' @param h2 heritability vector of the trait, every elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param effs a list of selected markers and their effects
#' @param sel.on whether to add selection
#' @param inner.env R environment of parameter "effs"
#' @param verbose whether to print detail
#'
#' @return phenotype of population
#' @export
#' @references Kao C and Zeng Z (2002) <https://www.genetics.org/content/160/3/1243.long>
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' geno <- geno.cvt(basepop.geno)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
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
#' info.eff <- cal.gnt(geno = geno, h2 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01), 
#'     effs = effs, sel.on = TRUE, inner.env = NULL, verbose = TRUE)
#' str(info.eff)
cal.gnt <- function(geno = NULL, h2 = NULL, effs = NULL, sel.on = TRUE, inner.env = NULL, verbose = TRUE) {
  if (is.null(geno)) return(NULL)
  
  mrk1 <- effs$mrk1
  qtn1 <- geno[mrk1, ]
  eff1 <- effs$eff1
  len.eff <- length(eff1)
  
  # change code from (0, 1, 2) to (-1, 0, 1)
  qtn.a <- qtn1 - 1
  eff.a <- eff1$eff.a
  ind.a <- as.vector(crossprod(qtn.a, eff.a))
  info.eff <- data.frame(ind.a = ind.a)
  var.add <- var(ind.a)
  if (verbose) logging.log("Total additive            variance:", var.add, "\n", verbose = verbose)
  
  if (h2[1] > 0 & h2[1] <=1) {
    var.pheno <- var.add / h2[1]
  } else if (h2[1] == 0){
    var.pheno <- 1
    ind.a <- ind.a * 0
    var.add <- 0
  } else {
    stop("Heritability of additive should be no more than 1 and no less than 0!")
  }
  
  # dominance effect
  if (len.eff >= 2) {
    # change code from (0, 1, 2) to (-0.5, 0.5, -0.5)
    qtn.d <- qtn1
    qtn.d[qtn.d == 2] <- 0
    qtn.d <- qtn.d - 0.5
    eff.d <- eff1$eff.d
    ind.d <- as.vector(crossprod(qtn.d, eff.d))
    var.dom <- var(ind.d)

    if (!sel.on) {
      if (var.dom != 0) {
        # adjust domimance effect according to ratio of additive variance and dominance variance
        ind.d <- ind.d * sqrt(var.pheno * h2[2] / var.dom)
        eff.d <- c(crossprod(ginv(qtn.d), ind.d))
        var.dom <- var(ind.d)
      }
      effs.adj <- get("effs", envir = inner.env)
      logging.log("adjust dominance effects of markers...\n", verbose = verbose)
      effs.adj$eff1$eff.d <- eff.d
      assign("effs", effs.adj, envir = inner.env)
    }
    
    info.eff$ind.d <- ind.d
    if (verbose) logging.log("Total dominance           variance:", var.dom, "\n", verbose = verbose)
  }
  
  # interaction effect
  if (len.eff == 6) {
    # two part of qtn
    ophalf <- 1:(nrow(qtn.a) %/% 2)
    edhalf <- ((nrow(qtn.a) %/% 2)+1):nrow(qtn.a)
    
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
        ind.aa <- ind.aa * sqrt(var.pheno * h2[3] / var.aa)
        eff.aa <- c(crossprod(ginv(qtn.aa), ind.aa))
        var.aa <- var(ind.aa)
      }
      if (var.ad != 0) {
        ind.ad <- ind.ad * sqrt(var.pheno * h2[4] / var.ad)
        eff.ad <- c(crossprod(ginv(qtn.ad), ind.ad))
        var.ad <- var(ind.ad)
      }
      if (var.da != 0) {
        ind.da <- ind.da * sqrt(var.pheno * h2[5] / var.da)
        eff.da <- c(crossprod(ginv(qtn.da), ind.da))
        var.da <- var(ind.da)
      }
      if (var.dd != 0) {
        ind.dd <- ind.dd * sqrt(var.pheno * h2[6] / var.dd)
        eff.dd <- c(crossprod(ginv(qtn.dd), ind.dd))
        var.dd <- var(ind.dd)
      } 
      effs.adj <- get("effs", envir = inner.env)
      logging.log("adjust interaction effects of markers...\n", verbose = verbose)
      effs.adj$eff1$eff.aa <- eff.aa
      effs.adj$eff1$eff.ad <- eff.ad
      effs.adj$eff1$eff.da <- eff.da
      effs.adj$eff1$eff.dd <- eff.dd
      assign("effs", effs.adj, envir = inner.env)
    }
    
    info.eff$ind.aa <- ind.aa
    info.eff$ind.ad <- ind.ad
    info.eff$ind.da <- ind.da
    info.eff$ind.dd <- ind.dd
    if (verbose) {
      logging.log("Total additiveXadditive   variance:", var.aa,  "\n", verbose = verbose)
      logging.log("Total dominanceXadditive  variance:", var.ad,  "\n", verbose = verbose)
      logging.log("Total dominanceXadditive  variance:", var.da,  "\n", verbose = verbose)
      logging.log("Total dominanceXdominance variance:", var.dd,  "\n", verbose = verbose)
    }
  }
  
  return(info.eff)
}

#' Calculate for fixed effects and random effects
#'
#' Build date: Nov 1, 2019
#' Last update: Nov 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop population information
#' @param FR list of fixed effects, random effects, and their combination
#' @param var.pheno phenotype variances of all traits
#' @param pop.env R environment of population information
#' @param verbose whether to print detail
#'
#' @return list of fixed effects and random effects
#' @export
#'
#' @examples
#' pop <- getpop(100, 1, 0.5)
#' a <- sample(c("a1", "a2", "a3"), 100, replace = TRUE)
#' b <- sample(c("b1", "b2", "b3"), 100, replace = TRUE)
#' pop$a <- a # load your fixed  effects
#' pop$b <- b # load your random effects
#' pop.env <- environment()
#' 
#' # combination of fixed effects
#' cmb.fix <- list(tr1 = c("mu", "gen", "sex"), # trait 1
#'                 tr2 = c("mu", "diet", "season")) # trait 2
#'             
#' fix <- list( 
#'          mu = list(level = "mu", eff = 2),  				      
#'         gen = FALSE,         
#'         sex = list(level = c("1", "2"), eff = c(0.5, 0.3)), 
#'        diet = list(level = c("d1", "d2", "d3"), eff = c(0.1, 0.2, 0.3)),
#'      season = list(level = c("s1", "s2", "s3", "s4"), eff = c(0.1, 0.2, 0.3, 0.2)), 
#'           a = list(level = c("a1", "a2", "a3"), eff = c(0.1, 0.2, 0.3))) 
#' 
#' # combination and ralation of random effects
#' tr1 <- list(rn = c("sir", "PE"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))
#' tr2 <- list(rn = c("dam", "litter"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))          
#' cmb.rand <- list(tr1 = tr1, tr2 = tr2)   
#'       
#' rand <- list(
#'         sir = list(mean = 0, sd = 1),		      
#'         dam = list(mean = 0, sd = 1),	       
#'          PE = list(level = c("p1", "p2", "p3"), eff = c(0.01, 0.02, 0.03)), 
#'      litter = list(level = c("l1", "l2"), eff = c(0.01, 0.02)),    
#'           b = list(level = c("b1", "b2", "b3"), eff = c(0.01, 0.02, 0.03)))
#'   
#' FR <- list(cmb.fix = cmb.fix, fix = fix, cmb.rand = cmb.rand, rand = rand)
#' fr <- cal.FR(pop = pop, FR = FR, var.pheno = c(10, 10), pop.env = pop.env, verbose = TRUE)
#' str(fr) 	          	         	          	         
cal.FR <- function(pop = NULL, FR, var.pheno = NULL, pop.env = NULL, verbose = TRUE) {
  
  nind <- nrow(pop)
  len.tr <- length(var.pheno)
  pop[is.na(pop)] <- "0"
  cmb.fix <- FR$cmb.fix
  fix <- FR$fix
  cmb.rand <- FR$cmb.rand
  rand <- FR$rand
  
  # generate fixed effects of individuals
  fn <- names(fix)
  if (is.null(fn)) {
    fes <- NULL
    
  } else {
    fes <- lapply(1:length(fn), function(i) { return(rep(0, nind)) })
    fes <- as.data.frame(fes)
    names(fes) <- fn
    for (i in 1:length(fix)) { # for
      if (!is.list(fix[[i]])) next
      lev.fix <- as.character(fix[[i]]$level)
      len.fix <- length(lev.fix)
      eff.fix <- fix[[i]]$eff
      if (fn[i] %in% names(pop)) {
        pt <- pop[, names(pop) == fn[i]]
        fe <- rep(0, nind)
        ele.pt <- as.character(unique(pt))
        if (len.fix != length(ele.pt)) { 
          stop("level length should be equal to effects length!")
        } else if (!setequal(lev.fix, ele.pt)) {
          stop(fn[i], " in fix should be corresponding to ", fn[i], " in pop!")
        } 
        for (j in 1:len.fix) {
          fe[pt == lev.fix[j]] <- eff.fix[j]
        }
        
      } else {
        sam <- sample(1:len.fix, nind, replace = TRUE, prob = rep(0.2, len.fix))
        while (length(unique(sam)) != len.fix) {
          sam <- sample(1:len.fix, nind, replace = TRUE, prob = rep(0.2, len.fix))
        }
        fl <- lev.fix[sam]
        pop.adj <- get("pop", envir = pop.env)
        logging.log("add", fn[i], "to population...\n", verbose = verbose)
        pop.adj <- cbind(pop.adj, fl)
        names(pop.adj)[names(pop.adj) == "fl"] <- fn[i]
        assign("pop", pop.adj, envir = pop.env)
        fe <- eff.fix[sam]
      }
 
      fes[[i]] <- fe
    } # end for
  }
  
  # combine fixed effects
  fix <- lapply(1:len.tr, function(i) { return(fes[cmb.fix[[i]]]) })
  names(fix) <- names(cmb.fix[1:len.tr])
  
  # generate random effects of individuals
  rn <- names(rand)
  if (is.null(rn)) {
    res <- NULL
    
  } else {
    res <- lapply(1:length(rn), function(i) { return(rep(0, nind)) })
    res <- as.data.frame(res)
    names(res) <- rn
    for (i in 1:length(rand)) { # for
      if (!is.list(rand[[i]])) next
      if (rn[i] == "sir" | rn[i] == "dam") {
        pt <- pop[, names(pop) == rn[i]]
        re <- rep(0, nind)
        ele.pt <- as.character(unique(pt))
        lev.rand <- as.character(unique(pt))
        len.rand <- length(lev.rand)
        if (len.rand == 1) {
          if (lev.rand == "0") {
            re <- rep(0, nind)
          } else {
            eff.rand <- rnorm(len.rand, mean = rand[[i]]$mean, sd = rand[[i]]$sd)
            for (j in 1:len.rand) {
              re[pt == lev.rand[j]] <- eff.rand[j]
            }
          }
          
        } else {
          eff.rand <- rnorm(len.rand, mean = rand[[i]]$mean, sd = rand[[i]]$sd)
          for (j in 1:len.rand) {
            re[pt == lev.rand[j]] <- eff.rand[j]
          }
        }
        
      } else {
        lev.rand <- as.character(rand[[i]]$level)
        len.rand <- length(lev.rand)
        eff.rand <- rand[[i]]$eff
        
        if (rn[i] %in% names(pop)) {
          pt <- pop[, names(pop) == rn[i]]
          re <- rep(0, nind)
          ele.pt <- as.character(unique(pt))
          if (len.rand != length(ele.pt)) { 
            stop("level length should be equal to effects length!")
          } else if (!setequal(lev.rand, ele.pt)) {
            stop(rn[i], " in rand should be corresponding to ", rn[i], " in pop!")
          } 
          for (j in 1:len.rand) {
            re[pt == lev.rand[j]] <- eff.rand[j]
          }
        
        } else {
          sam <- sample(1:len.rand, nind, replace = TRUE, prob = rep(0.2, len.rand))
          while (length(unique(sam)) != len.rand) {
            sam <- sample(1:len.rand, nind, replace = TRUE, prob = rep(0.2, len.rand))
          }
          rl <- lev.rand[sam]
          pop.adj <- get("pop", envir = pop.env)
          logging.log("add", rn[i], "to population...\n", verbose = verbose)
          pop.adj <- cbind(pop.adj, rl)
          names(pop.adj)[names(pop.adj) == "rl"] <- rn[i]
          assign("pop", pop.adj, envir = pop.env)
          re <- eff.rand[sam]
        }
      }
      res[[i]] <- re
    } # end for
  }
  
  # combine random effects
  rand <- lapply(1:len.tr, function(i) { 
    rn <- cmb.rand[[i]]$rn
    if (is.null(rn)) return(res[NULL])
    ratio <- cmb.rand[[i]]$ratio
    if (length(rn) != length(ratio))
      stop("phenotype variance ratio of random effects should be corresponding to random effects names!")
    rt <- res[rn]

    if (is.null(cmb.rand[[i]]$cr)) {
      cr <- diag(rep(1, length(rn)))
    } else {
      cr <- cmb.rand[[i]]$cr
      if (nrow(cr) != length(rn) | ncol(cr) != length(rn))
        stop("random correlation matrix should be corresponding to random effects names!")
    }
    if (length(ratio) == 1) ratio <- as.matrix(ratio)
    sd <- diag(sqrt(ratio * var.pheno[i]))
    Sigma <- sd %*% cr %*% sd
    mu <- colMeans(rt)
    # adjust random effects
    rt <- build.cov(rt, mu = mu, Sigma = Sigma, tol = 1e-06)
     
    var.r <- apply(rt, 2, var)
    logging.log("The variance of", names(rt), "of", paste0(names(cmb.rand)[i], ":"), var.r, "\n", verbose = verbose)
    return(rt)
  })
  names(rand) <- names(cmb.rand[1:len.tr])
  
  fr <- lapply(1:len.tr, function(i) { return(list(fix = fix[[i]], rand = rand[[i]])) })
  names(fr) <- names(cmb.fix[1:len.tr])
  return(fr)
}

#' Calculate for phenotype
#'
#' Build date: Nov 4, 2019
#' Last update: Nov 4, 2019
#'
#' @author Dong Yin
#'
#' @param fr list of fixed effects and random effects
#' @param info.eff list of phenotype decomposition
#' @param h2 heritability vector of the trait, every elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param num.ind population size
#' @param var.pheno phenotype variace
#' @param verbose whether to print detail
#'
#' @return list of phenotype
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' geno <- geno.cvt(pop.geno)
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
#' h2 <- c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01)             
#' info.eff <- cal.gnt(geno = geno, h2 = h2, effs = effs, 
#'     sel.on = TRUE, inner.env = NULL, verbose = TRUE)
#' 
#' # calculate for phenotype variance
#' ind.a <- info.eff$ind.a
#' var.pheno <- var(ind.a) / h2[1]
#' num.ind <- nrow(pop)
#'         
#' a <- sample(c("a1", "a2", "a3"), 100, replace = TRUE)
#' b <- sample(c("b1", "b2", "b3"), 100, replace = TRUE)
#' pop$a <- a # load your fixed  effects
#' pop$b <- b # load your random effects
#' pop.env <- environment()
#' 
#' # combination of fixed effects
#' cmb.fix <- list(tr1 = c("mu", "gen", "sex"), # trait 1
#'                 tr2 = c("mu", "diet", "season")) # trait 2
#'             
#' fix <- list( 
#'          mu = list(level = "mu", eff = 2),  				      
#'         gen = FALSE,         
#'         sex = list(level = c("1", "2"), eff = c(0.5, 0.3)), 
#'        diet = list(level = c("d1", "d2", "d3"), eff = c(0.1, 0.2, 0.3)),
#'      season = list(level = c("s1", "s2", "s3", "s4"), eff = c(0.1, 0.2, 0.3, 0.2)), 
#'           a = list(level = c("a1", "a2", "a3"), eff = c(0.1, 0.2, 0.3))) 
#' 
#' # combination and ralation of random effects
#' tr1 <- list(rn = c("sir", "PE"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))
#' tr2 <- list(rn = c("dam", "litter"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))          
#' cmb.rand <- list(tr1 = tr1, tr2 = tr2)   
#'       
#' rand <- list(
#'         sir = list(mean = 0, sd = 1),		      
#'         dam = list(mean = 0, sd = 1),	       
#'          PE = list(level = c("p1", "p2", "p3"), eff = c(0.01, 0.02, 0.03)), 
#'      litter = list(level = c("l1", "l2"), eff = c(0.01, 0.02)),    
#'           b = list(level = c("b1", "b2", "b3"), eff = c(0.01, 0.02, 0.03)))
#'   
#' FR <- list(cmb.fix = cmb.fix, fix = fix, cmb.rand = cmb.rand, rand = rand)
#' fr <- cal.FR(pop = pop, FR = FR, var.pheno = var.pheno, pop.env = pop.env, verbose = TRUE)
#' 
#' pheno.list <- cal.pheno(fr = fr, info.eff = info.eff, h2 = h2,
#'     num.ind = num.ind, var.pheno = var.pheno, verbose = TRUE)
#' str(pheno.list)
cal.pheno <- function(fr = NULL, info.eff = NULL, h2 = NULL, num.ind = NULL, var.pheno = NULL, verbose = TRUE) {
  if (!is.null(info.eff)) {
    ind.a <- info.eff$ind.a
    TBV <- ind.a
    TGV <- apply(info.eff, 1, sum)
    var.gnt <- apply(info.eff, 2, var)
  } else {
    ind.a <- rep(0, num.ind)
    TBV <- ind.a
    TGV <- ind.a
    var.gnt <- 0
  }
  
  if (!is.null(fr)) {
    info.eff <- cbind(fr[[1]], info.eff)
    var.fr <- apply(fr[[1]]$rand, 2, var)
  } else {
    var.fr <- 0
  }
  
  var.env <- var.pheno - sum(var.fr, var.gnt)

  if (var.env <= 0) 
    stop("please reduce your fixed variance, random variance or genetic variance to get a positive environmental variance!")
  
  ind.env <- rnorm(num.ind, 0, 1)
  ind.env <- ind.env * sqrt(var.env / var(ind.env))
  info.eff$ind.env <- ind.env
  if (!is.data.frame(info.eff)) info.eff <- as.data.frame(info.eff) 
 
  # get phenotype
  ind.pheno <- apply(info.eff, 1, sum)
  
  Vg <- var.gnt
  Ve <- var(ind.env)
  h2.new <- Vg / sum(Vg, var.fr, Ve)
  logging.log("Total environment         variance:", Ve, "\n", verbose = verbose)
  logging.log("Heritability:", h2.new, "\n", verbose = verbose)
  info.tr <- list(Vg = Vg, Ve = Ve, h2 = h2.new)
  info.pheno <- data.frame(TBV = TBV, TGV = TGV, pheno = ind.pheno)
  pheno.list <- list(info.tr = info.tr, info.eff = info.eff, info.pheno = info.pheno)
  return(pheno.list)  
}

#' Calculate for genetic effects list of selected markers
#'
#' Build date: Sep 11, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param pop.geno genotype of population, a individual has two columns
#' @param cal.model phenotype models with the options: "A", "AD", "ADI"
#' @param num.qtn.tr1 integer or integer vector, the number of QTN in a single trait
#' @param sd.tr1 standard deviation of different effects, the last 5 vector elements are corresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a
#' @param dist.qtn.tr1 distributions of the QTN effects with the options: "normal", "geometry" and "gamma", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param eff.unit.tr1 unit effect of geometric distribution of a single trait, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param shape.tr1 shape of gamma distribution of a single trait, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param scale.tr1 scale of gamma distribution of a single trait, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param multrait whether to apply multiple traits, TRUE represents applying, FALSE represents not
#' @param num.qtn.trn QTN distribution matrix, diagonal elements are total QTN number of the trait, non-diagonal elements are QTN number of overlap QTN between two traits
#' @param sd.trn a matrix with the standard deviation of the QTN effects
#' @param qtn.spot QTN probability in every block
#' @param maf Minor Allele Frequency, marker selection range is from  maf to 0.5
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
#'              num.qtn.tr1 = c(200, 200, 100),
#'              sd.tr1 = c(0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.03),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(400, 100, 100, 400), 2, 2),
#'              sd.trn = matrix(c(0.07, 0, 0, 0.07), 2, 2),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(effs)
cal.effs <-
    function(pop.geno = NULL,
             cal.model = "A",
             num.qtn.tr1 = c(2, 6, 10),
             sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
             dist.qtn.tr1 = rep("normal", 6),
             eff.unit.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
             shape.tr1 = c(1, 1, 1, 1, 1, 1),
             scale.tr1 = c(1, 1, 1, 1, 1, 1),
             multrait = FALSE,
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = matrix(c(1, 0, 0, 0.5), 2, 2),
             qtn.spot = rep(0.1, 10),
             maf = 0, 
             verbose = TRUE) {

# Start calculation

  if (is.null(pop.geno)) return(multrait)
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
    if (any(dim(num.qtn.trn) != dim(sd.trn))) 
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
      effs[[2*i]] <- list(eff.a = rnorm(length(effs[[2*i-1]]), 0, sd.trn[i, i]))
      logging.log("number of selected markers of trait", i, ":", length(effs[[2*i-1]]), "\n", verbose = verbose)
    }

  } else {
    num.qtn <- sum(num.qtn.tr1)
    sel.marker <- sort(sample(1:num.marker, num.qtn, replace = FALSE, prob = wt.marker))
    len.qtn <- length(num.qtn.tr1)
    logging.log("number of selected markers of trait 1:", num.qtn.tr1, "\n", verbose = verbose)
    
    if (cal.model == "A") {
      logging.log("Apply A model...\n", verbose = verbose)
      if (length(sd.tr1) < length(num.qtn.tr1))
        stop("The length of sd.tr1 should be no less than length of num.qtn.tr1!")
      eff.a <- cal.eff(num.qtn.tr1, sd.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff1 <- list(eff.a=eff.a)
      
    } else if (cal.model == "AD") {
      logging.log("Apply AD model...\n", verbose = verbose)
      eff.a <- cal.eff(num.qtn.tr1, sd.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff.d <- cal.eff(sum(num.qtn.tr1), sd.tr1[len.qtn+1], dist.qtn.tr1[2], eff.unit.tr1[2], shape.tr1[2], scale.tr1[2])
      eff1 <- list(eff.a=eff.a, eff.d=eff.d)
      
    } else if (cal.model == "ADI") {
      logging.log("Apply ADI model...\n", verbose = verbose)
      if (num.qtn %% 2 != 0) stop("the number of qtn should be even in the ADI model!")
      # the first part of qtn
      ophalf <- num.qtn %/% 2
      eff.a  <- cal.eff(num.qtn.tr1,      sd.tr1[1:len.qtn], dist.qtn.tr1[1], eff.unit.tr1[1], shape.tr1[1], scale.tr1[1])
      eff.d  <- cal.eff(sum(num.qtn.tr1), sd.tr1[len.qtn+1], dist.qtn.tr1[2], eff.unit.tr1[2], shape.tr1[2], scale.tr1[2])
      eff.aa <- cal.eff(ophalf,           sd.tr1[len.qtn+2], dist.qtn.tr1[3], eff.unit.tr1[3], shape.tr1[3], scale.tr1[3])
      eff.ad <- cal.eff(ophalf,           sd.tr1[len.qtn+3], dist.qtn.tr1[4], eff.unit.tr1[4], shape.tr1[4], scale.tr1[4])
      eff.da <- cal.eff(ophalf,           sd.tr1[len.qtn+4], dist.qtn.tr1[5], eff.unit.tr1[5], shape.tr1[5], scale.tr1[5])
      eff.dd <- cal.eff(ophalf,           sd.tr1[len.qtn+5], dist.qtn.tr1[6], eff.unit.tr1[6], shape.tr1[6], scale.tr1[6])
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
#' @param eff.sd standard deviation of different effects
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
#' eff.sd <- c(0.4, 0.2, 0.02) # three variance of qtn group
#' eff <- cal.eff(num.qtn = num.qtn, eff.sd = eff.sd, dist.qtn = "normal",
#'         eff.unit = 0.5, shape = 1, scale = 1)
#' str(eff)
#'
#' num.qtn <- sum(num.qtn)
#' eff.sd <- sum(eff.sd)
#' eff <- cal.eff(num.qtn = num.qtn, eff.sd = eff.sd, dist.qtn = "normal",
#'         eff.unit = 0.5, shape = 1, scale = 1)
#' str(eff)
cal.eff <- function(num.qtn, eff.sd, dist.qtn, eff.unit, shape, scale) {
  if (sum(num.qtn) == 0) return(0)
  # Judge which kind of distribution of QTN
  eff.qtn <- NULL
  if(dist.qtn == "normal") {
    for (nq in 1:length(num.qtn)) {
    	eff.qtn <- c(eff.qtn, rnorm(num.qtn[nq], 0, eff.sd[nq]))
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
  if (is.null(pop.geno)) return(NULL)
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
#' @param df data.frame without correlation
#' @param mu means of the variables 
#' @param Sigma covariance matrix of variables
#' @param tol tolerance (relative to largest variance) for numerical 
#'       lack of positive-definiteness in Sigma.
#'
#' @return data.frame with correlaion
#' @export
#' @references B. D. Ripley (1987) Stochastic Simulation. Wiley. Page 98
#'
#' @examples
#' Sigma <- matrix(c(14, 10, 10, 15), 2, 2)
#' Sigma
#' df <- cbind(rnorm(100), 0)
#' df <- as.data.frame(df)
#' names(df) <- paste0("tr", 1:ncol(df))
#' df.cov <- build.cov(df, Sigma = Sigma)
#' var(df.cov)
build.cov <- function(df = NULL, mu = rep(0, nrow(Sigma)), Sigma, tol = 1e-06) {
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

#' Set the phenotype of the population
#'
#' @param pop population information of generation, family ID, within-family ID, individual ID, paternal ID, maternal ID, and sex
#' @param pop.pheno phenotype information
#' @param sel.crit selection criteria with the options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", and "pheno"
#'
#' @return popluation information with phenotype
#' @export
#'
#' @examples
#' pop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' geno <- geno.cvt(pop.geno)
#' nind <- nrow(pop)
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
#' h2 <- c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01)             
#' info.eff <- cal.gnt(geno = geno, h2 = h2, effs = effs, 
#'     sel.on = TRUE, inner.env = NULL, verbose = TRUE)
#' 
#' # calculate for phenotype variance
#' ind.a <- info.eff$ind.a
#' var.pheno <- var(ind.a) / h2[1]
#'         
#' a <- sample(c("a1", "a2", "a3"), 100, replace = TRUE)
#' b <- sample(c("b1", "b2", "b3"), 100, replace = TRUE)
#' pop$a <- a # load your fixed  effects
#' pop$b <- b # load your random effects
#' pop.env <- environment()
#' 
#' # combination of fixed effects
#' cmb.fix <- list(tr1 = c("mu", "gen", "sex"), # trait 1
#'                 tr2 = c("mu", "diet", "season")) # trait 2
#'             
#' fix <- list( 
#'          mu = list(level = "mu", eff = 2),  				      
#'         gen = FALSE,         
#'         sex = list(level = c("1", "2"), eff = c(0.5, 0.3)), 
#'        diet = list(level = c("d1", "d2", "d3"), eff = c(0.1, 0.2, 0.3)),
#'      season = list(level = c("s1", "s2", "s3", "s4"), eff = c(0.1, 0.2, 0.3, 0.2)), 
#'           a = list(level = c("a1", "a2", "a3"), eff = c(0.1, 0.2, 0.3))) 
#' 
#' # combination and ralation of random effects
#' tr1 <- list(rn = c("sir", "PE"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))
#' tr2 <- list(rn = c("dam", "litter"), ratio = c(0.01, 0.03), 
#'             cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))          
#' cmb.rand <- list(tr1 = tr1, tr2 = tr2)   
#'       
#' rand <- list(
#'         sir = list(mean = 0, sd = 1),		      
#'         dam = list(mean = 0, sd = 1),	       
#'          PE = list(level = c("p1", "p2", "p3"), eff = c(0.01, 0.02, 0.03)), 
#'      litter = list(level = c("l1", "l2"), eff = c(0.01, 0.02)),    
#'           b = list(level = c("b1", "b2", "b3"), eff = c(0.01, 0.02, 0.03)))
#'   
#' FR <- list(cmb.fix = cmb.fix, fix = fix, cmb.rand = cmb.rand, rand = rand)
#' fr <- cal.FR(pop = pop, FR = FR, var.pheno = var.pheno, pop.env = pop.env, verbose = TRUE)
#' 
#' pheno.list <- cal.pheno(fr = fr, info.eff = info.eff, h2 = h2,
#'     num.ind = nind, var.pheno = var.pheno, verbose = TRUE)
#' str(pheno.list)
#' 
#' str(pop)
#' pop <- set.pheno(pop, pheno.list, sel.crit = "pheno")
#' str(pop)
set.pheno <- function(pop, pop.pheno, sel.crit) {
  f1 <- grep(pattern = "TBV|TGV|pheno|ebv", x = names(pop), value = FALSE)
  if (length(f1) != 0) pop <- pop[, -f1]
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
  pop <- cbind(pop, subset(pop.pheno$info.pheno, select = pn))
  return(pop)
}
