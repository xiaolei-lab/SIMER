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


#' Read and compare breeding plans
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 31, 2019
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters
#' @param simls a list with all simer result information
#' @param FR list of fixed effects, random effects, and their combination
#' @param index.wt economic weights of selection index method, its length should equals to the number of traits
#' @param decr whether to sorting with decreasing
#' @param selPath the path of breeding_plan
#' @param verbose whether to print details
#'
#' @return None
#' @export
#'
#' @examples
#' pop.env <- list(
#'   F1 = list(
#'     level = c("1", "2"),
#'     eff = list(tr1 = c(50, 30), tr2 = c(50, 30))
#'   ), 
#'   F2 = list(
#'     level = c("d1", "d2", "d3"),
#'     eff = list(tr1 = c(10, 20, 30), tr2 = c(10, 20, 30))
#'   ),
#'   R1 = list(
#'     level = c("l1", "l2", "l3"),
#'     ratio = list(tr1 = 0.1, tr2 = 0.1)
#'   )
#' )
#' SP <- param.simer(
#'   qtn.num = diag(10, 2),
#'   pop.env = pop.env,
#'   phe.model = list(
#'     tr1 = "T1 = A + F1 + F2 + R1 + E",
#'     tr2 = "T2 = A + F1 + F2 + R1 + E"
#'   )
#' )
#' 
#' SP <- simer(SP, outpath = "/home/yindong/Desktop")
#' 
#' selPath <- system.file("extdata", "04breeding_plan", package = "simer")
#' 
#' # follow code can be run only if you have installd hiblup package
#' # goal.plan <- complan(simls = simer.list, FR = FR, index.wt = c(0.5, 0.5), 
#' #                      decr = TRUE, selPath = selPath, verbose = TRUE)
#' # str(goal.plan)
complan <- function(SP = NULL, simls=NULL, FR=NULL, index.wt=c(0.5, 0.5), decr = TRUE, selPath=NULL, verbose=TRUE) {
  
  if (!dir.exists(selPath)) stop("Please input a right selection path!")
  
  verbose <- SP$global$verbose
  phe.var <- SP$pheno$phe.var
  pop <- do.call(rbind, SP$pheno$pop)
  pop.last <- SP$pheno$pop[[length(SP$pheno$pop)]]
  phe.name <- grep(pattern = "TBV", x = names(pop.last), value = TRUE)
  phe.name <- substr(phe.name, 1, nchar(phe.name) - 4)
  score.last <- colMeans(subset(pop.last, select = phe.name))
  
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
  plan.names <- paste("plan", 1:length(filenames), sep = "")
  
  idx4hi <- lapply(1:length(filenames), function(i) {
    filename <- filenames[i]
    filename.old <- filename
    filename <- file.path(selPath, filename)
    fileImage <- file(description=filename, open="r")
    
    logging.log(" Read breeding plan ", filenames[i], "...", sep = "", verbose = verbose)
    # Initialization for iteration within file
    inFile=TRUE
    while (inFile) {
      tt <- readLines(fileImage, n=1)
      if (length(tt) == 0) {
        inFile <- FALSE
        next
      }
      if (tt == '# start genotyping') {
        geno_tab <- NULL
        ttt <- readLines(fileImage, n=1)
        cn <- unlist(strsplit(ttt, split = '\t'))
        while (ttt != '') {
          ttt <- readLines(fileImage, n=1)
          if (ttt == '') {
            geno_tab <- as.data.frame(geno_tab, stringsAsFactors = FALSE)
            names(geno_tab) <- cn
            break
          }
          geno_tab <- rbind(geno_tab, unlist(strsplit(ttt, split = '\t')))
        }
        
      } else if (tt == '# start phenotyping') {
        pheno_tab <- NULL
        ttt <- readLines(fileImage, n=1)
        cn <- unlist(strsplit(ttt, split = '\t'))
        while (ttt != '') {
          ttt <- readLines(fileImage, n=1)
          if (ttt == '') {
            pheno_tab <- as.data.frame(pheno_tab, stringsAsFactors = FALSE)
            names(pheno_tab) <- cn
            break
          }
          pheno_tab <- rbind(pheno_tab, unlist(strsplit(ttt, split = '\t')))
        } # end inner while
        
      } else if (tt == "# fixed effect") {
        ttt <- readLines(fileImage, n=1)
        eff_fixed <- list()
        fj <- 0
        if (length(ttt) == 0) {
          inFile <- FALSE
          next
        }
        while (ttt != '') {
          fj <- fj + 1
          eff_fixed[[fj]] <- unlist(strsplit(ttt, split = "\t"))
          ttt <- readLines(fileImage, n=1)
          if (length(ttt) == 0) {
            inFile <- FALSE
            break
          }
          if (ttt == '') break
        }
        if (fj > 0) names(eff_fixed) <- paste('tr', 1:fj, sep = '')
        
      } else if (tt == "# random effect") {
        ttt <- readLines(fileImage, n=1)
        eff_random <- list()
        rj <- 0
        if (length(ttt) == 0) {
          inFile <- FALSE
          next
        }
        while (ttt != '') {
          rj <- rj + 1
          eff_random[[rj]] <- unlist(strsplit(ttt, split = "\t"))
          ttt <- readLines(fileImage, n=1)
          if (length(ttt) == 0) {
            inFile <- FALSE
            break
          }
          if (ttt == '') break
          
        }
        if (rj > 0) names(eff_random) <- paste('tr', 1:rj, sep = '')
      }
    } # end outer while
    close.connection(fileImage)
    logging.log("Done!\n", verbose = verbose)
    
    idx_geno <- sel.idx(geno_tab, pop)
    idx_pheno <- sel.idx(pheno_tab, pop)
    
    out.index <- list(idx_geno = idx_geno, idx_pheno = idx_pheno, eff_fixed = eff_fixed, eff_random = eff_random)
    return(out.index)
  })
  names(idx4hi) <- plan.names
  
  # extract fixed effects and random effects
  ls.fr <- lapply(idx4hi, function(idx) {
    fac_fixed <- NULL
    fac_random <- NULL
    if (length(idx$eff_fixed) != 0) {
      ntr <- length(idx$eff_fixed)
      trn <- paste0("tr", 1:ntr)
      for (i in 1:ntr) {
        fac_fixed[[i]] <- pop[, idx$eff_fixed[[i]], drop = FALSE]
      }
      names(fac_fixed) <- trn
    }
    if (length(idx$eff_random) != 0) {
      ntr <- length(idx$eff_random)
      trn <- paste0("tr", 1:ntr)
      for (i in 1:ntr) {
        fac_random[[i]] <- pop[, idx$eff_random[[i]], drop = FALSE]
      }
      names(fac_random) <- trn
    }
    ls_fr <- list(fac_fixed = fac_fixed, fac_random = fac_random)
    return(ls_fr)
  })
  names(ls.fr) <- plan.names
  
  # prepare for calculating ebv
  pn <- grep(pattern = "pheno", x = names(pop), value = TRUE)
  pheno_full <- subset(pop, select = c("index", pn))
  pedigree <- subset(pop, select = c("index", "sir", "dam"))
  map <- simls$map[, 1:3]
  mode <- ifelse("ind.d" %in% names(simls$trait$info.eff), "AD", "A")
  idx.ebv <- lapply(1:length(idx4hi), function(i) {
    logging.log(" \n", verbose = verbose)
    logging.log(" Call HIBLUP for analysing breeding plan", filenames[i], "\n", verbose = verbose)
    geno.id <- idx4hi[[i]]$idx_geno
    geno4hi <- geno.cvt1(simls$geno)
    geno <- geno4hi[, simls$genoid %in% geno.id]
    geno.id <- as.data.frame(geno.id)
    pn <- grep(pattern = "pheno", x = names(pop), value = TRUE)
    pheno <- pheno_full[idx4hi[[i]]$idx_pheno, ]
    fcf <- ls.fr[[i]]$fac_fixed
    fcr <- ls.fr[[i]]$fac_random
    CV <- NULL
    R <- NULL
    bivar.CV <- NULL
    bivar.R <- NULL
    if (ncol(pheno) != 2) {
      bivar.pos <- 2:ncol(pheno)
    } else {
      bivar.pos <- NULL
    }
    if (!is.null(fcf)) {
      CV.t <- lapply(1:length(fcf), function(jf) {
        cv.t <- model.matrix(~., data = fcf[[jf]])[idx4hi[[i]]$idx_pheno, ]
        if (!is.data.frame(cv.t)) cv.t <- as.data.frame(cv.t)
        return(cv.t)
      })
      if (ncol(pheno) == 2) {
        CV <- CV.t[[1]]
      } else {
        bivar.CV <- CV.t
      }
    }
    if (!is.null(fcr)) {
      R.t <- lapply(1:length(fcr), function(jr) {
        if (!is.data.frame(fcr[[jr]])) fcr[[jr]] <- as.data.frame(fcr[[jr]])
        r.t <- fcr[[jr]][idx4hi[[i]]$idx_pheno, ]
        if (!is.data.frame(r.t)) r.t <- as.data.frame(r.t)
        return(r.t)
      })
      if (ncol(pheno) == 2) {
        R <- R.t[[1]]
      } else {
        bivar.R <- R.t
      }
    }
    gebv <- NULL
    eval(parse(text = "tryCatch({
      if (!(\"hiblup\" %in% .packages())) suppressMessages(library(hiblup))
      gebv <- hiblup(pheno = pheno, bivar.pos = bivar.pos, geno = geno, map = map,
                     geno.id = geno.id, file.output = FALSE, pedigree = pedigree, mode = mode,
                     CV = CV, R = R, bivar.CV = bivar.CV, bivar.R = bivar.R, snp.solution = FALSE)
    }, error=function(e) {
      stop(\"Something wrong when running HIBLUP!\") })"))
    rm(geno4hi); rm(geno); gc()
    return(gebv$ebv[gebv$ebv[, 1] %in% pop.last$index, ])
  })
  names(idx.ebv) <- plan.names
  logging.log(" \n", verbose = verbose)

  # mating for the next generation
  ps <- rep(0.8, 2)
  num.prog <- 4
  ratio <- 0.5
  # extract genotype
  geno.ed <- pop.last$index[length(pop.last$index)] *2
  geno.op <- geno.ed - 2*nrow(pop.last) + 1
  pop.geno.last <- simls$geno[, geno.op:geno.ed]
  logging.log(" Mating according to analysises...", verbose = verbose)
  pop.gp <- lapply(idx.ebv, function(ebv) {
    ebv.sum <- ebv
    if (ncol(ebv.sum) != 2) ebv.sum <- apply(ebv.sum[, 2:ncol(ebv.sum)], 1, sum)
    ind.score.ordered <- ebv[order(ebv.sum, decreasing=TRUE), 1]
    count.sir <- sum(pop.last$sex == 1 | pop.last$sex == 0) * ps[1]
    count.dam <- sum(pop.last$sex == 2 | pop.last$sex == 0) * ps[2]
    ind.stay <- getsd(ind.score.ordered, pop.last, count.sir[1], count.dam[1])
    
    pop.gp.in <- # pop.gp with genotype and pop information
      reproduces(pop1 = pop.last,
                 pop1.geno.id = pop.last$index, 
                 pop1.geno = pop.geno.last,
                 ind.stay = ind.stay,
                 mtd.reprod = "randmate",
                 num.prog = 4,
                 ratio = 0.5)

    return(pop.gp.in)
  })
  names(pop.gp) <- plan.names
  logging.log("Done!\n", verbose = verbose)

  # calculate for phenotype
  effs <- simls$effs
  h2.tr1 <- 0.3
  sel.on <- TRUE
  inner.env <- environment()
  verbose <- TRUE
  if (length(simls$trait$info.tr$h2) == 1) {
    h2.tr1 <- simls$trait$info.tr$h2
    h2.trn <- NULL
    gnt.cov <- NULL
  } else {
    h2.tr1 <- NULL
    h2.trn <- simls$trait$info.tr$h2
    gnt.cov <- simls$trait$info.tr$Covg
  }
  pheno.curr <- lapply(1:length(pop.gp), function(i) {
    logging.log("\n", verbose = verbose)
    logging.log(" Generate phenotypes for breeding plan", filenames[i], "\n", verbose = verbose)
    gp <- pop.gp[[i]]
    pop1.pheno <-
      phenotype(effs = effs,
                FR = FR,
                pop = gp$pop,
                pop.geno = gp$geno,
                pos.map = NULL,
                var.pheno = var.pheno, 
                h2.tr1 = h2.tr1,
                gnt.cov = gnt.cov,
                h2.trn = h2.trn,
                sel.crit = "pheno",
                pop.total = rbind(pop[, 1:7], gp$pop),
                sel.on = TRUE,
                inner.env =  inner.env,
                verbose = verbose)
    pop.curr <- pop1.pheno$pop
    pheno <- subset(pop.curr, select = f1)

    if (ncol(pheno) > 1) {
      # calculate the weigth of index selection
      b <- cal.idx(pop.pheno = pop1.pheno, index.wt = index.wt)
    } else {
      b <- 1
    }
    pheno <- rowSums(sweep(pheno, 2, b, "*"))
    pheno.curr <- list(pheno = pheno, b = b)

    return(pheno.curr)
  })
  names(pheno.curr) <- plan.names

  score.curr <- unlist(lapply(1:length(pheno.curr), function(i) {
    pheno <- mean(pheno.curr[[i]]$pheno)
    score.curr <- pheno - sum(score.last * pheno.curr[[i]]$b)
    return(score.curr)
  }))
  names(score.curr) <- plan.names
  score.max <- max(score.curr)
  score.min <- min(score.curr)

  logging.log("\n", verbose = verbose)
  logging.log(" Genetic progress of every breeding plan:\n", verbose = verbose)
  for (i in 1:length(score.curr)) {
    if (score.curr[i] == score.max) {
      str0 <- " (max)"
    } else if (score.curr[i] == score.min) {
      str0 <- " (min)"
    } else {
      str0 <- ""
    }
    logging.log(" Plan", i, ": ", score.curr[i], str0, "\n", sep = "", verbose = verbose)
  }
  logging.log("\n", verbose = verbose)

  if (decr) {
    idx.goal <- which.max(score.curr)
  } else {
    idx.goal <- which.min(score.curr)
  }
  plan.goal <- list(plan_name = plan.names[idx.goal], plan_cont = idx4hi[[idx.goal]], plan_gains = score.curr)
  
  logging.log(" Individuals should be genotyping are:\n", verbose = verbose)
  logging.print(plan.goal$plan_cont$idx_geno, verbose = verbose)
  logging.log(" Individuals should be phenotyping are:\n", verbose = verbose)
  logging.print(plan.goal$plan_cont$idx_pheno, verbose = verbose)

  if (length(plan.goal$plan_cont$eff_fixed) > 0) {
    logging.log(" Fixed effects should be in the model are:\n", verbose = verbose)
    for (i in 1:length(plan.goal$plan_cont$eff_fixed)) {
      logging.log(" Trait", i, ":", plan.goal$plan_cont$eff_fixed[[i]], "\n", verbose = verbose)
    }
  }
  if (length(plan.goal$plan_cont$eff_random) > 0) {
    logging.log(" Random effects should be in the model are:\n", verbose = verbose)
    for (i in 1:length(plan.goal$plan_cont$eff_random)) {
      logging.log(" Trait", i, ":", plan.goal$plan_cont$eff_random[[i]], "\n", verbose = verbose)
    }
  }
  if (length(pheno.curr[[idx.goal]]$b) != 1) {
    plan.goal$b <-  pheno.curr[[idx.goal]]$b
    logging.log(" The weights of index selection are:\n", verbose = verbose)
    logging.log("", plan.goal$b, "\n", verbose = verbose)
  }

  rm(simls); rm(FR); rm(pop); rm(pop.last); rm(idx4hi); rm(ls.fr); rm(idx.ebv); 
  rm(pop.geno.last); rm(pop.gp); rm(pheno.curr); gc()
  return(plan.goal)
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

#' Select individual ID according to scheme
#'
#' Build date: Jan 13, 2020
#' Last update: Jan 13, 2020
#'
#' @author Dong Yin
#'
#' @param scheme breeding scheme for genotyping or phenotyping
#' @param pop population information 
#'
#' @return index vector
#' @export
#'
#' @examples
#' scheme <- data.frame(
#'     generation = c("1", "2", "3"), 
#'     family_index = rep("1:3,7",3),
#'     within_family_index = rep("all", 3), 
#'     sex = rep("1", 3), 
#'     stringsAsFactors = FALSE
#' )
#' pop <- getpop(nind = 100, from = 1, ratio = 0.5)
#' pop$gen <- rep(1:5, each = 20)
#' out.idx <- sel.idx(scheme = scheme, pop = pop)
sel.idx <- function(scheme, pop) {
  pop_sel <- NULL
  for(i in 1:nrow(scheme)) { # for2
    gen <- scheme$generation[i]
    pop_gen <- pop[pop$gen == gen, ]
    
    fam.index <- scheme$family_index[i]
    if (fam.index == "all") {
      pop_fam <- pop_gen
    } else {
      fam.index <- num2num(fam.index) + pop_gen$index[1] - 1
      pop_fam <- pop_gen[pop_gen$fam %in% fam.index, ]
    }
    
    infam.index <- scheme$within_family_index[i]
    if (infam.index == "all") {
      pop_infam <- pop_fam
    } else {
      infam.index <- num2num(infam.index)
      pop_infam <- pop_fam[pop_fam %in% infam.index, ]
    }
    
    sex.index <- scheme$sex[i]
    if (sex.index == "all") {
      pop_sex <- pop_infam
    } else {
      pop_sex <- pop_infam[pop_infam$sex == as.numeric(scheme$sex[i]), ]
    }
    
    pop_sel <- rbind(pop_sel, pop_sex)
  } # end for2
  if (nrow(pop_sel) == 0) stop("No individual left!")
  return(pop_sel$index)
}
