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


#' Main function of simer
#'
#' Build date: Jan 7, 2019
#' Last update: Oct 13, 2019
#'
#' @author Dong Yin, Lilin Yin, Haohao Zhang and Xiaolei Liu
#'
#' @param num.gen number of generations in simulation
#' @param replication replication index of simulation
#' @param verbose whether to print detail
#' @param mrk.dense whether markers are dense, it is TRUE when sequencing data
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param outcols the column number of an individual in the output genotype matrix, it can be 1 or 2 
#' @param out prefix of output file name
#' @param outpath path of output files
#' @param selPath the path of breeding_plan
#' @param out.format format of output, "numeric" or "plink"
#' @param seed.sim random seed of a simulation process
#' @param out.geno.gen indice of generations of output genotype
#' @param out.pheno.gen indice of generations of output phenotype
#' @param rawgeno1 extrinsic genotype matrix1
#' @param rawgeno2 extrinsic genotype matrix2
#' @param rawgeno3 extrinsic genotype matrix3
#' @param rawgeno4 extrinsic genotype matrix4
#' @param num.ind population size of the base population
#' @param prob weight of "0" and "1" in genotype matrix, the sum of elements in vector equal to 1
#' @param input.map map that should be input, the marker number should be consistent in both map file and genotype data
#' @param len.block length of every blocks
#' @param range.hot range of number of chromosome crossovers in a hot spot block
#' @param range.cold range of number of chromosome crossovers in a cold spot block
#' @param rate.mut mutation rate between 1e-8 and 1e-6
#' @param cal.model phenotype model with the options: "A", "AD", "ADI"
#' @param FR list of fixed effects, random effects, and their combination
#' @param cv list of population Coefficient of Variation or family Coefficient of Variation
#' @param var.pheno the phenotype variance, only used in single-trait simulation
#' @param h2.tr1 heritability vector of a single trait, every element are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param num.qtn.tr1 integer or integer vector, the number of QTN in a single trait
#' @param sd.tr1 standard deviation of different effects, the last 5 vector elements are corresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a
#' @param dist.qtn.tr1 distributions of the QTN effects with the options: "normal", "geometry", "gamma", and "beta", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively
#' @param prob.tr1 unit effect of geometric distribution of a single trait, its length should be same as dist.qtn.tr1
#' @param shape.tr1 shape of gamma distribution of a single trait, its length should be same as dist.qtn.tr1
#' @param scale.tr1 scale of gamma distribution of a single trait, its length should be same as dist.qtn.tr1
#' @param shape1.tr1 non-negative parameters of the Beta distribution, its length should be same as dist.qtn.tr1
#' @param shape2.tr1 non-negative parameters of the Beta distribution, its length should be same as dist.qtn.tr1
#' @param ncp.tr1 non-centrality parameter, its length should be same as dist.qtn.tr1
#' @param multrait whether to apply multiple traits, TRUE represents applying, FALSE represents not
#' @param num.qtn.trn QTN distribution matrix, diagonal elements are total QTN number of the trait, non-diagonal elements are QTN number of overlap QTN between two traits
#' @param sd.trn a matrix with the standard deviation of the QTN effects
#' @param gnt.cov genetic covariance matrix among all traits
#' @param h2.trn heritability among all traits
#' @param qtn.spot QTN probability in every block
#' @param maf Minor Allele Frequency, marker selection range is from maf to 0.5
#' @param sel.crit selection criteria with the options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", and "pheno"
#' @param sel.on whether to add selection
#' @param mtd.reprod different reproduction methods with the options: "clone", "dh", "selfpol", "singcro", "tricro", "doubcro", "backcro","randmate", "randexself", and "userped"
#' @param userped user-designed pedigree to control mating process
#' @param num.prog litter size of dams
#' @param ratio ratio of the males in all individuals
#' @param prog.tri litter size of the first single cross process in trible cross process
#' @param prog.doub litter size of the first two single cross process in double cross process
#' @param prog.back a vector with litter size in every generation of back-cross
#' @param refresh refresh ratio of core population of sires and dams, only used in ps > 1
#' @param keep.max.gen the max keep generation range in the selection for sires and dams, only used in ps > 1
#' @param ps if ps <= 1, fraction selected in selection of males and females; if ps > 1, ps is number of selected males and females
#' @param decr whether to sort by descreasing
#' @param sel.multi selection method of multiple traits with options: "tdm", "indcul" and "index"
#' @param index.wt economic weights of selection index method, its length should equals to the number of traits
#' @param index.tdm index represents which trait is being selected
#' @param goal.perc percentage of goal more than mean of scores of individuals
#' @param pass.perc percentage of expected excellent individuals
#' @param sel.sing selection method of single trait with options: "ind", "fam", "infam" and "comb"
#'
#' @return a list with population information, genotype matrix, map information, selection intensity
#' @export
#' @import bigmemory
#' @importFrom stats aov cor dnorm qnorm rgamma rnorm rbeta rgeom runif var shapiro.test
#' @importFrom utils write.table read.delim packageVersion
#' @importFrom methods getPackageName
#' @importFrom MASS mvrnorm ginv
#' @importFrom rMVP MVP.Data.MVP2Bfile
#'
#' @examples
#' \donttest{
#' # get map file, map is necessary
#' data(simdata)
#'
#' # run simer
#' simer.list <-
#'      simer(num.gen = 5,
#'            replication = 1,
#'            verbose = TRUE, 
#'            mrk.dense = TRUE,
#'            incols = 2, 
#'            outcols = 1, 
#'            out = "simer", 
#'            outpath = NULL,
#'            selPath = NULL, 
#'            out.format = "numeric",
#'            seed.sim = runif(1, 0, 100),
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
#'            cal.model = "A",
#'            FR = NULL, 
#'            cv = NULL, 
#'            h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'            num.qtn.tr1 = 500,
#'            sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02),
#'            dist.qtn.tr1 = rep("normal", 6),
#'            prob.tr1 = rep(0.5, 6),
#'            shape.tr1 = rep(1, 6),
#'            scale.tr1 = rep(1, 6),
#'            shape1.tr1 = rep(1, 6),
#'            shape2.tr1 = rep(1, 6),
#'            ncp.tr1 = rep(0, 6),
#'            multrait = FALSE,
#'            num.qtn.trn = matrix(c(400, 100, 100, 400), 2, 2),
#'            sd.trn = matrix(c(0.07, 0, 0, 0.07), 2, 2),
#'            gnt.cov = matrix(c(1, 2, 2, 16), 2, 2),
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
#'            ps = rep(0.8, 2),
#'            decr = TRUE,
#'            sel.multi = "index",
#'            index.wt = c(0.5, 0.5),
#'            index.tdm = 1,
#'            goal.perc = 0.1,
#'            pass.perc = 0.9, 
#'            sel.sing = "comb") 
#' pop <- simer.list$pop
#' effs <- simer.list$effs
#' trait <- simer.list$trait
#' geno <- simer.list$geno
#' genoid <- simer.list$genoid
#' map <- simer.list$map
#' si <- simer.list$si
#' head(pop)
#' str(effs)
#' str(trait)
#' geno[1:6, 1:6]
#' genoid[1:6]
#' str(map)
#' si
#' }
simer <-
    function(num.gen = 5,
             replication = 1,
             verbose = TRUE, 
             mrk.dense = TRUE,
             incols = 2, 
             outcols = 1, 
             out = "simer", 
             outpath = NULL,
             selPath = NULL, 
             out.format = "numeric",
             seed.sim = runif(1, 0, 100),
             out.geno.gen = (num.gen-2):num.gen,
             out.pheno.gen = 1:num.gen,
             rawgeno1 = NULL,
             rawgeno2 = NULL,
             rawgeno3 = NULL,
             rawgeno4 = NULL,
             num.ind = 100,
             prob = c(0.5, 0.5),
             input.map = NULL,
             len.block = 5e7,
             range.hot = 4:6,
             range.cold = 1:5,
             rate.mut = 1e-8,
             cal.model = "A",
             FR = NULL, 
             cv = NULL, 
             var.pheno = NULL, 
             h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
             num.qtn.tr1 = 500,
             sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02),
             dist.qtn.tr1 = rep("normal", 6),
             prob.tr1 = rep(0.5, 6),
             shape.tr1 = rep(1, 6),
             scale.tr1 = rep(1, 6),
             shape1.tr1 = rep(1, 6), 
             shape2.tr1 = rep(1, 6), 
             ncp.tr1 = rep(0, 6), 
             multrait = FALSE,
             num.qtn.trn = matrix(c(400, 100, 100, 400), 2, 2),
             sd.trn = matrix(c(1, 0, 0, 0.5), 2, 2),
             gnt.cov = matrix(c(1, 2, 2, 16), 2, 2),
             h2.trn = c(0.3, 0.5), 
             qtn.spot = rep(0.1, 10),
             maf = 0,
             sel.crit = "pheno",
             sel.on = TRUE, 
             mtd.reprod = "randmate",
             userped = NULL,
             num.prog = 2,
             ratio = 0.5,
             prog.tri = 2,
             prog.doub = 2,
             prog.back = rep(2, num.gen),
             refresh = c(1, 0.6),
             keep.max.gen = rep(3, 2),
             ps = rep(0.8, 2),
             decr = TRUE,
             sel.multi = "index",
             index.wt = c(0.5, 0.5),
             index.tdm = 1,
             goal.perc = 0.1,
             pass.perc = 0.9,
             sel.sing = "comb") {

# Start simer

# TODO: How to generate inbreeding sirs and uninbreeding dams
# TODO: optcontri.sel
# TODO: add MVP for output
# TODO: correct pedigree     
# TODO: add superior limit of homo   
# TODO: add multiple fix and random effects
# TODO: add summary() to population information
# TODO: add inbreeding coefficient
# TODO: update index selection
# TODO: add true block distribution  
# TODO: genomic mating for future
# TODO: inbreeding change in every generations
  
  simer.Version(width = 70, verbose = verbose)    
  
  inner.env <- environment()    
  # initialize logging
  if (!is.null(outpath)) {
    if (!dir.exists(outpath)) stop(paste0("Please check your output path: ", outpath))
    if (verbose) {
      logging.initialize("Simer", outpath = outpath)
    }
  }
  
	################### MAIN_FUNCTION_SETTING ###################
  logging.log("--------------------------- replication ", replication, "---------------------------\n", verbose = verbose)
  op <- Sys.time()
  logging.log(" SIMER BEGIN AT", as.character(op), "\n", verbose = verbose)
  set.seed(seed.sim)
  if (incols == 1) outcols <- 1

	################### BASE_POPULATION ###################
  # establish genotype of base population if there isn't by two ways:
  # 1. input rawgeno
  # 2. input num.marker and num.ind
  
  nmrk <- nrow(input.map)
  # combine genotype matrix
  if (is.list(rawgeno1)) {
    if (!(mtd.reprod == "randmate" || mtd.reprod == "randexself")) 
      stop("Only random matings support genotype list!")
    nsir <- ncol(rawgeno1$sir) / incols
    ndam <- ncol(rawgeno1$dam) / incols
    nind <- nsir + ndam
    basepop <- getpop(nind, 1, nsir/nind)
    rawgeno1 <- cbind(rawgeno1$sir[], rawgeno1$dam[])
  } else {
    # set base population information
    nind <- ifelse(is.null(rawgeno1), num.ind, ncol(rawgeno1) / incols)
    nsir <- nind * ratio
    ndam <- nind * (1-ratio)
    basepop <- getpop(nind, 1, ratio)
  }
  
  num.marker <- nrow(input.map)
  logging.log(" --- base population 1 ---\n", verbose = verbose)
  basepop.geno <-
      genotype(rawgeno = rawgeno1,
               incols = incols, 
               num.marker = num.marker,
               num.ind = num.ind,
               prob = prob, 
               verbose = verbose)

  # set block information and recombination information
  num.ind <- nind
  pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = len.block)
  blk.rg <- cal.blk(pos.map)
  recom.spot <- as.numeric(pos.map[blk.rg[, 1], 7])

  # calculate for marker information
  effs <-
    cal.effs(pop.geno = basepop.geno,
             incols = incols, 
             cal.model = cal.model,
             num.qtn.tr1 = num.qtn.tr1,
             sd.tr1 = sd.tr1,
             dist.qtn.tr1 = dist.qtn.tr1,
             prob.tr1 = prob.tr1,
             shape.tr1 = shape.tr1,
             scale.tr1 = scale.tr1,
             shape1.tr1 = shape1.tr1, 
             shape2.tr1 = shape2.tr1, 
             ncp.tr1 = ncp.tr1, 
             multrait = multrait,
             num.qtn.trn = num.qtn.trn,
             sd.trn = sd.trn,
             qtn.spot = qtn.spot,
             maf = maf, 
             verbose = verbose)

  # calculate phenotype according to genotype
  if (sel.on) {
    pop1.pheno <-
      phenotype(effs = effs,
                FR = FR, 
                cv = cv, 
                pop = basepop,
                pop.geno = basepop.geno,
                pos.map = pos.map,
                var.pheno = var.pheno, 
                h2.tr1 = h2.tr1,
                gnt.cov = gnt.cov,
                h2.trn = h2.trn, 
                sel.crit = sel.crit, 
                pop.total = basepop, 
                sel.on = sel.on, 
                inner.env =  inner.env, 
                verbose = verbose)
    basepop <- pop1.pheno$pop
    pop1.pheno$pop <- NULL
  }
  
  # only mutation in clone and doubled haploid
  if (mtd.reprod == "clone" || mtd.reprod == "dh" || mtd.reprod == "selfpol") {
    basepop$sex <- 0
    recom.spot <- NULL
    ratio <- 0
  }
  
  basepop.geno.em <-  # genotype matrix after Mutation
    genotype(geno = basepop.geno,
             incols = incols, 
             blk.rg = blk.rg,
             recom.spot = recom.spot,
             range.hot = range.hot,
             range.cold = range.cold,
             rate.mut = rate.mut, 
             verbose = verbose)

  if (mtd.reprod == "singcro" || mtd.reprod == "tricro" || mtd.reprod == "doubcro" || mtd.reprod == "backcro") {
    # set base population information
    basepop$sex <- 1

    if (is.null(rawgeno2)) {
      logging.log(" --- base population 2 ---\n", verbose = verbose)
      prob1 <- runif(1)
      prob <- c(prob1, 1 - prob1)
      pop2.geno <- genotype(incols = incols, num.marker = num.marker, num.ind = num.ind, prob = prob, verbose = verbose)
    } else {
      pop2.geno <- genotype(rawgeno = rawgeno2, verbose = verbose)
    }

    # set base population information
    nind2 <- ncol(pop2.geno) / incols
    pop2 <- getpop(nind2, nind+1, 0)
    
    # calculate phenotype according to genotype
    if (sel.on) {
      pop2.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop2,
                  pop.geno = pop2.geno,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop2, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop2 <- pop2.pheno$pop
      pop2.pheno$pop <- NULL
      
      # reset trait
      if (mtd.reprod != "backcro") {
        trait <- list()
        trait$pop.sir1 <- pop1.pheno
        if (mtd.reprod == "tricro") {
          trait$pop.sir2 <- pop2.pheno
        } else {
          trait$pop.dam1 <- pop2.pheno
        }
      }
    }
    
    pop2.geno.em <- # genotype matrix after Mutation
          genotype(geno = pop2.geno,
                   incols = incols, 
                   blk.rg = blk.rg,
                   recom.spot = recom.spot,
                   range.hot = range.hot,
                   range.cold = range.cold,
                   # recom.cri = "cri3",
                   rate.mut = rate.mut, 
                   verbose = verbose)
    pop3.geno.em <- NULL
    pop4.geno.em <- NULL
  }

  if (mtd.reprod == "tricro" || mtd.reprod == "doubcro") {
    if (is.null(rawgeno3)) {
      logging.log(" --- base population 3 ---\n", verbose = verbose)
      prob1 <- runif(1)
      prob <- c(prob1, 1 - prob1)
      pop3.geno <- genotype(incols = incols, num.marker = num.marker, num.ind = num.ind, prob = prob, verbose = verbose)
    } else {
      pop3.geno <- genotype(rawgeno = rawgeno3, verbose = verbose)
    }

    # set base population information
    nind3 <- ncol(pop3.geno) / incols
    pop3 <- getpop(nind3, nind+nind2+1, 1)
    
    # calculate phenotype according to genotype
    if (sel.on) {
      pop3.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop3,
                  pop.geno = pop3.geno,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop3, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop3 <- pop3.pheno$pop
      pop3.pheno$pop <- NULL
      
      if (mtd.reprod == "tricro") {
        trait$pop.dam1 <- pop3.pheno
      } else {
        trait$pop.sir2 <- pop3.pheno
      }
    }
    
    pop3.geno.em <- # genotype matrix after Mutation
          genotype(geno = pop3.geno,
                   incols = incols, 
                   blk.rg = blk.rg,
                   recom.spot = recom.spot,
                   range.hot = range.hot,
                   range.cold = range.cold,
                   # recom.cri = "cri3",
                   rate.mut = rate.mut, 
                   verbose = verbose)
    pop4.geno.em <- NULL
  }

  if (mtd.reprod == "doubcro") {
    logging.log(" --- base population 4 ---\n", verbose = verbose)
    if (is.null(rawgeno4)) {
      prob1 <- runif(1)
      prob <- c(prob1, 1 - prob1)
      pop4.geno <- genotype(incols = incols, num.marker = num.marker, num.ind = num.ind, prob = prob, verbose = verbose)
    } else {
      pop4.geno <- genotype(rawgeno = rawgeno4, verbose = verbose)
    }

    # set base population information
    nind4 <- ncol(pop4.geno) / incols
    pop4 <- getpop(nind4, nind+nind2+nind3+1, 0)
    
    # calculate phenotype according to genotype
    if (sel.on) {
      pop4.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop4,
                  pop.geno = pop4.geno,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop4, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop4 <- pop4.pheno$pop
      pop4.pheno$pop <- NULL
      trait$pop.dam2 <- pop4.pheno
    }
    
    pop4.geno.em <- # genotype matrix after Mutation
          genotype(geno = pop4.geno,
                   incols = incols, 
                   blk.rg = blk.rg,
                   recom.spot = recom.spot,
                   range.hot = range.hot,
                   range.cold = range.cold,
                   # recom.cri = "cri3",
                   rate.mut = rate.mut, 
                   verbose = verbose)
  }

  ################### SETTING_PROCESS ###################
  # 1. setting of number of progenies in every generation.
  # 2. setting of directory

  # adjust for genetic correlation
  if (!(all(ps <= 1) | all(ps > 1))) stop("Please input a correct ps!")
  ps[1] <- ifelse(sel.on, ps[1], 1)
  ps[2] <- ifelse(sel.on, ps[2], 1)
  # calculate number of individuals in every generation
  count.ind <- rep(nind, num.gen)
  count.sir <- count.dam <- NULL
  if (mtd.reprod == "clone" || mtd.reprod == "dh" || mtd.reprod == "selfpol" || mtd.reprod == "randmate" || mtd.reprod == "randexself") {
    if (num.gen > 1) {
      count.sir <- ifelse(all(ps <= 1), round(nsir * ps[1]), ps[1])
      count.dam <- ifelse(all(ps <= 1), round(ndam * ps[2]), ps[2])
      count.ind[2] <- count.dam * num.prog
      if (num.gen > 2) {
        for(i in 3:num.gen) {
          count.sir[i-1] <- ifelse(all(ps <= 1), round(count.ind[i-1] * ratio * ps[1]), ps[1])
          count.dam[i-1] <- ifelse(all(ps <= 1), round(count.ind[i-1] * (1-ratio) * ps[2]), ps[2])
          count.ind[i] <- count.dam[i-1] * num.prog
        }
      }
    }
    
  } else if (mtd.reprod == "singcro") {
    count.sir <- ifelse(all(ps <= 1), round(nrow(basepop) * ps[1]), ps[1])
    count.dam <- ifelse(all(ps <= 1), round(nrow(pop2) * ps[2]), ps[2])
    sing.ind <- count.dam * num.prog
    count.ind <- c(nrow(basepop), nrow(pop2), sing.ind)

  } else if (mtd.reprod == "tricro") {
    num.sir2 <- ifelse(all(ps <= 1), round(nrow(pop2) * ps[1]), ps[1])
    num.dam1 <- ifelse(all(ps <= 1), round(nrow(pop3) * ps[2]), ps[2])
    dam21.ind <- num.dam1 * prog.tri
    num.sir1 <- ifelse(all(ps <= 1), round(nrow(basepop) * ps[1]), ps[1])
    num.dam21 <- ifelse(all(ps <= 1), round(dam21.ind * (1-ratio) * ps[2]), ps[2])
    tri.ind <- num.dam21 * num.prog
    count.sir <- c(num.sir2, num.sir1)
    count.dam <- c(num.dam1, num.dam21)
    count.ind <- c(nrow(basepop), nrow(pop2), nrow(pop3), dam21.ind, tri.ind)
    
  } else if (mtd.reprod == "doubcro") {
    num.sir1 <- ifelse(all(ps <= 1), round(nrow(basepop) * ps[1]), ps[1])
    num.dam1 <- ifelse(all(ps <= 1), round(nrow(pop2) * ps[2]), ps[2])
    sir11.ind <- num.dam1 * prog.doub
    num.sir2 <- ifelse(all(ps <= 1), round(nrow(pop3) * ps[1]), ps[1])
    num.dam2 <- ifelse(all(ps <= 1), round(nrow(pop4) * ps[2]), ps[2])
    dam22.ind <- num.dam2 * prog.doub
    num.sir11 <- ifelse(all(ps <= 1), round(sir11.ind * ratio * ps[2]), ps[2])
    num.dam22 <- ifelse(all(ps <= 1), round(dam22.ind * (1-ratio) * ps[2]), ps[2])
    doub.ind <- num.dam22 * num.prog
    count.sir <- c(num.sir1, num.sir2, num.sir11)
    count.dam <- c(num.dam1, num.dam2, num.dam22)
    count.ind <- c(nrow(basepop), nrow(pop2), nrow(pop3), nrow(pop4), sir11.ind, dam22.ind, doub.ind)
    
  } else if (mtd.reprod == "backcro") {
    count.ind[1] <- nrow(basepop) + nrow(pop2)
    if (num.gen > 1) {
      count.sir[1] <- ifelse(all(ps <= 1), round(nrow(basepop) * ps[1]), ps[1])
      count.dam[1] <- ifelse(all(ps <= 1), round(nrow(pop2) * ps[2]), ps[2])
      count.ind[2] <- count.dam[1] * num.prog
      for(i in 3:num.gen) {
        count.sir[i-1] <- count.sir[i-2]
        count.dam[i-1] <- ifelse(all(ps <= 1), round(count.ind[i-1] * (1-ratio) * ps[2]), ps[2])
        count.ind[i] <- count.dam[i-1] * num.prog
      }
    }
  } # end if mtd.reprod

  if (mtd.reprod != "userped") {
    # Create a folder to save files
    if (!is.null(outpath)) {
      if (!dir.exists(outpath)) stop("Please check your outpath!")
      if (out.format == "numeric") {
        outpath = paste0(outpath, .Platform$file.sep, sum(count.ind), "_Simer_Data_numeric")
      } else if (out.format == "plink"){
        outpath = paste0(outpath, .Platform$file.sep, sum(count.ind), "_Simer_Data_plink")
      } else {
        stop("out.format should be 'numeric' or 'plink'!")
      }
      if (!dir.exists(outpath)) dir.create(outpath)
      
      directory.rep <- paste0(outpath, .Platform$file.sep, "replication", replication)
      if (dir.exists(directory.rep)) {
        remove_bigmatrix(file.path(directory.rep, out))
        unlink(directory.rep, recursive = TRUE)
      }
      dir.create(directory.rep)
    }
  }

  if (all(ps <= 1)) {
    # calculate selection intensity
    sel.i <- dnorm(qnorm(1 -ps)) / ps 
    logging.log(" --- selection intensity ---\n", verbose = verbose)
    logging.log(" Selection intensity is", sel.i, "for males and females\n", verbose = verbose)
  } else if (all(ps > 1)) {
    sel.i <- ps
    logging.log(" --- selected individuals number ---\n", verbose = verbose)
    logging.log(" Number of selected individuals is", sel.i, "for males and females in every generation\n", verbose = verbose)
  }
 
  ################### REPRODUCTION_PROCESS ###################
  # 1. Reproduction based on basepop and basepop.geno according
  #    to different reproduction method.
  logging.log(" --- start reproduction ---\n", verbose = verbose)
  # multi-generation: clone, dh, selpol, randmate, randexself
	geno.back <- paste0(out, ".geno.bin")
	geno.desc <- paste0(out, ".geno.desc")
	ind.stays <- ind.stay <- NULL
	core.stays <- core.stay <- NULL
  if (mtd.reprod == "clone" || mtd.reprod == "dh" || mtd.reprod == "selfpol" || mtd.reprod == "randmate" || mtd.reprod == "randexself") {
    out.geno.gen <- out.geno.gen[out.geno.gen > 0]
    out.pheno.gen <- out.pheno.gen[out.pheno.gen > 0]
    out.geno.index <- getindex(count.ind, out.geno.gen)
    out.pheno.index <- getindex(count.ind, out.pheno.gen)

    # store all genotype
    geno.total.temp <- big.matrix(
      nrow = num.marker,
      ncol = outcols*sum(count.ind),
      init = 3,
      type = 'char')

    if (!is.null(outpath)) {
      geno.total <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * sum(count.ind[out.geno.gen]),
        init = 3,
        type = 'char',
        backingpath = directory.rep,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.total <- big.matrix(
        nrow = num.marker,
        ncol = outcols * sum(count.ind[out.geno.gen]),
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }

    # set total population
    pop.total <- basepop

    gc <- basepop.geno
    if (incols == 2 & outcols == 1) gc <- geno.cvt1(gc)
    if (1 %in% out.geno.gen) {
      input.geno(geno.total, gc, outcols*count.ind[1], mrk.dense)
    }
    input.geno(geno.total.temp, gc, outcols*count.ind[1], mrk.dense)

    logging.log(" After generation 1 ,", count.ind[1], "individuals are generated...\n", verbose = verbose)

    if (num.gen > 1) {
      # add selection to generation1
      if (sel.on) {
        ind.ordered <-
          selects(pop = basepop,
                  decr = decr,
                  sel.multi = sel.multi,
                  index.wt = index.wt,
                  index.tdm = index.tdm,
                  goal.perc = goal.perc,
                  pass.perc = pass.perc,
                  sel.sing = sel.sing,
                  pop.total = basepop,
                  pop.pheno = pop1.pheno, 
                  verbose = verbose)
        index.tdm <- ind.ordered[1]
        ind.ordered <- ind.ordered[-1]
      } else {
        ind.ordered <- basepop$index
      }
      core.stays[[1]] <- core.stay <- ind.stays[[1]] <- ind.stay <- getsd(ind.ordered, basepop, count.sir[1], count.dam[1])

      pop.last <- basepop
      pop.geno.last <- basepop.geno.em
      
      pop.geno.core <- basepop.geno.em[, getgmt(c(ind.stay$sir, ind.stay$dam), incols = incols)]
      pop1.geno.id <- basepop$index
      for (i in 2:num.gen) {
        pop.gp <- # pop.gp with genotype and pop information
          reproduces(pop1 = pop.last,
                     pop1.geno.id = pop1.geno.id, 
                     pop1.geno = pop.geno.last,
                     incols = incols, 
                     ind.stay = ind.stay,
                     mtd.reprod = mtd.reprod,
                     num.prog = num.prog,
                     ratio = ratio)
        
        pop.geno.curr <- pop.gp$geno
        pop.curr <- pop.gp$pop
        pop1.geno.id <- pop.curr$index
        isd <- c(2, 5, 6)
   
        # input genotype
        gc <- pop.geno.curr
        if (incols == 2 & outcols == 1) gc <- geno.cvt1(gc)
        if (i %in% out.geno.gen) {
          out.gg <- out.geno.gen[1:which(out.geno.gen == i)]
          input.geno(geno.total, gc, outcols * sum(count.ind[out.gg]), mrk.dense)
        }
        input.geno(geno.total.temp, gc, outcols*sum(count.ind[1:i]), mrk.dense)
        
        pop.total.temp <- rbind(pop.total[1:sum(count.ind[1:(i-1)]), isd], pop.curr[, isd])
        if (sel.on) {
          pop.pheno <-
            phenotype(effs = effs,
                      FR = FR,  
                      cv = cv, 
                      pop = pop.curr,
                      pop.geno = pop.geno.curr,
                      pos.map = pos.map,
                      var.pheno = var.pheno, 
                      h2.tr1 = h2.tr1,
                      gnt.cov = gnt.cov,
                      h2.trn = h2.trn, 
                      sel.crit = sel.crit, 
                      pop.total = pop.total.temp, 
                      sel.on = sel.on, 
                      inner.env =  inner.env, 
                      verbose = verbose)
          pop.curr <- pop.pheno$pop
          pop.pheno$pop <- NULL
        }
        
        pop.total <- rbind(pop.total, pop.curr)
       
        logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)

        if (i == num.gen) break
        
        # output index.tdm and ordered individuals indice
        if (sel.on) {
          ind.ordered <-
            selects(pop = pop.curr,
                    decr = decr,
                    sel.multi = sel.multi,
                    index.wt = index.wt,
                    index.tdm = index.tdm,
                    goal.perc = goal.perc,
                    pass.perc = pass.perc,
                    sel.sing = sel.sing,
                    pop.total = pop.total.temp,
                    pop.pheno = pop.pheno, 
                    verbose = verbose)
          index.tdm <- ind.ordered[1]
          ind.ordered <- ind.ordered[-1]
        } else {
          ind.ordered <- pop.curr$index
        }
        core.stays[[i]] <- core.stay <- ind.stays[[i]] <- ind.stay <- getsd(ind.ordered, pop.curr, count.sir[i], count.dam[i])

        pop.geno.last <-  # genotype matrix after Exchange and Mutation
          genotype(geno = pop.geno.curr,
                   incols = incols, 
                   blk.rg = blk.rg,
                   recom.spot = recom.spot,
                   range.hot = range.hot,
                   range.cold = range.cold,
                   # recom.cri = "cri3",
                   rate.mut = rate.mut, 
                   verbose = verbose)
        pop.last <- pop.curr
        
        if (sel.on & all(ps > 1)) {
          info.core <- sel.core(ind.stay, core.stay, refresh, keep.max.gen, 
              incols, pop.total = pop.total, pop.geno.curr, pop.geno.core)
          core.stays[[i]] <- core.stay <- info.core$core.stay
          pop.geno.last <- pop.geno.core <- info.core$core.geno
          pop1.geno.id <- c(core.stay$sir, core.stay$dam)       
          ind.stay <- core.stay
        }
      }  # end for
    }
    if(num.gen > 1) {
      names(ind.stays) <- paste0("gen", 1:(num.gen-1))
      names(core.stays) <- paste0("gen", 1:(num.gen-1))
    } 
    
    # if traits have genetic correlation
    # generate phenotype at last
    pop.pheno <-
      phenotype(effs = effs,
                FR = FR,  
                cv = cv, 
                pop = pop.total,
                pop.geno = geno.total.temp,
                pos.map = pos.map,
                var.pheno = var.pheno, 
                h2.tr1 = h2.tr1,
                gnt.cov = gnt.cov,
                h2.trn = h2.trn, 
                sel.crit = sel.crit, 
                pop.total = pop.total, 
                sel.on = FALSE, 
                inner.env =  inner.env, 
                verbose = verbose)
    pop.total <- pop.pheno$pop
    pop.pheno$pop <- NULL
    trait <- pop.pheno
    
    if (!is.null(outpath)) {
      # write files
      logging.log(" --- write files of total population ---\n", verbose = verbose)
      write.file(pop.total, geno.total, pos.map, out.geno.index, out.pheno.index, out, directory.rep, out.format, verbose)
      flush(geno.total)
    }
    
    if (num.gen > 1) {
      rm(pop.gp); rm(pop.curr); rm(pop.geno.curr); rm(pop.last); 
      rm(pop.geno.last); rm(pop.total.temp); 
    }
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(geno.total.temp); gc()
     
    # certain-generation: singcro, tricro, doubcro
  } else if (mtd.reprod == "singcro") {
    out.geno.index <- 1:sum(count.ind)
    logging.log(" After generation", 1, ",", sum(count.ind[1:2]), "individuals are generated...\n", verbose = verbose)
    
    if (!is.null(outpath)) {
      dir.sir <- paste0(directory.rep, .Platform$file.sep, count.ind[1], "_sir")
      dir.dam <- paste0(directory.rep, .Platform$file.sep, count.ind[2], "_dam")
      dir.sgc <- paste0(directory.rep, .Platform$file.sep, count.ind[3], "_single_cross")
      if (dir.exists(dir.sir)) { unlink(dir.sir, recursive = TRUE) }
      if (dir.exists(dir.dam)) { unlink(dir.dam, recursive = TRUE) }
      if (dir.exists(dir.sgc)) { unlink(dir.sgc, recursive = TRUE) }
      dir.create(dir.sir)
      dir.create(dir.dam)
      dir.create(dir.sgc)

      geno.sir <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char',
        backingpath = dir.sir,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char',
        backingpath = dir.dam,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.singcro <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char',
        backingpath = dir.sgc,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.sir <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char')
      geno.dam <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char')
      geno.singcro <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }
    
    if (sel.on) {
      # output index.tdm and ordered individuals indice
      ind.ordered <-
        selects(pop = basepop,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = basepop,
                pop.pheno = pop1.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered1 <- ind.ordered[-1]
      ind.ordered <-
        selects(pop = pop2,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop2,
                pop.pheno = pop2.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered2 <- ind.ordered[-1]
    } else {
      ind.ordered1 <- basepop$index
      ind.ordered2 <- pop2$index
    }
    ind.stays[[1]] <- getsd(ind.ordered1, basepop, count.sir, 0)
    ind.stays[[2]] <- getsd(ind.ordered2, pop2, 0, count.dam)
    names(ind.stays) <- c("basepop", "pop2")
    ind.stay$sir <- ind.stays[[1]]$sir
    ind.stay$dam <- ind.stays[[2]]$dam
    core.stays[[1]] <- ind.stay
    names(core.stays) <- "gen1"
      
    pop1.geno.id <- basepop$index
    pop2.geno.id <- pop2$index
    
    pop.gp <-
        reproduces(pop1 = basepop,
                   pop2 = pop2,
                   pop1.geno.id = basepop$index, 
                   pop2.geno.id = pop2$index, 
                   pop1.geno = basepop.geno.em,
                   pop2.geno = pop2.geno.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = mtd.reprod,
                   num.prog = num.prog,
                   ratio = ratio)

    pop.geno.singcro <- pop.gp$geno
    pop.singcro <- pop.gp$pop
    isd <- c(2, 5, 6)
    pop.total.temp <- rbind(basepop[, isd], pop2[, isd], pop.singcro[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.singcro,
                  pop.geno = pop.geno.singcro,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.singcro <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.singcro <- pop.pheno
    }
    
    logging.log(" After generation", 2, ",", sum(count.ind[1:3]), "individuals are generated...\n", verbose = verbose)

    gc.sir <- basepop.geno
    gc.dam <- pop2.geno
    gc.singcro <- pop.geno.singcro
    if (incols == 2 & outcols == 1) {
      gc.sir <- geno.cvt1(gc.sir)
      gc.dam <- geno.cvt1(gc.dam)
      gc.singcro <- geno.cvt1(gc.singcro)
    }
    input.geno(geno.sir, gc.sir, ncol(geno.sir), mrk.dense)
    input.geno(geno.dam, gc.dam, ncol(geno.dam), mrk.dense)
    input.geno(geno.singcro, gc.singcro, ncol(geno.singcro), mrk.dense)
    
    # if traits have genetic correlation
    # generate phenotype at last
    if (!sel.on) {
      pop.total <- rbind(basepop, pop2, pop.singcro)
      geno.total <- cbind(basepop.geno, pop2.geno, pop.geno.singcro)
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.total,
                  pop.geno = geno.total,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.total <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait <- pop.pheno
      basepop <- pop.total[1:nind, ]
      pop2 <- pop.total[(nind+1):(nind+nind2), ]
      pop.singcro <- pop.total[(nind+nind2+1):(nind+nind2+nrow(pop.singcro)), ]
    }
    
    if (!is.null(outpath)) {
      flush(geno.sir)
      flush(geno.dam)
      flush(geno.singcro)
      # write files
      logging.log(" --- write files of sirs ---\n", verbose = verbose)
      write.file(basepop, geno.sir, pos.map, 1:nrow(basepop), 1:nrow(basepop), out, dir.sir, out.format, verbose)
      logging.log(" --- write files of dams ---\n", verbose = verbose)
      write.file(pop2, geno.dam, pos.map, 1:nrow(pop2), 1:nrow(pop2), out, dir.dam, out.format, verbose)
      logging.log(" --- write files of progenies ---\n", verbose = verbose)
      write.file(pop.singcro, geno.singcro, pos.map, 1:nrow(pop.singcro), 1:nrow(pop.singcro), out, dir.sgc, out.format, verbose)
    }
    
    # set total information of population and genotype
    pop.total <- list(pop.sir1 = basepop, pop.dam1 = pop2, pop.singcro = pop.singcro)
    geno.total <- list(geno.sir1 = gc.sir, geno.dam1 = gc.dam, geno.singcro = gc.singcro)
    
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(pop2); rm(pop2.geno); rm(pop2.geno.em);
    rm(geno.sir); rm(geno.dam); rm(geno.singcro); rm(pop.gp); rm(pop.singcro); rm(pop.geno.singcro); 
    rm(gc.sir); rm(gc.dam); rm(gc.singcro); rm(pop.total.temp); gc()
    
  } else if (mtd.reprod == "tricro") {
    out.geno.index <- 1:sum(count.ind)
    logging.log(" After generation", 1, ",", sum(count.ind[1:3]), "individuals are generated...\n", verbose = verbose)
    
    if (!is.null(outpath)) {
      dir.sir1  <- paste0(directory.rep, .Platform$file.sep, count.ind[1], "_sir1")
      dir.dam1  <- paste0(directory.rep, .Platform$file.sep, count.ind[2], "_dam1")
      dir.sir2  <- paste0(directory.rep, .Platform$file.sep, count.ind[3], "_sir2")
      dir.dam21 <- paste0(directory.rep, .Platform$file.sep, count.ind[4], "_dam21")
      dir.trc   <- paste0(directory.rep, .Platform$file.sep, count.ind[5], "_three-ways_cross")
      if (dir.exists(dir.sir1))  { unlink(dir.sir1, recursive = TRUE) }
      if (dir.exists(dir.dam1))  { unlink(dir.dam1, recursive = TRUE) }
      if (dir.exists(dir.sir2))  { unlink(dir.sir2, recursive = TRUE) }
      if (dir.exists(dir.dam21)) { unlink(dir.dam21, recursive = TRUE) }
      if (dir.exists(dir.trc))   { unlink(dir.trc, recursive = TRUE) }
      dir.create(dir.sir1)
      dir.create(dir.dam1)
      dir.create(dir.sir2)
      dir.create(dir.dam21)
      dir.create(dir.trc)
    
      geno.sir1 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char',
        backingpath = dir.sir1,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam1 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char',
        backingpath = dir.dam1,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.sir2 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char',
        backingpath = dir.sir2,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam21 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[4],
        init = 3,
        type = 'char',
        backingpath = dir.dam21,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.tricro <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[5],
        init = 3,
        type = 'char',
        backingpath = dir.trc,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.sir1 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char')
      geno.dam1 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char')
      geno.sir2 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char')
      geno.dam21 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[4],
        init = 3,
        type = 'char')
      geno.tricro <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[5],
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }
    
    # correct the sex
    pop2$sex <- 1
    pop3$sex <- 2
    
    if (sel.on) {
      # add selection to generation1
      ind.ordered <-
        selects(pop = pop2,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop2,
                pop.pheno = pop2.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered1 <- ind.ordered[-1]
      ind.ordered <-
        selects(pop = pop3,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop3,
                pop.pheno = pop3.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered2 <- ind.ordered[-1]
    } else {
      ind.ordered1 <- pop2$index
      ind.ordered2 <- pop3$index
    }
    core.stays[[1]] <- core.stay <- ind.stays[[1]] <- getsd(ind.ordered1, pop2, count.sir[1], 0)
    core.stays[[2]] <- core.stay <- ind.stays[[2]] <- getsd(ind.ordered2, pop3, 0, count.dam[1])
    ind.stay$sir <- ind.stays[[1]]$sir
    ind.stay$dam <- ind.stays[[2]]$dam
    core.stays[[1]] <- ind.stay
    
    # the first generation to the second generation
    pop.gp <-
        reproduces(pop1 = pop2,
                   pop2 = pop3,
                   pop1.geno.id = pop2$index, 
                   pop2.geno.id = pop3$index, 
                   pop1.geno = pop2.geno.em,
                   pop2.geno = pop3.geno.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = "singcro",
                   num.prog = prog.tri,
                   ratio = ratio)

    pop.geno.dam21 <- pop.gp$geno
    pop.dam21 <- pop.gp$pop
    isd <- c(2, 5, 6)
    pop.total.temp <- rbind(basepop[, isd], pop2[, isd], pop3[, isd], pop.dam21[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.dam21,
                  pop.geno = pop.geno.dam21,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.dam21 <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.dam21 <- pop.pheno
      
      # output index.tdm and ordered individuals indice
      ind.ordered <-
        selects(pop = basepop,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = basepop,
                pop.pheno = pop1.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered1 <- ind.ordered[-1]
      ind.ordered <-
        selects(pop = pop.dam21,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop.total.temp,
                pop.pheno = pop.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered2 <- ind.ordered[-1]
    } else {
      ind.ordered1 <- basepop$index
      ind.ordered2 <- pop.dam21$index
    }
    ind.stays[[3]] <- getsd(ind.ordered1, basepop, count.sir[2], 0)
    ind.stays[[4]] <- getsd(ind.ordered2, pop.dam21, 0, count.dam[2])
    names(ind.stays) <- c("pop2", "pop3", "basepop", "pop.dam21")
    ind.stay$sir <- ind.stays[[3]]$sir
    ind.stay$dam <- ind.stays[[4]]$dam
    core.stays[[2]] <- ind.stay
    names(core.stays) <- c("gen1", "gen2")
    
    logging.log(" After generation", 2, ",", sum(count.ind[1:4]), "individuals are generated...\n", verbose = verbose)
    
    pop.geno.dam21.em <-  # genotype matrix after Exchange and Mutation
        genotype(geno = pop.geno.dam21,
                 incols = incols, 
                 blk.rg = blk.rg,
                 recom.spot = recom.spot,
                 range.hot = range.hot,
                 range.cold = range.cold,
                 # recom.cri = "cri3",
                 rate.mut = rate.mut, 
                 verbose = verbose)

    # the second generation to the third generation
    pop.gp <-
        reproduces(pop1 = basepop,
                   pop2 = pop.dam21,
                   pop1.geno.id = basepop$index, 
                   pop2.geno.id = pop.dam21$index, 
                   pop1.geno = basepop.geno.em,
                   pop2.geno = pop.geno.dam21.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = "singcro",
                   num.prog = num.prog,
                   ratio = ratio)

    pop.geno.tricro <- pop.gp$geno
    pop.tricro <- pop.gp$pop
    isd <- c(2, 5, 6)
    pop.total.temp <- rbind(pop.total.temp, pop.tricro[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.tricro,
                  pop.geno = pop.geno.tricro,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.tricro <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.tricro <- pop.pheno
    }
    
    logging.log(" After generation", 3, ",", sum(count.ind[1:5]), "individuals are generated...\n", verbose = verbose)

    gc.sir1 <- basepop.geno
    gc.sir2 <- pop2.geno
    gc.dam1 <- pop3.geno
    gc.dam21 <- pop.geno.dam21
    gc.tricro <- pop.geno.tricro
    if (incols == 2 & outcols == 1) {
      gc.sir1 <- geno.cvt1(gc.sir1)
      gc.sir2 <- geno.cvt1(gc.sir2)
      gc.dam1 <- geno.cvt1(gc.dam1)
      gc.dam21 <- geno.cvt1(gc.dam21)
      gc.tricro <- geno.cvt1(gc.tricro)
    }
    input.geno(geno.sir1, gc.sir1, ncol(geno.sir1), mrk.dense)
    input.geno(geno.sir2, gc.sir2, ncol(geno.dam1), mrk.dense)
    input.geno(geno.dam1, gc.dam1, ncol(geno.sir2), mrk.dense)
    input.geno(geno.dam21, gc.dam21, ncol(geno.dam21), mrk.dense)
    input.geno(geno.tricro, gc.tricro, ncol(geno.tricro), mrk.dense)
    
    # if traits have genetic correlation
    # generate phenotype at last
    if (!sel.on) {
      pop.total <- rbind(basepop, pop2, pop3, pop.dam21, pop.tricro)
      geno.total <- cbind(basepop.geno, pop2.geno, pop3.geno, pop.geno.dam21[], pop.geno.tricro[])
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.total,
                  pop.geno = geno.total,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.total <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait <- pop.pheno
      basepop <- pop.total[1:nind, ]
      pop2 <- pop.total[(nind+1):(nind+nind2), ]
      pop3 <- pop.total[(nind+nind2+1):(nind+nind2+nind3), ]
      pop.dam21 <- pop.total[(nind+nind2+nind3+1):(nind+nind2+nind3+nrow(pop.dam21)), ]
      pop.tricro <- pop.total[(nind+nind2+nind3+nrow(pop.dam21)+1):(nind+nind2+nind3+nrow(pop.dam21)+nrow(pop.tricro)), ]
    }
    
    if (!is.null(outpath)) {
      flush(geno.sir1)
      flush(geno.dam1)
      flush(geno.sir2)
      flush(geno.dam21)
      flush(geno.tricro)
      # write files
      logging.log(" --- write files of sir1s ---\n", verbose = verbose)
      write.file(basepop, geno.sir1, pos.map, 1:nrow(basepop), 1:nrow(basepop), out, dir.sir1, out.format, verbose)
      logging.log(" --- write files of sir2s ---\n", verbose = verbose)
      write.file(pop2, geno.sir2, pos.map, 1:nrow(pop2), 1:nrow(pop2), out, dir.sir2, out.format, verbose)
      logging.log(" --- write files of dam1s ---\n", verbose = verbose)
      write.file(pop3, geno.dam1, pos.map, 1:nrow(pop3), 1:nrow(pop3), out, dir.dam1, out.format, verbose)
      logging.log(" --- write files of dam21s ---\n", verbose = verbose)
      write.file(pop.dam21, geno.dam21, pos.map, 1:nrow(pop.dam21), 1:nrow(pop.dam21), out, dir.dam21, out.format, verbose)
      logging.log(" --- write files of progenies ---\n", verbose = verbose)
      write.file(pop.tricro, geno.tricro, pos.map, 1:nrow(pop.tricro), 1:nrow(pop.tricro), out, dir.trc, out.format, verbose)
    }
    
    # set total information of population and genotype
    pop.total <- list(pop.sir1 = basepop, pop.sir2 = pop2, pop.dam1 = pop3, pop.dam21 = pop.dam21, pop.tricro = pop.tricro)
    geno.total <- list(geno.sir1 = gc.sir1, geno.sir2 = gc.sir2, geno.dam1 = gc.dam1, geno.dam21 = gc.dam21, geno.tricro = gc.tricro)
    
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(pop2); rm(pop2.geno); rm(pop2.geno.em);
    rm(pop3); rm(pop3.geno); rm(pop3.geno.em); rm(geno.sir1); rm(geno.dam1); rm(geno.sir2);
    rm(pop.gp); rm(pop.dam21); rm(geno.dam21); rm(pop.geno.dam21); rm(pop.geno.dam21.em);
    rm(gc.sir1); rm(gc.sir2); rm(gc.dam1); rm(gc.dam21); rm(gc.tricro);
    rm(pop.tricro); rm(geno.tricro); rm(pop.total.temp); gc()

  } else if (mtd.reprod == "doubcro") {
    out.geno.index <- 1:sum(count.ind)
    logging.log(" After generation", 1, ",", sum(count.ind[1:4]), "individuals are generated...\n", verbose = verbose)
    
    if (!is.null(outpath)) {
      dir.sir1  <- paste0(directory.rep, .Platform$file.sep, count.ind[1], "_sir1")
      dir.dam1  <- paste0(directory.rep, .Platform$file.sep, count.ind[2], "_dam1")
      dir.sir2  <- paste0(directory.rep, .Platform$file.sep, count.ind[3], "_sir2")
      dir.dam2  <- paste0(directory.rep, .Platform$file.sep, count.ind[4], "_dam2")
      dir.sir11 <- paste0(directory.rep, .Platform$file.sep, count.ind[5], "_sir11")
      dir.dam22 <- paste0(directory.rep, .Platform$file.sep, count.ind[6], "_dam22")
      dir.dbc   <- paste0(directory.rep, .Platform$file.sep, count.ind[7], "_double_cross")
      if (dir.exists(dir.sir1))  { unlink(dir.sir1, recursive = TRUE) }
      if (dir.exists(dir.dam1))  { unlink(dir.dam1, recursive = TRUE) }
      if (dir.exists(dir.sir2))  { unlink(dir.sir2, recursive = TRUE) }
      if (dir.exists(dir.dam2))  { unlink(dir.dam2, recursive = TRUE) }
      if (dir.exists(dir.sir11)) { unlink(dir.sir11, recursive = TRUE) }
      if (dir.exists(dir.dam22)) { unlink(dir.dam22, recursive = TRUE) }
      if (dir.exists(dir.dbc))   { unlink(dir.dbc, recursive = TRUE) }
      dir.create(dir.sir1)
      dir.create(dir.dam1)
      dir.create(dir.sir2)
      dir.create(dir.dam2)
      dir.create(dir.sir11)
      dir.create(dir.dam22)
      dir.create(dir.dbc)
    
      geno.sir1 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char',
        backingpath = dir.sir1,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam1 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char',
        backingpath = dir.dam1,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.sir2 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char',
        backingpath = dir.sir2,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam2 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[4],
        init = 3,
        type = 'char',
        backingpath = dir.dam2,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.sir11 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[5],
        init = 3,
        type = 'char',
        backingpath = dir.sir11,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.dam22 <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[6],
        init = 3,
        type = 'char',
        backingpath = dir.dam22,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      geno.doubcro <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[7],
        init = 3,
        type = 'char',
        backingpath = dir.dbc,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.sir1 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[1],
        init = 3,
        type = 'char')
      geno.dam1 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[2],
        init = 3,
        type = 'char')
      geno.sir2 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[3],
        init = 3,
        type = 'char')
      geno.dam2 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[4],
        init = 3,
        type = 'char')
      geno.sir11 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[5],
        init = 3,
        type = 'char')
      geno.dam22 <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[6],
        init = 3,
        type = 'char')
      geno.doubcro <- big.matrix(
        nrow = num.marker,
        ncol = outcols * count.ind[7],
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }

    if (sel.on) {
      # add selection to generation1
      ind.ordered <-
        selects(pop = basepop,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = basepop,
                pop.pheno = pop1.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered1 <- ind.ordered[-1]
      ind.ordered <-
        selects(pop = pop2,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop2,
                pop.pheno = pop2.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered2 <- ind.ordered[-1]
    } else {
      ind.ordered1 <- basepop$index
      ind.ordered2 <- pop2$index
    }
    ind.stays[[1]]  <- getsd(ind.ordered1, basepop, count.sir[1], 0)
    ind.stays[[2]]  <- getsd(ind.ordered2, pop2, 0, count.dam[1])
    ind.stay$sir <- ind.stays[[1]]$sir
    ind.stay$dam <- ind.stays[[2]]$dam
    core.stays[[1]] <- ind.stay
    
    # the first generation to the second generation(the first two populations)
    pop.gp <-
        reproduces(pop1 = basepop,
                   pop2 = pop2,
                   pop1.geno.id = basepop$index, 
                   pop2.geno.id = pop2$index,
                   pop1.geno = basepop.geno.em,
                   pop2.geno = pop2.geno.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = "singcro",
                   num.prog = prog.doub,
                   ratio = ratio)

    pop.geno.sir11 <- pop.gp$geno
    pop.sir11 <- pop.gp$pop
    pop.sir11$index <- pop.sir11$index - pop.sir11$index[1] + 1 + pop4$index[length(pop4$index)]
    isd <- c(2, 5, 6)
    pop.total.temp <- rbind(basepop[, isd], pop2[, isd], pop3[, isd], pop4[, isd], pop.sir11[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.sir11,
                  pop.geno = pop.geno.sir11,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.sir11 <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.sir11 <- pop.pheno
      
      # output index.tdm and ordered individuals indice
      ind.ordered <-
        selects(pop = pop.sir11,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop.total.temp,
                pop.pheno = pop.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered.sir11 <- ind.ordered[-1]
    } else {
      ind.ordered.sir11 <- pop.sir11$index
    }
    ind.stays[[5]] <- getsd(ind.ordered.sir11, pop.sir11, count.sir[3], 0)

    pop.geno.sir11.em <-  # genotype matrix after Exchange and Mutation
        genotype(geno = pop.geno.sir11,
                 incols = incols, 
                 blk.rg = blk.rg,
                 recom.spot = recom.spot,
                 range.hot = range.hot,
                 range.cold = range.cold,
                 # recom.cri = "cri3",
                 rate.mut = rate.mut, 
                 verbose = verbose)

    if (sel.on) {
      # add selection to generation1
      ind.ordered <-
        selects(pop = pop3,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop3,
                pop.pheno = pop3.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered1 <- ind.ordered[-1]
      ind.ordered <-
        selects(pop = pop4,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop4,
                pop.pheno = pop4.pheno, 
                verbose = verbose)
      index.tdm <- ind.ordered[1]
      ind.ordered2 <- ind.ordered[-1]
    } else {
      ind.ordered1 <- pop3$index
      ind.ordered2 <- pop4$index
    }
    ind.stays[[3]] <- getsd(ind.ordered1, pop3, count.sir[2], 0)
    ind.stays[[4]] <- getsd(ind.ordered2, pop4, 0, count.dam[2])
    ind.stay$sir <- ind.stays[[3]]$sir
    ind.stay$dam <- ind.stays[[4]]$dam
    core.stays[[2]] <- ind.stay
    
    # the first generation to the second generation(the last two populations)
    pop.gp <-
        reproduces(pop1 = pop3,
                   pop2 = pop4,
                   pop1.geno.id = pop3$index, 
                   pop2.geno.id = pop4$index, 
                   pop1.geno = pop3.geno.em,
                   pop2.geno = pop4.geno.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = "singcro",
                   num.prog = prog.doub,
                   ratio = ratio)

    pop.geno.dam22 <- pop.gp$geno
    pop.dam22 <- pop.gp$pop
    pop.dam22$index <- pop.dam22$index - pop.dam22$index[1] + 1 + pop.sir11$index[length(pop.sir11$index)]
    pop.total.temp <- rbind(pop.total.temp, pop.dam22[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.dam22,
                  pop.geno = pop.geno.dam22,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.dam22 <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.dam22 <- pop.pheno
      
      # output index.tdm and ordered individuals indice
      ind.ordered <-
        selects(pop = pop.dam22,
                decr = decr,
                sel.multi = sel.multi,
                index.wt = index.wt,
                index.tdm = index.tdm,
                goal.perc = goal.perc,
                pass.perc = pass.perc,
                sel.sing = sel.sing,
                pop.total = pop.total.temp,
                pop.pheno = pop.pheno, 
                verbose = verbose)
      # index.tdm <- ind.ordered[1]
      ind.ordered.dam22 <- ind.ordered[-1]
    } else {
      ind.ordered.dam22 <- pop.dam22$index
    }
    ind.stays[[6]] <- getsd(ind.ordered.dam22, pop.dam22, 0, count.dam[3])
    names(ind.stays) <- c("basepop", "pop2", "pop3", "pop4", "pop.sir11", "pop.dam22")
    ind.stay$sir <- ind.stays[[5]]$sir
    ind.stay$dam <- ind.stays[[6]]$dam
    core.stays[[3]] <- ind.stay
    names(core.stays) <- c("gen1.0", "gen1.5", "gen2")
    
    pop.geno.dam22.em <-  # genotype matrix after Exchange and Mutation
        genotype(geno = pop.geno.dam22,
                 incols = incols, 
                 blk.rg = blk.rg,
                 recom.spot = recom.spot,
                 range.hot = range.hot,
                 range.cold = range.cold,
                 # recom.cri = "cri3",
                 rate.mut = rate.mut, 
                 verbose = verbose)

    logging.log(" After generation", 2, ",", sum(count.ind[1:6]), "individuals are generated...\n", verbose = verbose)
   
    # the second generation to the third generation
    pop.gp <-
        reproduces(pop1 = pop.sir11,
                   pop2 = pop.dam22,
                   pop1.geno.id = pop.sir11$index, 
                   pop2.geno.id = pop.dam22$index, 
                   pop1.geno = pop.geno.sir11.em,
                   pop2.geno = pop.geno.dam22.em,
                   incols = incols, 
                   ind.stay = ind.stay,
                   mtd.reprod = "singcro",
                   num.prog = num.prog,
                   ratio = ratio)

    pop.geno.doubcro <- pop.gp$geno
    pop.doubcro <- pop.gp$pop
    pop.total.temp <- rbind(pop.total.temp, pop.doubcro[, isd])
    
    if (sel.on) {
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.doubcro,
                  pop.geno = pop.geno.doubcro,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total.temp, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.doubcro <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait$pop.doubcro <- pop.pheno
    }
    
    logging.log(" After generation", 3, ",", sum(count.ind[1:7]), "individuals are generated...\n", verbose = verbose)

    gc.sir1 <- basepop.geno
    gc.dam1 <- pop2.geno
    gc.sir2 <- pop3.geno
    gc.dam2 <- pop4.geno
    gc.sir11 <- pop.geno.sir11
    gc.dam22 <- pop.geno.dam22
    gc.doubcro <- pop.geno.doubcro
    if (incols == 2 & outcols == 1) {
      gc.sir1 <- geno.cvt1(gc.sir1)
      gc.dam1 <- geno.cvt1(gc.dam1)
      gc.sir2 <- geno.cvt1(gc.sir2)
      gc.dam2 <- geno.cvt1(gc.dam2)
      gc.sir11 <- geno.cvt1(gc.sir11)
      gc.dam22 <- geno.cvt1(gc.dam22)
      gc.doubcro <- geno.cvt1(gc.doubcro)
    }
    input.geno(geno.sir1, gc.sir1, ncol(geno.sir1), mrk.dense)
    input.geno(geno.dam1, gc.dam1, ncol(geno.dam1), mrk.dense)
    input.geno(geno.sir2, gc.sir2, ncol(geno.sir2), mrk.dense)
    input.geno(geno.dam2, gc.dam2, ncol(geno.dam2), mrk.dense)
    input.geno(geno.sir11, gc.sir11, ncol(geno.sir11), mrk.dense)
    input.geno(geno.dam22, gc.dam22, ncol(geno.dam22), mrk.dense)
    input.geno(geno.doubcro, gc.doubcro, ncol(geno.doubcro), mrk.dense)
    
    # if traits have genetic correlation
    # generate phenotype at last
    if (!sel.on) {
      pop.total <- rbind(basepop, pop2, pop3, pop4, pop.sir11, pop.dam22, pop.doubcro)
      geno.total <- cbind(basepop.geno, pop2.geno, pop3.geno, pop4.geno, pop.geno.sir11[], pop.geno.dam22[], pop.geno.doubcro[])
      pop.pheno <-
        phenotype(effs = effs,
                  FR = FR,  
                  cv = cv, 
                  pop = pop.total,
                  pop.geno = geno.total,
                  pos.map = pos.map,
                  var.pheno = var.pheno, 
                  h2.tr1 = h2.tr1,
                  gnt.cov = gnt.cov,
                  h2.trn = h2.trn, 
                  sel.crit = sel.crit, 
                  pop.total = pop.total, 
                  sel.on = sel.on, 
                  inner.env =  inner.env, 
                  verbose = verbose)
      pop.total <- pop.pheno$pop
      pop.pheno$pop <- NULL
      trait <- pop.pheno
      basepop <- pop.total[1:nind, ]
      pop2 <- pop.total[(nind+1):(nind+nind2), ]
      pop3 <- pop.total[(nind+nind2+1):(nind+nind2+nind3), ]
      pop4 <- pop.total[(nind+nind2+nind3+1):(nind+nind2+nind3+nind4), ]
      pop.sir11 <- pop.total[(nind+nind2+nind3+nind4+1):(nind+nind2+nind3+nind4+nrow(pop.sir11)), ]
      pop.dam22 <- pop.total[(nind+nind2+nind3+nind4+nrow(pop.sir11)+1):(nind+nind2+nind3+nind4+nrow(pop.sir11)+nrow(pop.dam22)), ]
      pop.doubcro <- pop.total[(nind+nind2+nind3+nrow(pop.sir11)+nrow(pop.dam22)+1):(nind+nind2+nind3+nind4+nrow(pop.sir11)+nrow(pop.dam22)+nrow(pop.doubcro)), ]
    }
  
    if (!is.null(outpath)) {
      flush(geno.sir1)
      flush(geno.dam1)
      flush(geno.sir2)
      flush(geno.dam2)
      flush(geno.sir11)
      flush(geno.dam22)
      flush(geno.doubcro)
    
      # write files
      logging.log(" --- write files of sir1s ---\n", verbose = verbose)
      write.file(basepop, geno.sir1, pos.map, 1:nrow(basepop), 1:nrow(basepop), out, dir.sir1, out.format, verbose)
      logging.log(" --- write files of dam1s ---\n", verbose = verbose)
      write.file(pop2, geno.dam1, pos.map, 1:nrow(pop2), 1:nrow(pop2), out, dir.dam1, out.format, verbose)
      logging.log(" --- write files of sir2s ---\n", verbose = verbose)
      write.file(pop3, geno.sir2, pos.map, 1:nrow(pop3), 1:nrow(pop3), out, dir.sir2, out.format, verbose)
      logging.log(" --- write files of dam2s ---\n", verbose = verbose)
      write.file(pop4, geno.dam2, pos.map, 1:nrow(pop4), 1:nrow(pop4), out, dir.dam2, out.format, verbose)
      logging.log(" --- write files of sir11s ---\n", verbose = verbose)
      write.file(pop.sir11, geno.sir11, pos.map, 1:nrow(pop.sir11), 1:nrow(pop.sir11), out, dir.sir11, out.format, verbose)
      logging.log(" --- write files of dam22s ---\n", verbose = verbose)
      write.file(pop.dam22, geno.dam22, pos.map, 1:nrow(pop.dam22), 1:nrow(pop.dam22), out, dir.dam22, out.format, verbose)
      logging.log(" --- write files of progenies ---\n", verbose = verbose)
      write.file(pop.doubcro, geno.doubcro, pos.map, 1:nrow(pop.doubcro), 1:nrow(pop.doubcro), out, dir.dbc, out.format, verbose)
    }
  
    # set total information of population and genotype
    pop.total <- list(pop.sir1 = basepop, pop.dam1 = pop2, pop.sir2 = pop3, pop.dam2 = pop4, pop.sir11 = pop.sir11, pop.dam22 = pop.dam22, pop.doubcro = pop.doubcro)
    geno.total <- list(geno.sir1 = gc.sir1, geno.dam1 = gc.dam1, geno.sir2 = gc.sir2, geno.dam2 = gc.dam2, geno.sir11 = gc.sir11, geno.dam22 = gc.dam22, geno.doubcro = gc.doubcro)
    
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(geno.sir1);
    rm(pop2); rm(pop2.geno); rm(pop2.geno.em); rm(geno.dam1);
    rm(pop3); rm(pop3.geno); rm(pop3.geno.em); rm(geno.sir2);
    rm(pop4); rm(pop4.geno); rm(pop4.geno.em); rm(geno.dam2);
    rm(pop.sir11); rm(pop.geno.sir11); rm(pop.geno.sir11.em);
    rm(pop.dam22); rm(pop.geno.dam22); rm(pop.geno.dam22.em);
    rm(pop.gp); rm(pop.doubcro); rm(pop.geno.doubcro); 
    rm(gc.sir1); rm(gc.dam1); rm(gc.sir2); rm(gc.dam2);
    rm(gc.sir11); rm(gc.dam22); rm(gc.doubcro);
    rm(pop.total.temp); gc()

  } else if (mtd.reprod == "backcro") {
    if (num.gen != length(prog.back))
      stop(" Number of generation should equal to the length of prog.back!")
    out.geno.gen <- out.geno.gen[out.geno.gen > 0]
    out.pheno.gen <- out.pheno.gen[out.pheno.gen > 0]
    out.geno.index <- getindex(count.ind, out.geno.gen)
    out.pheno.index <- getindex(count.ind, out.pheno.gen)

    # store all genotype
    geno.total.temp <- big.matrix(
      nrow = num.marker,
      ncol = outcols * sum(count.ind),
      init = 3,
      type = 'char')
    
    if (!is.null(outpath)) {
      geno.total <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = outcols * sum(count.ind[out.geno.gen]),
        init = 3,
        type = 'char',
        backingpath = directory.rep,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.total <- big.matrix(
        nrow = num.marker,
        ncol = outcols * sum(count.ind[out.geno.gen]),
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }
 
    # set total population
    pop.total <- rbind(basepop, pop2)
    gc.base <- basepop.geno
    gc.pop2 <- pop2.geno
    if (incols == 2 & outcols == 1) {
      gc.base <- geno.cvt1(gc.base)
      gc.pop2 <- geno.cvt1(gc.pop2)
    }
    if (1 %in% out.geno.gen) {
      input.geno(geno.total, gc.base, outcols * nrow(basepop), mrk.dense)
      input.geno(geno.total, gc.pop2, outcols * count.ind[1], mrk.dense)
    }
    if (!sel.on) {
      input.geno(geno.total.temp, gc.base, outcols*nrow(basepop), mrk.dense)
      input.geno(geno.total.temp, gc.pop2, outcols*count.ind[1], mrk.dense)
    }
    logging.log(" After generation 1 ,", count.ind[1], "individuals are generated...\n", verbose = verbose)

    if (num.gen > 1) {
      if (sel.on) {
        # add selection to generation1
        ind.ordered <-
          selects(pop = basepop,
                  decr = decr,
                  sel.multi = sel.multi,
                  index.wt = index.wt,
                  index.tdm = index.tdm,
                  goal.perc = goal.perc,
                  pass.perc = pass.perc,
                  sel.sing = sel.sing,
                  pop.total = basepop,
                  pop.pheno = pop1.pheno, 
                  verbose = verbose)
        index.tdm <- ind.ordered[1]
        ind.ordered1 <- ind.ordered[-1] 
        ind.ordered <-
          selects(pop = pop2,
                  decr = decr,
                  sel.multi = sel.multi,
                  index.wt = index.wt,
                  index.tdm = index.tdm,
                  goal.perc = goal.perc,
                  pass.perc = pass.perc,
                  sel.sing = sel.sing,
                  pop.total = pop2,
                  pop.pheno = pop2.pheno, 
                  verbose = verbose)
        index.tdm <- ind.ordered[1]
        ind.ordered2 <- ind.ordered[-1]
      } else {
        ind.ordered1 <- basepop$index
        ind.ordered2 <- pop2$index
      }
      ind.stays[[1]] <- getsd(ind.ordered1, basepop, count.sir[1], 0)
      ind.stays[[2]] <- getsd(ind.ordered2, pop2, 0, count.dam[1])
      ind.stay$sir <- ind.stays[[1]]$sir
      ind.stay$dam <- ind.stays[[2]]$dam
      core.stays[[1]] <- ind.stay
      
      for (i in 2:num.gen) {
        pop.gp <-
          reproduces(pop1 = basepop,
                     pop2 = pop2,
                     pop1.geno.id = basepop$index, 
                     pop2.geno.id = pop2$index, 
                     pop1.geno = basepop.geno.em,
                     pop2.geno = pop2.geno.em,
                     incols = incols, 
                     ind.stay = ind.stay,
                     mtd.reprod = "singcro",
                     num.prog = num.prog,
                     ratio = ratio)
        
        pop.geno.curr <- pop.gp$geno
        pop.curr <- pop.gp$pop
        
        if (i %in% out.geno.gen) {
          gc <- pop.geno.curr
          if (incols == 2 & outcols == 1) gc <- geno.cvt1(gc)
          out.gg <- out.geno.gen[1:which(out.geno.gen == i)]
          input.geno(geno.total, gc, outcols * sum(count.ind[out.gg]), mrk.dense)
        }
        input.geno(geno.total.temp, pop.geno.curr, outcols*sum(count.ind[1:i]), mrk.dense)
        
        isd <- c(2, 5, 6)
        pop.total.temp <- rbind(pop.total[1:sum(count.ind[1:(i-1)]), isd], pop.curr[, isd])
        if (sel.on) {
          pop.pheno <-
            phenotype(effs = effs,
                      FR = FR,  
                      cv = cv, 
                      pop = pop.curr,
                      pop.geno = pop.geno.curr,
                      pos.map = pos.map,
                      var.pheno = var.pheno, 
                      h2.tr1 = h2.tr1,
                      gnt.cov = gnt.cov,
                      h2.trn = h2.trn, 
                      sel.crit = sel.crit, 
                      pop.total = pop.total.temp, 
                      sel.on = sel.on, 
                      inner.env =  inner.env, 
                      verbose = verbose)
          pop.curr <- pop.pheno$pop
          pop.pheno$pop <- NULL
        }
        
        pop.total <- rbind(pop.total, pop.curr)
        
        logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)
        
        if (i == num.gen) break
        
        if (sel.on) {
          # output index.tdm and ordered individuals indice
          ind.ordered <-
            selects(pop = pop.curr,
                    decr = decr,
                    sel.multi = sel.multi,
                    index.wt = index.wt,
                    index.tdm = index.tdm,
                    goal.perc = goal.perc,
                    pass.perc = pass.perc,
                    sel.sing = sel.sing,
                    pop.total = pop.total,
                    pop.pheno = pop.pheno, 
                    verbose = verbose)
          index.tdm <- ind.ordered[1]
          ind.ordered2 <- ind.ordered[-1]
        } else {
          ind.ordered2 <- pop.curr$index
        }
        ind.stays[[i+1]] <- getsd(ind.ordered2, pop.curr, 0, count.dam[i])
        ind.stay$sir <- ind.stays[[1]]$sir
        ind.stay$dam <- ind.stays[[i+1]]$dam
        core.stays[[i]] <- ind.stay
        
        pop2.geno.em <-  # genotype matrix after Exchange and Mutation
          genotype(geno = pop.geno.curr,
                   incols = incols, 
                   blk.rg = blk.rg,
                   recom.spot = recom.spot,
                   range.hot = range.hot,
                   range.cold = range.cold,
                   # recom.cri = "cri3",
                   rate.mut = rate.mut, 
                   verbose = verbose)
        pop2 <- pop.curr
      }  # end for
    }
    if(num.gen > 1) {
      names(ind.stays) <- c("basepop", "pop2", paste0("gen", 2:(num.gen-1)))
      names(core.stays) <- paste0("gen", 1:(num.gen-1))
    } 
    
    names(core.stays) <- paste0("gen", 1:(num.gen-1))
    
    # if traits have genetic correlation
    # generate phenotype at last
    pop.pheno <-
      phenotype(effs = effs,
                FR = FR,  
                cv = cv, 
                pop = pop.total,
                pop.geno = geno.total.temp,
                pos.map = pos.map,
                var.pheno = var.pheno, 
                h2.tr1 = h2.tr1,
                gnt.cov = gnt.cov,
                h2.trn = h2.trn, 
                sel.crit = sel.crit, 
                pop.total = pop.total, 
                sel.on = FALSE, 
                inner.env =  inner.env, 
                verbose = verbose)
    pop.total <- pop.pheno$pop
    pop.pheno$pop <- NULL
    trait <- pop.pheno
    
    if (!is.null(outpath)) {
      flush(geno.total)
      # write files
      logging.log(" --- write files of total population ---\n", verbose = verbose)
      write.file(pop.total, geno.total, pos.map, out.geno.index, out.pheno.index, out, directory.rep, out.format, verbose)
    }
    
    if (num.gen > 1) {
      rm(pop.gp); rm(pop.curr); rm(pop.geno.curr); rm(pop.total.temp)
    }
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(pop2); rm(pop2.geno); 
    rm(pop2.geno.em); rm(geno.total.temp); gc()

  } else if (mtd.reprod == "userped") {

    pop1.geno.copy <- basepop.geno

    if (is.null(userped)) {
      stop(" Please input pedigree in the process userped!")
    }
    rawped <- userped
    rawped[is.na(rawped)] <- "0"
    if (as.numeric(rawped[1, 2]) < basepop$index[1]) {
      stop(" The index of the first sir should be in index of pop1!")
    }

    # Thanks to YinLL for sharing codes of pedigree sorting
    pedx <- as.matrix(rawped)
    pedx0 <- c(setdiff(pedx[, 2],pedx[, 1]), setdiff(pedx[, 3],pedx[, 1]))

    if(length(pedx0) != 0){
      pedx <- rbind(cbind(pedx0, "0", "0"), pedx)
    }

    pedx <- pedx[pedx[, 1] != "0", ]
    pedx <- pedx[!duplicated(pedx), ]
    pedx <- pedx[!duplicated(pedx[, 1]), ]

    pedx1 <- cbind(1:(ncol(basepop.geno)/2), "0", "0")
    pedx2 <- pedx[!(pedx[, 2] == "0" & pedx[, 3] == "0"), ]
    go = TRUE
    i <- 1
    count.ind <- nrow(pedx1)
    logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)
    while(go == TRUE) {
      i <- i + 1
      Cpedx <- c(pedx1[, 1])
      idx <- (pedx2[, 2] %in% Cpedx) & (pedx2[, 3] %in% Cpedx)
      if (sum(idx) == 0) {
        logging.log(" Some individuals in pedigree are not in mating process!\n They are", verbose = verbose)
        simer.print(pedx2[, 1], verbose = verbose)
        pedx2 <- pedx2[-c(1:nrow(pedx2)), ]
      } else {
        index.sir <- as.numeric(pedx2[idx, 2])
        index.dam <- as.numeric(pedx2[idx, 3])
        pop.geno.curr <- mate(pop.geno = pop1.geno.copy, index.sir = index.sir, index.dam = index.dam)
        pop1.geno.copy <- cbind(pop1.geno.copy[], pop.geno.curr[])
        pedx1 <- rbind(pedx1, pedx2[idx, ])
        pedx2 <- pedx2[!idx, ]
        count.ind <- c(count.ind, length(index.dam))
        logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)
      }
      if (class(pedx2) == "character") pedx2 <- matrix(pedx2, 1)
      if (dim(pedx2)[1] == 0) go = FALSE
    }
    ped <- pedx1
    rm(pedx1);rm(pedx2);gc()

    # Create a folder to save files
    if (!is.null(outpath)) {
      if (!dir.exists(outpath)) stop("Please check your outpath!")
      if (out.format == "numeric") {
        outpath = paste0(outpath, .Platform$file.sep, sum(count.ind), "_Simer_Data_numeric")
      } else if (out.format == "plink"){
        outpath = paste0(outpath, .Platform$file.sep, sum(count.ind), "_Simer_Data_plink")
      } else {
        stop("out.format should be 'numeric' or 'plink'!")
      }
      if (!dir.exists(outpath)) { dir.create(outpath) }
      
      directory.rep <- paste0(outpath, .Platform$file.sep, "replication", replication)
      if (dir.exists(directory.rep)) {
        remove_bigmatrix(file.path(directory.rep, "simer"))
        unlink(directory.rep, recursive = TRUE)
      }
      dir.create(directory.rep)
    }

    index <- ped[, 1]
    out.geno.index <- index
    ped.sir <- ped[, 2]
    ped.dam <- ped[, 3]
    sex <- rep(0, length(index))
    sex[index %in% unique(ped.sir)] <- 1
    sex[index %in% unique(ped.dam)] <- 2
    sex[sex == 0] <- sample(1:2, sum(sex == 0), replace = TRUE)
    fam.temp <- getfam(ped.sir, ped.dam, 1, "pm")
    gen <- rep(1:length(count.ind), count.ind)
    pop.total <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    gc <- pop1.geno.copy
    if (incols == 2 & outcols == 1) gc <- geno.cvt1(gc)
    if (!is.null(outpath)) {
      geno.total <- filebacked.big.matrix(
        nrow = num.marker,
        ncol = ncol(gc),
        init = 3,
        type = 'char',
        backingpath = directory.rep,
        backingfile = geno.back,
        descriptorfile = geno.desc)
      options(bigmemory.typecast.warning=FALSE)
    } else {
      geno.total <- big.matrix(
        nrow = num.marker,
        ncol = ncol(gc),
        init = 3,
        type = 'char')
      options(bigmemory.typecast.warning=FALSE)
    }
    input.geno(geno.total, gc, ncol(geno.total), mrk.dense)

    isd <- c(2, 5, 6)
    pop.pheno <-
      phenotype(effs = effs,
                FR = FR,  
                cv = cv, 
                pop = pop.total,
                pop.geno = pop1.geno.copy,
                pos.map = pos.map,
                var.pheno = var.pheno, 
                h2.tr1 = h2.tr1,
                gnt.cov = gnt.cov,
                h2.trn = h2.trn, 
                sel.crit = sel.crit, 
                pop.total = pop.total[, isd], 
                sel.on = sel.on, 
                inner.env =  inner.env, 
                verbose = verbose)
    pop.total <- pop.pheno$pop
    pop.pheno$pop <- NULL
    trait <- pop.pheno
    
    if (!is.null(outpath)) {
      flush(geno.total)
      logging.log(" --- write files of total population ---\n", verbose = verbose)
      write.file(pop.total, geno.total, pos.map, index, index, out, directory.rep, out.format, verbose)
    }
    
    rm(basepop); rm(basepop.geno); rm(basepop.geno.em); rm(userped); rm(rawped); rm(ped); gc()

  } else {
    stop("Please input correct reproduction method!")
  }

  # total information list
  simer.list <- list(pop = pop.total, effs = effs, trait = trait, geno = geno.total, genoid = out.geno.index, map = pos.map, si = sel.i, ind.stays = ind.stays, core.stays = core.stays)
  rm(ind.stays); rm(effs); rm(trait); rm(pop.total); rm(geno.total); rm(input.map); rm(pos.map); gc()
  
  if (!is.null(selPath)) {
    goal.plan <- complan(simls = simer.list, FR = FR, index.wt = index.wt, decr = decr, selPath = selPath, verbose = verbose)
    simer.list$goal.plan <- goal.plan
    rm(goal.plan); gc()
  }

  print_accomplished(width = 70, verbose = verbose)
  # Return the last directory
  ed <- Sys.time()
  logging.log(" SIMER DONE WITHIN TOTAL RUN TIME:", format_time(as.numeric(ed)-as.numeric(op)), "\n", verbose = verbose)
  return(simer.list)
}
