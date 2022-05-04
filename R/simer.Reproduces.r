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


#' Reproduction
#' 
#' Population reproduction by different mate design.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 29, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "randmate")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run reproduction
#' SP <- reproduces(SP)
reproduces <- function(SP, ncpus = 0, verbose = TRUE) {

### Start reproduction

  reprod.way <- SP$reprod$reprod.way
  
  if (reprod.way == "clone") {
    SP <- mate.clone(SP, ncpus = ncpus, verbose = verbose)

  } else if (reprod.way == "dh") {
    SP <- mate.dh(SP, ncpus = ncpus, verbose = verbose)

  } else if (reprod.way == "selfpol") {
    SP <- mate.selfpol(SP, ncpus = ncpus, verbose = verbose)
  
  } else if (reprod.way == "randmate") {
    SP <- mate.randmate(SP, ncpus = ncpus, verbose = verbose)

  } else if (reprod.way == "randexself") {
    SP <- mate.randexself(SP, ncpus = ncpus, verbose = verbose)
  
  } else if (reprod.way == "2waycro") {
    SP <- mate.2waycro(SP, ncpus = ncpus, verbose = verbose)
    
  } else if (reprod.way == "3waycro") {
    SP <- mate.3waycro(SP, ncpus = ncpus, verbose = verbose)
    
  } else if (reprod.way == "4waycro") {
    SP <- mate.4waycro(SP, ncpus = ncpus, verbose = verbose)
  
  } else if (reprod.way == "backcro") {
    SP <- mate.backcro(SP, ncpus = ncpus, verbose = verbose)
  
  } else if (reprod.way == "userped") {
    SP <- mate.userped(SP, ncpus = ncpus, verbose = verbose)
    
  } else {
    stop("'reprod.way' should be 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro' or 'userped'!")
  }

  return(SP)
}

#' Mate
#' 
#' Mating according to the indice of sires and dams.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param pop.geno the genotype data.
#' @param index.sir the indice of sires.
#' @param index.dam the indice of dams.
#' @param incols '1':one-column genotype represents an individual; '2': two-column genotype represents an individual.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#'
#' @return a genotype matrix after mating
#' @export
#'
#' @examples
#' # Generate the genotype data
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)
#' SP <- genotype(SP)
#' pop.geno <- SP$geno$pop.geno$gen1
#' 
#' # The mating design
#' index.sir <- rep(1:50, each = 2)
#' index.dam <- rep(51:100, each = 2)
#' 
#' # Mate according to mating design
#' geno.curr <- mate(pop.geno = pop.geno, index.sir = index.sir,
#'                   index.dam = index.dam)
#' geno.curr[1:5, 1:5]
mate <- function(pop.geno, index.sir, index.dam, incols = 1, ncpus = 0) {
  
  pop.marker <- nrow(pop.geno)
  pop.geno.curr <- big.matrix(
      nrow = pop.marker,
      ncol = length(index.dam) * incols,
      init = 3,
      type = "char"
  )

  if (incols == 2) {
    s1 <- sample(c(0, 1), size = length(index.dam), replace=TRUE)
    s2 <- sample(c(0, 1), size = length(index.dam), replace=TRUE)
    gmt.sir <- index.sir * 2 - s1
    gmt.dam <- index.dam * 2 - s2
    gmt.comb <- c(gmt.sir, gmt.dam)
    gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
    gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
    
    BigMat2BigMat(pop.geno.curr@address, pop.geno@address, colIdx = gmt.comb, threads = ncpus)
    
  } else {
    # mix parent genotype to generate progeny genotype
    GenoMixer(pop.geno.curr@address, pop.geno@address, index.sir, index.dam, threads = ncpus)
  }
  
  return(pop.geno.curr)
}

#' Clone
#' 
#' Produce individuals by clone.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "clone")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run clone
#' SP <- mate.clone(SP)
mate.clone <- function(SP, ncpus = 0, verbose = TRUE) {
  
  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  for (i in 1:pop.gen) {
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- pop[, 1]
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
    incols <- SP$geno$incols
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
    if (is.null(pop.sel)) {
      ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
      ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      pop.sel <- list(sir = ind.sir , dam = ind.dam)
    }
    prog <- SP$reprod$prog
    
    pop.marker <- nrow(pop.geno)
    if (all(pop.sel$sir == pop.sel$dam)) {
      ped.dam <- pop.sel$dam
    } else {
      ped.dam <- c(pop.sel$sir, pop.sel$dam)
    }
    gmt.dam <- match(ped.dam, pop.geno.id)
    num.2ind <- length(ped.dam) * incols
    
    # pop.geno.curr <- matrix(3, nrow = pop.marker, ncol = num.2ind*prog)
    pop.geno.curr <- big.matrix(
      nrow = pop.marker,
      ncol = num.2ind*prog,
      init = 3,
      type = "char")
    
    if (incols == 2) {
      gmt.dam <- gmt.dam * 2
      gmt.sir <- gmt.dam - 1
      gmt.comb <- c(gmt.sir, gmt.dam)
      gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
      gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
    } else if (incols == 1) {
      gmt.comb <- gmt.dam
    } else {
      stop("Please input a correct incols!")
    }
    
    gmt.comb <- rep(gmt.comb, times = prog)
    BigMat2BigMat(pop.geno.curr@address, pop.geno@address, colIdx = gmt.comb, threads = ncpus)
    
    ped.sir <- rep(ped.dam, times = prog)
    ped.dam <- rep(ped.dam, times = prog)
    sex <- rep(0, length(ped.dam))
    index <- seq(pop$index[length(pop$index)]+1, length.out = length(ped.dam))
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, length(ped.dam))
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### phenotype and selection ###
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
  }
  
  return(SP)
}

#' Doubled haploid
#' 
#' Produce individuals by doubled haploid.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "dh")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run doubled haploid
#' SP <- mate.dh(SP)
mate.dh <- function(SP, ncpus = 0, verbose = TRUE) {

  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  for (i in 1:pop.gen) {
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- pop[, 1]
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
    incols <- SP$geno$incols
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
    if (is.null(pop.sel)) {
      ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
      ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      pop.sel <- list(sir = ind.sir , dam = ind.dam)
    }
    prog <- SP$reprod$prog
    
    pop.marker <- nrow(pop.geno)
    if (all(pop.sel$sir == pop.sel$dam)) {
      ped.dam <- pop.sel$dam
    } else {
      ped.dam <- c(pop.sel$sir, pop.sel$dam)
    }
    gmt.dam <- match(ped.dam, pop.geno.id)
    num.2ind <- length(ped.dam) * incols
    if (prog %% 2 != 0) {
      stop("prog should be an even in dh option!")
    }
    
    # pop.geno.curr <- matrix(3, nrow = pop.marker, ncol = num.2ind*prog)
    pop.geno.curr <- big.matrix(
      nrow = pop.marker,
      ncol = num.2ind*prog,
      init = 3,
      type = "char")
    
    gmt.dam <- match(ped.dam, pop.geno.id)
    if (incols == 2) {
      gmt.dam <- gmt.dam * 2
      gmt.sir <- gmt.dam - 1
      gmt.comb <- c(gmt.sir, gmt.dam)
      gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
      gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
    } else if (incols == 1) {
      gmt.comb <- gmt.dam
    } else {
      stop("Please input a correct incols!")
    }
    
    gmt.comb <- rep(rep(gmt.comb, each = 2), times = prog/2)
    BigMat2BigMat(pop.geno.curr@address, pop.geno@address, colIdx = gmt.comb, threads = ncpus)
    
    ped.sir <- rep(rep(ped.dam, each = 2), times = prog/2)
    ped.dam <- rep(rep(ped.dam, each = 2), times = prog/2)
    sex <- rep(0, length(ped.dam))
    index <- seq(pop$index[length(pop$index)]+1, length.out = length(ped.dam))
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, length(ped.dam))
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### phenotype and selection ###
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
  }
  
  return(SP)
}

#' Self-pollination
#' 
#' Produce individuals by self-pollination.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "selfpol")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run self-pollination
#' SP <- mate.selfpol(SP)
mate.selfpol <- function(SP, ncpus = 0, verbose = TRUE) {

  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  for (i in 1:pop.gen) {
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- pop[, 1]
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
    incols <- SP$geno$incols
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
    if (is.null(pop.sel)) {
      ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
      ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      pop.sel <- list(sir = ind.sir , dam = ind.dam)
    }
    prog <- SP$reprod$prog
    
    if (all(pop.sel$sir == pop.sel$dam)) {
      ped.dam <- pop.sel$dam
    } else {
      ped.dam <- c(pop.sel$sir, pop.sel$dam)
    }
    ped.sir <- rep(ped.dam, each = prog)
    ped.dam <- rep(ped.dam, each = prog)
    index.sir <- match(ped.dam, pop.geno.id)
    index.dam <- index.sir
    
    pop.geno.curr <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
    
    sex <- rep(0, length(index.dam))
    index <- seq(pop$index[length(pop$index)]+1, length.out = length(ped.dam))
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, length(ped.dam))
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### phenotype and selection ###
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
  }
  
  return(SP)
}

#' Random mating
#' 
#' Produce individuals by random-mating.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "randmate")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run random mating
#' SP <- mate.randmate(SP)
mate.randmate <- function(SP, ncpus = 0, verbose = TRUE) {
  
  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  for (i in 1:pop.gen) {
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- pop[, 1]
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
    incols <- SP$geno$incols
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
    if (is.null(pop.sel)) {
      ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
      ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      pop.sel <- list(sir = ind.sir , dam = ind.dam)
    }
    sex.rate <- SP$reprod$sex.rate
    prog <- SP$reprod$prog
    
    ped.sir <- pop.sel$sir
    ped.dam <- pop.sel$dam
    if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
    ped.sir <- sample(ped.sir, size = length(ped.dam), replace = TRUE)
    ped.sir <- rep(ped.sir, each = prog)
    ped.dam <- rep(ped.dam, each = prog)
    pop.ind <- length(ped.dam)
    
    index.sir <- match(ped.sir, pop.geno.id)
    index.dam <- match(ped.dam, pop.geno.id)
    
    pop.geno.curr <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
    
    if (all(pop$sex == 0)) {
      sex <- rep(0, length(ped.dam))
    } else {
      sex <- rep(2, pop.ind)
      sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
    }
    index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, pop.ind)
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### genotype, phenotype, and selection ###
    SP <- genotype(SP)
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
  }
  
  return(SP)
}

#' Random mating excluding self-pollination
#'
#' Produce individuals by random mating excluding self-pollination.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "randexself")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run random mating excluding self-pollination
#' SP <- mate.randexself(SP)
mate.randexself <- function(SP, ncpus = 0, verbose = TRUE) {
  
  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  for (i in 1:pop.gen) {
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- pop[, 1]
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
    incols <- SP$geno$incols
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
    if (is.null(pop.sel)) {
      ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
      ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      pop.sel <- list(sir = ind.sir , dam = ind.dam)
    }
    sex.rate <- SP$reprod$sex.rate
    prog <- SP$reprod$prog
    
    ped.sir <- pop.sel$sir
    ped.dam <- pop.sel$dam
    if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
    ped.sir <- sample(ped.sir, size = length(ped.dam), replace = TRUE)
    
    # To make sure ped.sir is different from ped.dam
    idt.flag <- ped.sir == ped.dam
    sum.idt <- sum(idt.flag)
    while (sum.idt != 0) {
      ped.sir[idt.flag] <- sample(ped.sir, size = sum.idt, replace = TRUE)
      idt.flag <- ped.sir == ped.dam
      sum.idt <- sum(idt.flag)
    }
    
    ped.sir <- rep(ped.sir, each = prog)
    ped.dam <- rep(ped.dam, each = prog)
    pop.ind <- length(ped.dam)
    
    index.sir <- match(ped.sir, pop.geno.id)
    index.dam <- match(ped.dam, pop.geno.id)
    
    pop.geno.curr <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
    
    if (all(pop$sex == 0)) {
      sex <- rep(0, length(ped.dam))
    } else {
      sex <- rep(2, pop.ind)
      sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
    }
    index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, pop.ind)
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### genotype, phenotype, and selection ###
    SP <- genotype(SP)
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
  }
  
  return(SP)
}

#' Two-way cross
#' 
#' Produce individuals by two-way cross.
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "2waycro")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run two-way cross
#' SP <- mate.2waycro(SP)
mate.2waycro <- function(SP, ncpus = 0, verbose = TRUE) {
  
  count.ind <- NULL
  
  # reproduction parameters
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  incols <- SP$geno$incols
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  if (is.null(pop.sel)) {
    ind.sir <- pop$index[pop$sex == 1]
    ind.dam <- pop$index[pop$sex == 2]
    pop.sel <- list(sir = ind.sir , dam = ind.dam)
  }
  sex.rate <- SP$reprod$sex.rate
  prog <- SP$reprod$prog
  
  count.ind <- c(count.ind, nrow(pop))
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  
  ped.sir <- pop.sel$sir
  ped.dam <- pop.sel$dam
  if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
  ped.sir <- sample(ped.sir, size = length(ped.dam), replace = TRUE)
  ped.sir <- rep(ped.sir, each = prog)
  ped.dam <- rep(ped.dam, each = prog)
  pop.ind <- length(ped.dam)
  
  index.sir <- match(ped.sir, pop.geno.id)
  index.dam <- match(ped.dam, pop.geno.id)
  
  pop.geno.curr <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
  
  sex <- rep(2, pop.ind)
  sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
  index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
  fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
  gen <- rep(pop$gen[1]+1, pop.ind)
  pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
  
  SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
  SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  
  ### genotype, phenotype, and selection ###
  SP <- genotype(SP)
  SP <- phenotype(SP)
  SP <- selects(SP)
  count.ind <- c(count.ind, nrow(pop.curr))
  logging.log(" After generation", 2, ",", sum(count.ind[1:2]), "individuals are generated...\n", verbose = verbose)
  
  return(SP)
}

#' Three-way cross
#'
#' Produce individuals by three-way cross.
#'
#' Build date: Apr 11, 2022
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "3waycro")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # three different breeds are cut by sex
#' SP$pheno$pop$gen1$sex <- rep(c(1, 2, 1), c(30, 30, 40))
#' # Run selection
#' SP <- selects(SP)
#' # Run three-way cross
#' SP <- mate.3waycro(SP)
mate.3waycro <- function(SP, ncpus = 0, verbose = TRUE) {
  
  count.ind <- NULL
  
  # reproduction parameters
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  incols <- SP$geno$incols
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  if (is.null(pop.sel)) {
    ind.sir <- pop$index[pop$sex == 1]
    ind.dam <- pop$index[pop$sex == 2]
    pop.sel <- list(sir = ind.sir , dam = ind.dam)
  }
  sex.rate <- SP$reprod$sex.rate
  prog <- SP$reprod$prog
  
  count.ind <- c(count.ind, nrow(pop))
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  
  sex1 <- pop$sex
  sex2 <- c(1, sex1[-length(sex1)])
  sex.op <- which(sex1 != sex2)
  if (length(sex.op) != 2) {
    stop("Something wrong in the format of three-way cross data!")
  }
  
  ### the first two-way cross ###
  ped.sir1 <- intersect(pop[1:(sex.op[1]-1), 1], pop.sel$sir)
  ped.dam1 <- intersect(pop[sex.op[1]:(sex.op[2]-1), 1], pop.sel$dam)
  ped.sir2 <- intersect(pop[sex.op[2]:nrow(pop), 1], pop.sel$sir)
  if (length(ped.sir1) == 1) ped.sir1 <- rep(ped.sir1, 2)
  ped.sir1 <- sample(ped.sir1, size = length(ped.dam1), replace = TRUE)
  ped.sir1 <- rep(ped.sir1, each = prog)
  ped.dam1 <- rep(ped.dam1, each = prog)
  pop.ind <- length(ped.dam1)
  
  index.sir1 <- match(ped.sir1, pop.geno.id)
  index.dam1 <- match(ped.dam1, pop.geno.id)
  
  pop.geno.dam2 <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir1, index.dam = index.dam1, ncpus = ncpus)
  
  sex <- rep(2, pop.ind)
  sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
  index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
  fam.temp <- getfam(ped.sir1, ped.dam1, pop$fam[length(pop$fam)]+1, "pm")
  gen <- rep(pop$gen[1]+1, pop.ind)
  pop.dam2 <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir1, dam = ped.dam1, sex = sex)
  
  SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.dam2
  SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.dam2
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  
  ### genotype, phenotype, and selection ###
  SP <- genotype(SP)
  SP <- phenotype(SP)
  SP <- selects(SP)
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- c(pop.geno.id, pop[, 1])
  pop.geno.curr <- big.matrix(
    nrow = nrow(pop.geno),
    ncol = length(pop.geno.id),
    init = 3,
    type = "char")
  BigMat2BigMat(pop.geno.curr@address, pop.geno@address, colIdx = 1:ncol(pop.geno), threads = ncpus)
  BigMat2BigMat(pop.geno.curr@address, pop.geno.dam2@address, colIdx = 1:ncol(pop.geno.dam2), op = ncol(pop.geno)+1, threads = ncpus)
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  
  count.ind <- c(count.ind, nrow(pop))
  logging.log(" After generation", 2, ",", sum(count.ind[1:2]), "individuals are generated...\n", verbose = verbose)
  
  ### the second two-way cross ###
  ped.dam2 <- pop.sel$dam
  if (length(ped.sir2) == 1) ped.sir2 <- rep(ped.sir2, 2)
  ped.sir2 <- sample(ped.sir2, size = length(ped.dam2), replace = TRUE)
  ped.sir2 <- rep(ped.sir2, each = prog)
  ped.dam2 <- rep(ped.dam2, each = prog)
  pop.ind <- length(ped.dam2)
  
  index.sir2 <- match(ped.sir2, pop.geno.id)
  index.dam2 <- match(ped.dam2, pop.geno.id)
  
  pop.geno.curr <- mate(pop.geno = pop.geno.curr, incols = incols, index.sir = index.sir2, index.dam = index.dam2, ncpus = ncpus)
  
  sex <- rep(2, pop.ind)
  sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
  index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
  fam.temp <- getfam(ped.sir2, ped.dam2, pop$fam[length(pop$fam)]+1, "pm")
  gen <- rep(pop$gen[1]+1, pop.ind)
  pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir2, dam = ped.dam2, sex = sex)
  
  SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
  SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  
  ### genotype, phenotype, and selection ###
  SP <- genotype(SP)
  SP <- phenotype(SP)
  SP <- selects(SP)
  count.ind <- c(count.ind, nrow(pop.curr))
  logging.log(" After generation", 3, ",", sum(count.ind[1:3]), "individuals are generated...\n", verbose = verbose)
  
  return(SP)
}

#' Four-way cross process
#'
#' Produce individuals by four-way cross.
#'
#' Build date: Apr 11, 2022
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "4waycro")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # four different breeds are cut by sex
#' SP$pheno$pop$gen1$sex <- rep(c(1, 2, 1, 2), c(25, 25, 25, 25))
#' # Run selection
#' SP <- selects(SP)
#' # Run four-way cross
#' SP <- mate.4waycro(SP)
mate.4waycro <- function(SP, ncpus = 0, verbose = TRUE) {
  
  count.ind <- NULL
  
  # reproduction parameters
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  incols <- SP$geno$incols
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  if (is.null(pop.sel)) {
    ind.sir <- pop$index[pop$sex == 1]
    ind.dam <- pop$index[pop$sex == 2]
    pop.sel <- list(sir = ind.sir , dam = ind.dam)
  }
  sex.rate <- SP$reprod$sex.rate
  prog <- SP$reprod$prog
  
  count.ind <- c(count.ind, nrow(pop))
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  
  sex1 <- pop$sex
  sex2 <- c(1, sex1[-length(sex1)])
  sex.op <- which(sex1 != sex2)
  if (length(sex.op) != 3) {
    stop("Something wrong in the format of four-way cross data!")
  }
  
  ### the first two two-way crosses ###
  ped.sir1 <- intersect(pop[1:(sex.op[1]-1), 1], pop.sel$sir)
  ped.dam1 <- intersect(pop[sex.op[1]:(sex.op[2]-1), 1], pop.sel$dam)
  ped.sir2 <- intersect(pop[sex.op[2]:(sex.op[3]-1), 1], pop.sel$sir)
  ped.dam2 <- intersect(pop[sex.op[3]:nrow(pop), 1], pop.sel$dam)
  if (length(ped.sir1) == 1) ped.sir1 <- rep(ped.sir1, 2)
  if (length(ped.sir2) == 1) ped.sir2 <- rep(ped.sir2, 2)
  ped.sir1 <- sample(ped.sir1, size = length(ped.dam1), replace = TRUE)
  ped.sir2 <- sample(ped.sir2, size = length(ped.dam2), replace = TRUE)
  ped.sir1 <- rep(ped.sir1, each = prog)
  ped.dam1 <- rep(ped.dam1, each = prog)
  ped.sir2 <- rep(ped.sir2, each = prog)
  ped.dam2 <- rep(ped.dam2, each = prog)
  pop.ind <- length(ped.dam1) + length(ped.dam2)
  
  index.sir1 <- match(ped.sir1, pop.geno.id)
  index.dam1 <- match(ped.dam1, pop.geno.id)
  index.sir2 <- match(ped.sir2, pop.geno.id)
  index.dam2 <- match(ped.dam2, pop.geno.id)
  
  pop.geno.sir11 <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir1, index.dam = index.dam1, ncpus = ncpus)
  pop.geno.dam22 <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir2, index.dam = index.dam2, ncpus = ncpus)
  pop.geno.curr <- big.matrix(
    nrow = nrow(pop.geno),
    ncol = ncol(pop.geno.sir11) + ncol(pop.geno.dam22),
    init = 3,
    type = "char")
  BigMat2BigMat(pop.geno.curr@address, pop.geno.sir11@address, colIdx = 1:ncol(pop.geno.sir11), threads = ncpus)
  BigMat2BigMat(pop.geno.curr@address, pop.geno.dam22@address, colIdx = 1:ncol(pop.geno.dam22), op = ncol(pop.geno.sir11)+1, threads = ncpus)
  
  sex <- rep(2, pop.ind)
  sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
  index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
  fam.temp <- getfam(c(ped.sir1, ped.sir2), c(ped.dam1, ped.dam2), pop$fam[length(pop$fam)]+1, "pm")
  gen <- rep(pop$gen[1]+1, pop.ind)
  pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = c(ped.sir1, ped.sir2), dam = c(ped.dam1, ped.dam2), sex = sex)
  
  SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
  SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  
  ### genotype, phenotype, and selection ###
  SP <- genotype(SP)
  SP <- phenotype(SP)
  SP <- selects(SP)
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  
  count.ind <- c(count.ind, nrow(pop))
  logging.log(" After generation", 2, ",", sum(count.ind[1:2]), "individuals are generated...\n", verbose = verbose)
  
  ### the third two-way cross ###
  ped.sir11 <- pop.sel$sir
  ped.dam22 <- pop.sel$dam
  if (length(ped.sir11) == 1) ped.sir11 <- rep(ped.sir11, 2)
  ped.sir11 <- sample(ped.sir11, size = length(ped.dam22), replace = TRUE)
  ped.sir11 <- rep(ped.sir11, each = prog)
  ped.dam22 <- rep(ped.dam22, each = prog)
  pop.ind <- length(ped.dam22)
  
  index.sir11 <- match(ped.sir11, pop.geno.id)
  index.dam22 <- match(ped.dam22, pop.geno.id)
  
  pop.geno.curr <- mate(pop.geno = pop.geno.curr, incols = incols, index.sir = index.sir11, index.dam = index.dam22, ncpus = ncpus)
  
  sex <- rep(2, pop.ind)
  sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
  index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
  fam.temp <- getfam(ped.sir11, ped.dam22, pop$fam[length(pop$fam)]+1, "pm")
  gen <- rep(pop$gen[1]+1, pop.ind)
  pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir11, dam = ped.dam22, sex = sex)
  
  SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
  SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
    names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
  
  ### genotype, phenotype, and selection ###
  SP <- genotype(SP)
  SP <- phenotype(SP)
  SP <- selects(SP)
  count.ind <- c(count.ind, nrow(pop.curr))
  logging.log(" After generation", 3, ",", sum(count.ind[1:3]), "individuals are generated...\n", verbose = verbose)
  
  return(SP)
}

#' Back cross
#' 
#' Produce individuals by back cross.
#'
#' Build date: Apr 12, 2022
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.gen}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "backcro")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run back cross
#' SP <- mate.backcro(SP)
mate.backcro <- function(SP, ncpus = 0, verbose = TRUE) {
  
  # reproduction parameters
  pop.gen <- SP$reprod$pop.gen - 1
  count.ind <- nrow(SP$pheno$pop[[length(SP$pheno$pop)]])
  logging.log(" After generation", 1, ",", sum(count.ind[1:1]), "individuals are generated...\n", verbose = verbose)
  if (pop.gen == 0) return(SP)
  
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  incols <- SP$geno$incols
  pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  if (is.null(pop.sel)) {
    ind.sir <- pop$index[pop$sex == 1]
    ind.dam <- pop$index[pop$sex == 2]
    pop.sel <- list(sir = ind.sir , dam = ind.dam)
  }
  sex.rate <- SP$reprod$sex.rate
  prog <- SP$reprod$prog
  
  # it is used in every generation
  ped.sir.ori <- pop.sel$sir
  pop.geno.id.ori <- pop.geno.id
  pop.geno.ori <- pop.geno
  
  for (i in 1:pop.gen) {
    ped.sir <- ped.sir.ori
    ped.dam <- pop.sel$dam
    if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
    ped.sir <- sample(ped.sir, size = length(ped.dam), replace = TRUE)
    ped.sir <- rep(ped.sir, each = prog)
    ped.dam <- rep(ped.dam, each = prog)
    pop.ind <- length(ped.dam)
    
    index.sir <- match(ped.sir, pop.geno.id)
    index.dam <- match(ped.dam, pop.geno.id)
    
    pop.geno.curr <- mate(pop.geno = pop.geno, incols = incols, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
    
    sex <- rep(2, pop.ind)
    sex[sample(1:pop.ind, pop.ind * sex.rate)] <- 1
    index <- seq(pop$index[length(pop$index)]+1, length.out = pop.ind)
    fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
    gen <- rep(pop$gen[1]+1, pop.ind)
    pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
    
    SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
    SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
    names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
      names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
    
    ### genotype, phenotype, and selection ###
    SP <- genotype(SP)
    SP <- phenotype(SP)
    SP <- selects(SP)
    count.ind <- c(count.ind, nrow(pop.curr))
    logging.log(" After generation", i + 1, ",", sum(count.ind[1:(i + 1)]), "individuals are generated...\n", verbose = verbose)
    
    pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
    pop.geno.id <- c(pop.geno.id.ori, pop[, 1])
    pop.geno <- big.matrix(
      nrow = nrow(pop.geno),
      ncol = ncol(pop.geno.ori) + ncol(pop.geno.curr),
      init = 3,
      type = "char")
    BigMat2BigMat(pop.geno@address, pop.geno.ori@address, colIdx = 1:ncol(pop.geno.ori), threads = ncpus)
    BigMat2BigMat(pop.geno@address, pop.geno.curr@address, colIdx = 1:ncol(pop.geno.curr), op = ncol(pop.geno.ori)+1, threads = ncpus)
    pop.sel <- SP$sel$pop.sel[[length(SP$sel$pop.sel)]]
  }
  
  return(SP)
}

#' User-specified pedigree mating
#'
#' Produce individuals by user-specified pedigree mating.
#'
#' Build date: Apr 12, 2022
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$reprod$pop.sel}{the generations of simulated population.}
#' \item{$reprod$reprod.way}{reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.}
#' \item{$reprod$sex.rate}{the sex ratio of simulated population.}
#' \item{$reprod$prog}{the progeny number of an individual.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = 10)
#' # Generate genotype simulation parameters
#' SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
#' # Generate phenotype simulation parameters
#' SP <- param.pheno(SP = SP, pop.ind = 100)
#' # Generate selection parameters
#' SP <- param.sel(SP = SP, sel.single = "comb")
#' # Generate reproduction parameters
#' SP <- param.reprod(SP = SP, reprod.way = "userped")
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
#' # Run genotype simulation
#' SP <- genotype(SP)
#' # Run phenotype simulation
#' SP <- phenotype(SP)
#' # Run selection
#' SP <- selects(SP)
#' # Run user-specified pedigree mating
#' SP <- mate.userped(SP)
mate.userped <- function(SP, ncpus = 0, verbose = TRUE) {
  
  # reproduction parameters
  pop <- SP$pheno$pop[[length(SP$pheno$pop)]]
  pop.geno.id <- pop[, 1]
  pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  incols <- SP$geno$incols
  userped <- SP$reprod$userped
  
  # thanks to YinLL for sharing codes of pedigree sorting
  userped[is.na(userped)] <- "0"
  pedx <- as.matrix(userped)
  pedx0 <- c(setdiff(pedx[, 2],pedx[, 1]), setdiff(pedx[, 3],pedx[, 1]))
  pedx0 <- pedx0[pedx0 != 0]
  
  if(length(pedx0) != 0){
    pedx <- rbind(cbind(pedx0, "0", "0"), pedx)
  }
  
  pedx <- pedx[!duplicated(pedx), ]
  pedx <- pedx[!duplicated(pedx[, 1]), ]
  
  pedx1 <- pedx[ (pedx[, 2] == "0" & pedx[, 3] == "0"), ]
  pedx2 <- pedx[!(pedx[, 2] == "0" & pedx[, 3] == "0"), ]
  go <- TRUE
  i <- 1
  count.ind <- nrow(pedx1)
  logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)
  while(go == TRUE) {
    i <- i + 1
    Cpedx <- c(pedx1[, 1])
    idx <- (pedx2[, 2] %in% Cpedx) & (pedx2[, 3] %in% Cpedx)
    if (sum(idx) == 0) {
      logging.log(" Some individuals in pedigree are not in mating process!\n They are", verbose = verbose)
      logging.print(pedx2[, 1], verbose = verbose)
      pedx2 <- pedx2[-c(1:nrow(pedx2)), ]
      
    } else {
      index <- pedx2[, 1]
      ped.sir <- pedx2[, 2]
      ped.dam <- pedx2[, 3]
      index.sir <- match(ped.sir, pop.geno.id)
      index.dam <- match(ped.dam, pop.geno.id)
      pop.ind <- length(index)
      
      pop.geno.curr <- mate(pop.geno = pop.geno, index.sir = index.sir, index.dam = index.dam, ncpus = ncpus)
      
      sex <- rep(0, length(index))
      sex[index %in% unique(pedx[, 2])] <- 1
      sex[index %in% unique(pedx[, 3])] <- 2
      sex[sex == 0] <- sample(1:2, sum(sex == 0), replace = TRUE)
      fam.temp <- getfam(ped.sir, ped.dam, pop$fam[length(pop$fam)]+1, "pm")
      gen <- rep(pop$gen[1]+1, pop.ind)
      pop.curr <- data.frame(index = index, gen = gen, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)
      
      SP$geno$pop.geno[[length(SP$geno$pop.geno) + 1]] <- pop.geno.curr
      SP$pheno$pop[[length(SP$pheno$pop) + 1]] <- pop.curr
      names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- 
        names(SP$pheno$pop)[length(SP$pheno$pop)] <- paste0("gen", length(SP$pheno$pop))
      
      ### genotype, phenotype, and selection ###
      SP <- genotype(SP)
      SP <- phenotype(SP)
      
      pedx1 <- rbind(pedx1, pedx2[idx, ])
      pedx2 <- pedx2[!idx, ]
      count.ind <- c(count.ind, pop.ind)
      logging.log(" After generation", i, ",", sum(count.ind[1:i]), "individuals are generated...\n", verbose = verbose)
    }
    if (any(class(pedx2) == "character")) pedx2 <- matrix(pedx2, 1)
    if (dim(pedx2)[1] == 0) go = FALSE
  }
  
  return(SP)
}

#' Family index and within-family index
#' 
#' Get indice of family and within-family
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param sir the indice of sires.
#' @param dam the indice of dams.
#' @param fam.op the initial index of family indice.
#' @param mode "pat": paternal mode; "mat": maternal mode; "pm": paternal and maternal mode.
#' 
#' @return a matrix with family indice and within-family indice.
#' @export
#'
#' @examples
#' s <- c(0, 0, 0, 0, 1, 3, 3, 1, 5, 7, 5, 7, 1, 3, 5, 7)
#' d <- c(0, 0, 0, 0, 2, 4, 4, 2, 6, 8, 8, 6, 6, 8, 4, 8)
#' fam <- getfam(sir = s, dam = d, fam.op = 1, mode = "pm")
#' fam
getfam <- function(sir, dam, fam.op, mode = c("pat", "mat", "pm")) {
  fam <- rep(0, length(sir))
  infam <- rep(0, length(sir))
  if (mode == "pat") {
    uni.sir <- unique(sir)
    for (i in 1:length(uni.sir)) {
      flag <- sir == uni.sir[i]
      fam[flag] <- fam.op + i - 1
      infam[flag] <- 1:sum(flag)
    }
  } else if (mode == "mat") {
    uni.dam <- unique(dam)
    for (i in 1:length(uni.dam)) {
      flag <- dam == uni.dam[i]
      fam[flag] <- fam.op + i - 1
      infam[flag] <- 1:sum(flag)
    }
  } else if (mode == "pm") {
    uni.sd <- unique(cbind(sir, dam))
    for (i in 1:nrow(uni.sd)) {
      flag1 <- sir == uni.sd[i, 1]
      flag2 <- dam == uni.sd[i, 2]
      flag <- flag1 & flag2
      fam[flag] <- fam.op + i - 1
      infam[flag] <- 1:sum(flag)
    }
  } else {
    stop("Please input right mode!")
  }

  return(cbind(fam, infam))
}

#' Individual number per generation
#' 
#' Calculate the individual number per generation.
#'
#' Build date: Apr 12, 2022
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param pop the population information containing environmental factors and other effects.
#' @param pop.gen the generations of simulated population.
#' @param ps if ps <= 1, fraction selected in selection of males and females; if ps > 1, ps is number of selected males and females.
#' @param reprod.way reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.
#' @param sex.rate the sex ratio of simulated population.
#' @param prog the progeny number of an individual.
#'
#' @return the vector containing the individual number per generation.
#' @export
#'
#' @examples
#' pop <- generate.pop(pop.ind = 100)
#' count.ind <- IndPerGen(pop)
IndPerGen <- function(pop, pop.gen = 2, ps = c(0.8, 0.8), reprod.way = "randmate", sex.rate = 0.5, prog = 2) {
  
  if (!(all(ps <= 1) | all(ps > 1))) {
    stop("Please input a correct ps!")
  }
  
  count.ind <- rep(nrow(pop), pop.gen)
  if (reprod.way == "clone" || reprod.way == "dh" || reprod.way == "selfpol" || reprod.way == "randmate" || reprod.way == "randexself") {
    if (pop.gen > 1) {
      ped.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      nDam <- ifelse(all(ps <= 1), round(length(ped.dam) * ps[2]), ps[2])
      if (reprod.way == "clone" || reprod.way == "dh" || reprod.way == "selfpol") {
        nDam <- ifelse(all(ps <= 1), round(nrow(pop) * ps[2]), ps[2])
        sex.rate <- 0
      }
      count.ind[2] <- nDam * prog
      if (pop.gen > 2) {
        for(i in 3:pop.gen) {
          nDam <- ifelse(all(ps <= 1), round(round(count.ind[i-1] * (1-sex.rate)) * ps[2]), ps[2])
          count.ind[i] <- nDam * prog
        }
      }
    }
    
  } else if (reprod.way == "2waycro") {
    ped.dam <- pop$index[pop$sex == 2]
    nDam <- ifelse(all(ps <= 1), round(length(ped.dam) * ps[2]), ps[2])
    count.ind <- c(nrow(pop), nDam * prog)
    
  } else if (reprod.way == "3waycro") {
    sex1 <- pop$sex
    sex2 <- c(1, sex1[-length(sex1)])
    sex.op <- which(sex1 != sex2)
    if (length(sex.op) != 2) {
      stop("Something wrong in the format of three-way cross data!")
    }
    ped.dam <- pop$index[pop$sex == 2]
    nDam <- ifelse(all(ps <= 1), round(length(ped.dam) * ps[2]), ps[2])
    count.ind <- c(nrow(pop), nDam * prog)
    nDam <- ifelse(all(ps <= 1), round(round(count.ind[2] * sex.rate) * ps[2]), ps[2])
    count.ind <- c(count.ind, nDam * prog)
    
  } else if (reprod.way == "4waycro") {
    sex1 <- pop$sex
    sex2 <- c(1, sex1[-length(sex1)])
    sex.op <- which(sex1 != sex2)
    if (length(sex.op) != 3) {
      stop("Something wrong in the format of four-way cross data!")
    }
    ped.dam1 <- pop[sex.op[1]:(sex.op[2]-1), 1]
    ped.dam2 <- pop[sex.op[3]:nrow(pop), 1]
    nDam1 <- ifelse(all(ps <= 1), round(length(ped.dam1) * ps[2]), ps[2])
    nDam2 <- ifelse(all(ps <= 1), round(length(ped.dam2) * ps[2]), ps[2])
    count.ind <- c(nrow(pop), (nDam1 + nDam2) * prog)
    nDam <- ifelse(all(ps <= 1), round(round(count.ind[2] * sex.rate) * ps[2]), ps[2])
    count.ind <- c(count.ind, nDam * prog)
    
  } else if (reprod.way == "backcro") {
    if (pop.gen > 1) {
      ped.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
      nDam <- ifelse(all(ps <= 1), round(length(ped.dam) * ps[2]), ps[2])
      count.ind[2] <- nDam * prog
      if (pop.gen > 2) {
        for(i in 3:pop.gen) {
          nDam <- ifelse(all(ps <= 1), round(round(count.ind[i-1] * (1-sex.rate)) * ps[2]), ps[2])
          count.ind[i] <- nDam * prog
        }
      }
    }
    
  } else if (reprod.way == "userped") {
    count.ind <- 0
    
  } else {
    stop("'reprod.way' should be 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro' or 'userped'!")
  }
  
  return(count.ind)
}
