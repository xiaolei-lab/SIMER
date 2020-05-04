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


#' Do reproducion by different mate design
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop2 population information of population2
#' @param pop1.geno.id genotype id of population1
#' @param pop2.geno.id genotype id of population2
#' @param pop1.geno genotype matrix of population1
#' @param pop2.geno genotype matrix of population2
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param mtd.reprod different reproduction methods with the options: "clone", "dh", "selfpol", "singcro", "randmate", and "randexself"
#' @param num.prog litter size of dams
#' @param ratio ratio of the males in all individuals
#'
#' @return population information and genotype matrix of current population
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' pop.list <- reproduces(pop1 = basepop,
#'                        pop1.geno.id = basepop$index,  
#'                        pop1.geno = basepop.geno,
#'                        ind.stay = list(sir=1:10, dam=11:90),
#'                        mtd.reprod = "randmate",
#'                        num.prog = 2,
#'                        ratio = 0.5)
#' pop <- pop.list$pop
#' geno <- pop.list$geno
#' str(pop)
#' str(geno)
reproduces <-
    function(pop1 = NULL,
             pop2 = NULL,
             pop1.geno.id = NULL,
             pop2.geno.id = NULL,
             pop1.geno = NULL,
             pop2.geno = NULL,
             incols = 2, 
             ind.stay = NULL,
             mtd.reprod = "randmate",
             num.prog = 2,
             ratio = 0.5) {

# Start reproduction

  if (mtd.reprod == "clone") {
    pop <- mate.clone(pop1, pop1.geno.id, pop1.geno, incols, ind.stay, num.prog)

  } else if (mtd.reprod == "dh") {
    pop <- mate.dh(pop1, pop1.geno.id, pop1.geno, incols, ind.stay, num.prog)

  } else if (mtd.reprod == "selfpol") {
    pop <- mate.selfpol(pop1, pop1.geno.id, pop1.geno, incols, ind.stay, num.prog, ratio)

  } else if (mtd.reprod == "singcro") {
    pop <- mate.singcro(pop1, pop2, pop1.geno.id, pop2.geno.id, pop1.geno, pop2.geno, incols, ind.stay, num.prog, ratio)

  } else if (mtd.reprod == "randmate") {
    pop <- mate.randmate(pop1, pop1.geno.id, pop1.geno, incols, ind.stay, num.prog, ratio)

  } else if (mtd.reprod == "randexself") {
    pop <- mate.randexself(pop1, pop1.geno.id, pop1.geno, incols, ind.stay, num.prog, ratio)

  } else {
    stop("Please input a right option within mtd.reprod!")
  }

  return(pop)
}

# mate process
#' Mating according to indice of sires and dams
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop.geno genotype matrix of population
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param index.sir indice of sires
#' @param index.dam indice of dams
#'
#' @return genotype matrix after mating
#' @export
#'
#' @examples
#' pop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' index.sir <- rep(c(1, 3, 4, 5, 7, 9), each = 2)
#' index.dam <- rep(c(11, 12, 13, 14, 16, 18), each = 2)
#' geno.curr <- mate(pop.geno = pop.geno, index.sir = index.sir,
#'                  index.dam = index.dam)
#' str(geno.curr)
mate <- function(pop.geno, incols = 2, index.sir, index.dam) {
  num.marker <- nrow(pop.geno)
  pop.geno.curr <- matrix(0, nrow = num.marker, ncol = length(index.dam) * incols)
  # pop.geno.curr <- big.matrix(
  #     nrow = num.marker,
  #     ncol = length(index.dam) * incols,
  #     type = "char")
  # options(bigmemory.typecast.warning=FALSE)

  if (incols == 2) {
    s1 <- sample(c(0, 1), size = length(index.dam), replace=TRUE)
    s2 <- sample(c(0, 1), size = length(index.dam), replace=TRUE)
    gmt.sir <- index.sir * 2 - s1
    gmt.dam <- index.dam * 2 - s2
    gmt.comb <- c(gmt.sir, gmt.dam)
    gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
    gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
    
    pop.geno.curr <- pop.geno[, gmt.comb]
    
  } else {
    # calculate weight of every marker
    num.block <- 100
    len.block <- num.marker %/% num.block
    tail.block <- num.marker %% num.block + len.block
    num.inblock <- c(rep(len.block, (num.block-1)), tail.block)
    accum.block <- Reduce("+", num.inblock, accumulate = TRUE)
    for (i in 1:100) {
      ed <- accum.block[i]
      op <- ed - num.inblock[i] + 1
      judpar <- sample(c(0, 1), length(index.dam), replace = TRUE)
      index.prog <- judpar * index.sir + (1-judpar) * index.dam
      pop.geno.curr[op:ed, ] <- pop.geno[op:ed, index.prog]
    }
  }
  
  return(pop.geno.curr)
}

#' Clone process
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop1.geno.id genotype id of population1
#' @param pop1.geno genotype matrix of population1
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#'
#' @return population information and genotype matrix of population after clone procecss
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir =  basepop$index, dam = basepop$index)
#' pop.clone <- mate.clone(pop1 = basepop, pop1.geno.id = basepop$index, 
#'     pop1.geno = basepop.geno, ind.stay = ind.stay, num.prog = 2)
#' pop <- pop.clone$pop
#' geno <- pop.clone$geno
#' str(pop)
#' str(geno)
mate.clone <- function(pop1, pop1.geno.id, pop1.geno, incols = 2, ind.stay, num.prog) {

  num.marker <- nrow(pop1.geno)
  ped.dam <- ind.stay$dam
  num.2ind <- length(ped.dam) * incols

  pop.geno.curr <- matrix(3, nrow = num.marker, ncol = num.2ind*num.prog)
  # pop.geno.curr <- big.matrix(
  #     nrow = num.marker,
  #     ncol = num.2ind*num.prog,
  #     type = "char")
  # options(bigmemory.typecast.warning=FALSE)
  
  if (incols == 2) {
    gmt.dam <- cal.genoloc(ped.dam, pop1.geno.id)
    gmt.dam <- gmt.dam * 2
    gmt.sir <- gmt.dam - 1
    gmt.comb <- c(gmt.sir, gmt.dam)
    gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
    gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
  } else {
    gmt.comb <- cal.genoloc(ped.dam, pop1.geno.id)
  }
  pop.geno.adj <- pop1.geno[, gmt.comb]
  
  for (i in 1:num.prog) {
    ed <- i * num.2ind
    op <- ed - num.2ind + 1
    pop.geno.curr[, op:ed] <- pop.geno.adj
    # input.geno(pop.geno.curr, pop.geno.adj, i * num.2ind, TRUE)
  }
  
  ped.sir <- rep(ped.dam, times = num.prog)
  ped.dam <- rep(ped.dam, times = num.prog)
  sex <- rep(0, length(ped.dam))
  index <- seq(pop1$index[length(pop1$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop1$fam[length(pop1$fam)]+1, "pm")
  gen <- rep(pop1$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.clone <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.clone)
}

#' Doubled haploid process
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop1.geno.id genotype id of population1
#' @param pop1.geno genotype matrix of population1
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#'
#' @return population information and genotype matrix of population after doubled haploid procecss
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir = basepop$index, dam = basepop$index)
#' pop.dh <- mate.dh(pop1 = basepop, pop1.geno.id = basepop$index, 
#'     pop1.geno = basepop.geno, ind.stay = ind.stay, num.prog = 2)
#' pop <- pop.dh$pop
#' geno <- pop.dh$geno
#' str(pop)
#' str(geno)
mate.dh <- function(pop1, pop1.geno.id, pop1.geno, incols = 2, ind.stay, num.prog) {

  num.marker <- nrow(pop1.geno)
  ped.dam <- ind.stay$dam
  num.2ind <- length(ped.dam) * incols
  if (num.prog %% 2 != 0) {
    stop("num.prog should be an even in dh option!")
  }

  pop.geno.curr <- matrix(3, nrow = num.marker, ncol = num.2ind*num.prog)
  # pop.geno.curr <- big.matrix(
  #    nrow = num.marker,
  #    ncol = num.2ind*num.prog,
  #    type = "char")
  # options(bigmemory.typecast.warning=FALSE)

  if (incols == 2) {
    gmt.dam <- cal.genoloc(ped.dam, pop1.geno.id)
    gmt.dam <- gmt.dam * 2
    gmt.sir <- gmt.dam - 1
    gmt.comb <- c(gmt.sir, gmt.dam)
    gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
    gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
  } else {
    gmt.comb <- cal.genoloc(ped.dam, pop1.geno.id)
  }
  pop.geno.adj <- pop1.geno[, gmt.comb]
  
  pop.geno.comb <- cbind(pop.geno.adj, pop.geno.adj)
  pop.geno.comb[, seq(1, ncol(pop.geno.comb), 2)] <- pop.geno.adj
  pop.geno.comb[, seq(2, ncol(pop.geno.comb), 2)] <- pop.geno.adj
  
  for (i in 2*1:(num.prog/2)) {
    ed <- i * num.2ind
    op <- ed - 2 * num.2ind + 1
    pop.geno.curr[, op:ed] <- pop.geno.comb
    # input.geno(pop.geno.curr, pop.geno.comb, i * num.2ind, TRUE)
  }
  
  ped.sir <- rep(rep(ped.dam, each = 2), times = num.prog/2)
  ped.dam <- rep(rep(ped.dam, each = 2), times = num.prog/2)
  sex <- rep(0, length(ped.dam))
  index <- seq(pop1$index[length(pop1$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop1$fam[length(pop1$fam)]+1, "pm")
  gen <- rep(pop1$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.dh <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.dh)
}

#' Self-pollination
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop1.geno.id genotype id of population1
#' @param pop1.geno genotype matrix of population1
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#' @param ratio ratio of males in all individuals
#'
#' @return population information and genotype matrix of population after self-pollination process
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir = basepop$index, dam = basepop$index)
#' pop.selfpol <- mate.selfpol(pop1 = basepop, pop1.geno.id = basepop$index, 
#'     pop1.geno = basepop.geno, ind.stay = ind.stay, num.prog = 2, ratio = 0.5)
#' pop <- pop.selfpol$pop
#' geno <- pop.selfpol$geno
#' str(pop)
#' str(geno)
mate.selfpol <- function(pop1, pop1.geno.id, pop1.geno, incols = 2, ind.stay, num.prog, ratio) {

  if (floor(num.prog * ratio) != num.prog * ratio) {
    stop("The product of num.prog and ratio should be a integer!")
  }

  ped.dam <- ind.stay$dam
  ped.sir <- rep(ped.dam, each = num.prog)
  ped.dam <- rep(ped.dam, each = num.prog)
  index.sir <- cal.genoloc(ped.dam, pop1.geno.id)
  index.dam <- index.sir

  pop.geno.curr <- mate(pop.geno = pop1.geno, incols = incols, index.sir = index.sir, index.dam = index.dam)

  sex <- rep(0, length(index.dam))
  index <- seq(pop1$index[length(pop1$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop1$fam[length(pop1$fam)]+1, "pm")
  gen <- rep(pop1$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.selfpol <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.selfpol)
}

#' Single cross process
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop2 population information of population2
#' @param pop1.geno.id genotype id of population1
#' @param pop2.geno.id genotype id of population2
#' @param pop1.geno genotype matrix of population1
#' @param pop2.geno genotype matrix of population2
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#' @param ratio ratio of males in all individuals
#'
#' @return population information and genotype matrix of population after single cross process
#' @export
#' 
#' @examples
#' pop1 <- getpop(nind = 100, from = 1, ratio = 0.1)
#' pop2 <- getpop(nind = 100, from = 101, ratio = 0.1)
#' pop1.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' pop2.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir = pop1$index, dam = pop2$index)
#' pop.singcro <- mate.singcro(pop1 = pop2, pop2 = pop2,
#'     pop1.geno.id = pop1$index, pop2.geno.id = pop2$index, 
#'     pop1.geno = pop1.geno, pop2.geno = pop2.geno,
#'     ind.stay = ind.stay, num.prog = 2, ratio = 0.5)
#' pop <- pop.singcro$pop
#' geno <- pop.singcro$geno
#' str(pop)
#' str(geno)
mate.singcro <- function(pop1, pop2, pop1.geno.id, pop2.geno.id, pop1.geno, pop2.geno, incols = 2, ind.stay, num.prog, ratio) {

  if (is.null(pop1) || is.null(pop2) || is.null(pop1.geno) || is.null(pop2.geno)) {
    stop("Only two breeds are needed in the single cross!")
  }
  if (nrow(pop1.geno) != nrow(pop2.geno)) {
    stop("Rows of genotype matrixs should be equal!")
  }
  if (floor(num.prog * ratio) != num.prog * ratio) {
    stop("The product of num.prog and ratio should be a integer!")
  }

  ped.sir <- ind.stay$sir
  ped.dam <- ind.stay$dam
  num.dam <- length(ped.dam)
  if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
  ped.sir <- sample(ped.sir, size = num.dam, replace = TRUE)
  ped.sir <- rep(ped.sir, each = num.prog)
  ped.dam <- rep(ped.dam, each = num.prog)

  index.sir <- cal.genoloc(ped.sir, pop1.geno.id)
  index.dam <- cal.genoloc(ped.dam, pop2.geno.id)
  pop12.geno <- cbind(pop1.geno[], pop2.geno[])
  
  pop.geno.curr <- mate(pop.geno = pop12.geno, incols = incols, index.sir = index.sir, index.dam = index.dam)

  sex <- rep(c(rep(1, num.prog*ratio), rep(2, num.prog*(1-ratio))), num.dam)
  index <- seq(pop2$index[length(pop2$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop2$fam[length(pop2$fam)]+1, "pm")
  gen <- rep(pop2$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.singcro <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.singcro)
}

#' Random mating process
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop1.geno.id genotype id of population1
#' @param pop1.geno genotype matrix of population1
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#' @param ratio ratio of males in all individuals
#'
#' @return population information and genotype matrix of population after random mating process
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir = basepop$index[basepop$sex == 1], 
#'                  dam = basepop$index[basepop$sex == 2])
#' pop.randmate <- mate.randexself(pop1 = basepop,
#'     pop1.geno.id = basepop$index, pop1.geno = basepop.geno, 
#'     ind.stay = ind.stay, num.prog = 2, ratio = 0.5)
#' pop <- pop.randmate$pop
#' geno <- pop.randmate$geno
#' str(pop)
#' str(geno)
mate.randmate <- function(pop1, pop1.geno.id, pop1.geno, incols = 2, ind.stay, num.prog, ratio) {

  if (floor(num.prog * ratio) != num.prog * ratio) {
    stop("The product of num.prog and ratio should be a integer!")
  }

  ped.sir <- ind.stay$sir
  ped.dam <- ind.stay$dam
  num.dam <- length(ped.dam)
  if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
  ped.sir <- sample(ped.sir, size = num.dam, replace = TRUE)
  ped.sir <- rep(ped.sir, each = num.prog)
  ped.dam <- rep(ped.dam, each = num.prog)
  index.sir <- cal.genoloc(ped.sir, pop1.geno.id)
  index.dam <- cal.genoloc(ped.dam, pop1.geno.id)

  pop.geno.curr <- mate(pop.geno = pop1.geno, incols = incols, index.sir = index.sir, index.dam = index.dam)

  if (all(pop1$sex == 0)) {
    sex <- rep(0, num.prog*num.dam)
  } else {
    sex <- rep(c(rep(1, num.prog*ratio), rep(2, num.prog*(1-ratio))), num.dam)
  }
  index <- seq(pop1$index[length(pop1$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop1$fam[length(pop1$fam)]+1, "pm")
  gen <- rep(pop1$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.randmate <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.randmate)
}

#' Random mating excluding self-pollination process
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param pop1 population information of population1
#' @param pop1.geno.id genotype id of population1
#' @param pop1.geno genotype matrix of population1
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param ind.stay selected individuals regarded as parents
#' @param num.prog litter size of dams
#' @param ratio ratio of males in all individuals
#'
#' @return population information and genotype matrix of population after random mating excluding self-pollination process
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- genotype(num.marker = 48353, num.ind = 100, verbose = TRUE)
#' ind.stay <- list(sir = basepop$index[basepop$sex == 1], 
#'                  dam = basepop$index[basepop$sex == 2])
#' pop.randexself <- mate.randexself(pop1 = basepop,
#'     pop1.geno.id = basepop$index, pop1.geno = basepop.geno, 
#'     ind.stay = ind.stay, num.prog = 2, ratio = 0.5)
#' pop <- pop.randexself$pop
#' geno <- pop.randexself$geno
#' str(pop)
#' str(geno)
mate.randexself <- function(pop1, pop1.geno.id, pop1.geno, incols = 2, ind.stay, num.prog, ratio) {

  if (floor(num.prog * ratio) != num.prog * ratio) {
    stop("The product of num.prog and ratio should be a integer!")
  }

  ped.sir <- ind.stay$sir
  ped.dam <- ind.stay$dam
  num.dam <- length(ped.dam)
  if (length(ped.sir) == 1) ped.sir <- rep(ped.sir, 2)
  ped.sir <- sample(ped.sir, size = num.dam, replace = TRUE)

  # To make sure ped.sir is different from ped.dam
	idt.flag <- ped.sir == ped.dam
  sum.idt <- sum(idt.flag)
	while (sum.idt != 0) {
  	ped.sir[idt.flag] <- sample(ped.sir, size = sum.idt, replace = TRUE)
  	idt.flag <- ped.sir == ped.dam
    sum.idt <- sum(idt.flag)
	}

  ped.sir <- rep(ped.sir, each = num.prog)
  ped.dam <- rep(ped.dam, each = num.prog)
  index.sir <- cal.genoloc(ped.sir, pop1.geno.id)
  index.dam <- cal.genoloc(ped.dam, pop1.geno.id)

  pop.geno.curr <- mate(pop.geno = pop1.geno, incols = incols, index.sir = index.sir, index.dam = index.dam)

  if (all(pop1$sex == 0)) {
    sex <- rep(0, num.prog*num.dam)
  } else {
    sex <- rep(c(rep(1, num.prog*ratio), rep(2, num.prog*(1-ratio))), num.dam)
  }
  index <- seq(pop1$index[length(pop1$index)]+1, length.out = length(ped.dam))
  fam.temp <- getfam(ped.sir, ped.dam, pop1$fam[length(pop1$fam)]+1, "pm")
  gen <- rep(pop1$gen[1]+1, length(ped.dam))
  pop.curr <- data.frame(gen = gen, index = index, fam = fam.temp[, 1], infam = fam.temp[, 2], sir = ped.sir, dam = ped.dam, sex = sex)

  list.randexself <- list(geno = pop.geno.curr, pop = pop.curr)
  return(list.randexself)
}

#' Get indice of family and within-family
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param sir indice of sires
#' @param dam indice of dams
#' @param fam.op initial index of family indice
#' @param mode mode to get indice with "pat", "mat" and "pm"
#'
#' @return a matrix with family indice and within-family indice
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

#' Generate population according to number of individauls
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param nind number of the individuals in a population
#' @param from initial index of the population
#' @param ratio ratio of males in a population
#' @param gen generation ID of the population
#'
#' @return population information
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' str(basepop)
getpop <- function(nind, from, ratio, gen = 1) {
  pop <- data.frame(
    gen   = rep(gen, nind),
    index = seq(from = from, length.out = nind),
    fam   = seq(from = from, length.out = nind),
    infam = seq(from = from, length.out = nind),
    sir   = rep(0, nind),
    dam   = rep(0, nind),
    sex   = rep(1:2, c(floor(ratio*nind), (nind - floor(ratio*nind))))
  )
  return(pop)
}

#' Get output index of phenotype and pedigree
#'
#' Build date: Nov 14, 2018
#' Last update: Aug 1, 2019
#'
#' @author Dong Yin
#'
#' @param count.ind a vector with population size of every generations
#' @param out.gen a vector with generations needing ouputting
#'
#' @return individuals indice need outputting
#' @export
#' 
#' @examples
#' count.ind <- c(100, 200, 400, 800)
#' out.gen <- 2:4
#' idx <- getindex(count.ind = count.ind, out.gen = out.gen)
#' str(idx)
getindex <- function(count.ind, out.gen) {
  out.index <- NULL
  for (i in 1:length(out.gen)) {
      if (out.gen[i] == 1) {
        op <- 1
      } else {
        op <- sum(count.ind[1:(out.gen[i]-1)]) + 1
      }
      ed <- op + count.ind[out.gen[i]] - 1
      out.index <- c(out.index, op:ed)
  }
  return(out.index)
}

#' Get selected sires and dams
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 30, 2020
#'
#' @author Dong Yin
#'
#' @param ind.ordered individual ID after sorting
#' @param pop population information
#' @param num.sir number of sires
#' @param num.dam number of dams
#'
#' @return sire ID and dam ID after selection
#' @export
#'
#' @examples
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.5)
#' ind.ordered <- sample.int(100)
#' 
#' ind.stay <- getsd(ind.ordered = ind.ordered, pop = basepop, 
#'                   num.sir = 10, num.dam = 40)
#' str(ind.stay)
getsd <- function(ind.ordered, pop, num.sir, num.dam) {
  ind.sir <- pop$index[pop$sex == 1 | pop$sex == 0]
  ind.dam <- pop$index[pop$sex == 2 | pop$sex == 0]
  sir.ordered <- intersect(ind.ordered, ind.sir)
  dam.ordered <- intersect(ind.ordered, ind.dam)
  if (length(sir.ordered) == 0 | num.sir == 0) {
    sir.stay <- NULL
  } else {
    if (num.sir <= length(ind.sir)) {
      sir.stay <- sir.ordered[1:num.sir]
    } else {
      sir.stay <- c(sir.ordered, sample(sir.ordered, num.sir-length(ind.sir), replace = TRUE))
    }
  }
  if (length(dam.ordered) == 0 | num.dam == 0) {
    dam.stay <- NULL
  } else {
    if (num.dam <= length(ind.dam)) {
      dam.stay <- dam.ordered[1:num.dam]
    } else {
      dam.stay <- c(dam.ordered, sample(dam.ordered, num.dam-length(ind.dam), replace = TRUE))
    }
  }
  ind.stay <- list(sir = sir.stay, dam = dam.stay)
  return(ind.stay)
}
