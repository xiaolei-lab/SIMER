
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


#' Generate genotype and editing genotype
#'
#' Build date: Nov 14, 2018
#' Last update: Oct 13, 2019
#'
#' @author Dong Yin
#'
#' @param rawgeno extrinsic genotype matrix
#' @param geno genotype matrix need dealing with
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#' @param num.marker number of the markers
#' @param num.ind population size of the base population
#' @param prob weight of "0" and "1" in genotype matrix, the sum of elements in vector equal to 1
#' @param blk.rg it represents the starting position and the ending position of a block
#' @param recom.spot whether to consider recombination in every block
#' @param range.hot range of number of chromosome crossovers in a hot spot block
#' @param range.cold range of number of chromosome crossovers in a cold spot block
#' @param rate.mut mutation rate between 1e-8 and 1e-6
#' @param verbose whether to print detail
#'
#' @return a genotype matrix with block and map information
#' @export
#'
#' @examples
#' # get map file and create genotype matrix
#' data(simdata)
#' num.marker = nrow(input.map)
#' num.ind = 100
#' incols = 2
#' basepop.geno <- genotype(num.marker = num.marker, num.ind = num.ind, verbose = TRUE)
#' basepop.geno[1:5, 1:5]
#'
#' # get block information and recombination information
#' nmrk <- nrow(basepop.geno)
#' nind <- ncol(basepop.geno) / incols
#' pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
#' blk.rg <- cal.blk(pos.map)
#' recom.spot <- as.numeric(pos.map[blk.rg[, 1], 7])
#'
#' # genotype matrix after Exchange and Mutation
#' basepop.geno.em <-
#' genotype(geno = basepop.geno,
#'          incols = 2, 
#'          blk.rg = blk.rg,
#'          recom.spot = recom.spot,
#'          range.hot = 4:6,
#'          range.cold = 1:5,
#'          rate.mut = 1e-8, 
#'          verbose = TRUE)
#' basepop.geno.em[1:5, 1:5]
genotype <-
    function(rawgeno = NULL,
             geno = NULL,
             incols = 2, 
             num.marker = NULL,
             num.ind = NULL,
             prob = c(0.5, 0.5),
	           blk.rg = NULL,
             recom.spot = NULL,
             range.hot = 4:6,
             range.cold = 1:5,
             rate.mut = 1e-8, 
             verbose = TRUE) {

# Start genotype

  if (is.null(rawgeno) & is.null(geno) & is.null(num.marker) & is.null(num.ind)) {
    stop("Please input information of genotype!")
  }

  if (!is.null(rawgeno)) {
    logging.log(" Input outer genotype matrix...\n", verbose = verbose)
    if (!is.matrix(rawgeno)) {
      rawgeno <- as.matrix(rawgeno)
    }
    outgeno <- rawgeno

  } else if (!is.null(num.marker) && !is.null(num.ind)){
    logging.log(" Establish genotype matrix of base-population...\n", verbose = verbose)
    if (incols == 2) {
      codes <- c(0, 1)
    } else if (incols == 1) {
      codes <- c(0, 1, 2)
      if (!is.null(prob))
        prob <- c(prob[1]^2, 2*prob[1]*prob[2], prob[2]^2)
    } else {
      stop("incols should only be 1 or 2!")
    }
    outgeno <- matrix(sample(codes, num.marker*num.ind*incols, prob = prob, replace = TRUE), num.marker, incols*num.ind)
    
  } else if (!is.null(geno) & !is.null(recom.spot) & incols == 2) {
    logging.log(" Chromosome exchange and mutation on genotype matrix...\n", verbose = verbose)
    num.marker <- nrow(geno)
    num.ind <- ncol(geno) / incols
    # outgeno <- deepcopy(geno) # deepcopy() in bigmemory
    # ind.swap <- sample(c(0, 1), num.ind, replace = TRUE)

    geno.swap <- function(ind) {
      # if (ind.swap[ind] == 0) {
      #   t1 <- geno[, 2*ind-1]
      #   t2 <- geno[, 2*ind]
      #   return(cbind(t1, t2))
      # }

      # find swap range in the chromosome
      swap.rg <- do.call(rbind, lapply(1:nrow(blk.rg), function(blk) {
        num.swap <- ifelse(recom.spot[blk], sample(range.hot, 1), sample(range.cold, 1))
        spot.swap.raw <- sort(sample((blk.rg[blk, 1]+1):(blk.rg[blk, 2]-1), num.swap, replace = TRUE))
        spot.swap <- c(blk.rg[blk, 1], spot.swap.raw)
        # simplify the recombination process
        flag1.swap <- !(num.swap %% 2 == 0) # if the first block need swapping
        flag.swap <- rep(flag1.swap, (num.swap+1))
        flag.swap[seq(2, (num.swap+1), 2)] <- !flag1.swap
        flag.swap.raw <- flag.swap[-length(flag.swap)]
        swap.op <- spot.swap[flag.swap]
        swap.ed <- spot.swap.raw[flag.swap.raw]
        return(cbind(swap.op, swap.ed))
      }))

      geno.t <- geno[, (2*ind-1):(2*ind)]
      for (swap in 1:nrow(swap.rg)) {
        op <- swap.rg[swap, 1]
        ed <- swap.rg[swap, 2]
        geno.t[op:ed, ] <- geno.t[op:ed, 2:1]
      }

      return(geno.t)
    }# end geno.swap.ind
    outgeno <- do.call(cbind, lapply(1:num.ind, geno.swap))
    
    # mutation
    spot.total <- num.marker * incols * num.ind
    num.mut <- ceiling(spot.total * rate.mut)
    row.mut <- sample(1:num.marker, num.mut)
    col.mut <- sample(1:(incols*num.ind), num.mut)
    for (i in 1:length(row.mut)) {
      if (geno[row.mut[i], col.mut[i]] == 1) {
        outgeno[row.mut[i], col.mut[i]] <- 0
      } else {
        outgeno[row.mut[i], col.mut[i]] <- 1
      }
    }

  } else if ((!is.null(geno) & incols == 1)|(!is.null(geno) & is.null(recom.spot))) {
    logging.log(" Mutation on genotype matrix...\n", verbose = verbose)
    num.marker <- nrow(geno)
    num.ind <- ncol(geno) / incols
    outgeno <- geno
    # outgeno <- deepcopy(geno)

    # mutation
    spot.total <- num.marker * incols * num.ind
    num.mut <- ceiling(spot.total * rate.mut)
    row.mut <- sample(1:num.marker, num.mut)
    col.mut <- sample(1:(incols*num.ind), num.mut)
    for (i in 1:length(row.mut)) {
      if (geno[row.mut[i], col.mut[i]] == 1) {
        outgeno[row.mut[i], col.mut[i]] <- 0
      } else {
        outgeno[row.mut[i], col.mut[i]] <- 1
      }
    }

  } else {
    stop("Please input the correct genotype matrix!")
  }

  return(outgeno)
}

#' Get map containing block and recombination information
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param input.map map that should be input, the marker number should be consistent in both map file and genotype data
#' @param num.marker number of markers of genotype matrix
#' @param len.block length of every blocks
#'
#' @return a map containing block and recombination information
#' @export
#'
#' @examples
#' data(simdata)
#' nmrk <- nrow(input.map)
#' pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
#' str(pos.map)
check.map <- function(input.map = NULL, num.marker = NULL, len.block = 5e7) {
  if (is.null(input.map))
    stop("Please input a map file!")
  if (num.marker != nrow(input.map))
    stop("The number of markers should be equal between genotype file and map file!")
  if (!is.data.frame(input.map))
    input.map <- is.data.frame(input.map)
  chrs <- unique(input.map[, 2])
  pos.mrk <- as.numeric(input.map[, 3])
  nc.map <- ncol(input.map)
  if (nc.map == 7) {
    map <- input.map
  } else {
    temp.recom <- lapply(chrs, function(chr) {
    f <- input.map[, 2] == chr
    if (nc.map == 5) {
      block <- (pos.mrk[f] %/% len.block) + 1
    } else if (nc.map == 6) {
      block <- input.map[f, 6]
    } else {
      stop("Please check the format of map!")
    }
    ub <- unique(block)
    tb <- table(block)
    if (tb[length(tb)] < 0.5*mean(tb[1:(length(tb)-1)])) {
      tb[length(ub)-1] <- tb[length(tb)-1] + tb[length(tb)]
      tb <- tb[-length(tb)]
      block[block == ub[length(ub)]] <- ub[length(ub)-1]
    }
    s1 <- length(tb)
    s2 <- s1 %/% 3
    r1 <- rep(c(1, 0, 1), c(s2, (s1-2*s2), s2))
    recom <- rep(r1, tb)
    t <- cbind(block, recom)
    return(t)
    })
    recom.spot <- do.call(rbind, temp.recom)
    map <- cbind(input.map[, 1:5], recom.spot)
  }

  return(map)
}


#' Calculate for block ranges in map
#'
#' Build date: Aug 15, 2019
#' Last update: Aug 15, 2019
#'
#' @author Dong Yin
#'
#' @param pos.map map with block information and recombination information
#'
#' @return block ranges
#' @export
#'
#' @examples
#' # get map with block and recombination information
#' data(simdata)
#' nmrk <- nrow(input.map)
#' pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
#'
#' # calculate for block ranges
#' blk.rg <- cal.blk(pos.map)
#' dim(blk.rg)
#' head(blk.rg)
cal.blk <- function(pos.map) {
  chr.uni <- unique(pos.map[, 2])
  chr.tab <- sapply(1:length(chr.uni), function(chr) { return(sum(pos.map[, 2] == chr.uni[chr])) })
  chr.ed <- sapply(1:length(chr.tab), function(chr) { return(sum(chr.tab[1:chr])) })
  chr.op <- chr.ed - chr.tab + 1
  blk.rg <- do.call(rbind, lapply(1:length(chr.tab), function(chr) {
    blk.tab <- table(as.numeric(pos.map[chr.op[chr]:chr.ed[chr], 6]))
    blk.ed <- sapply(1:length(blk.tab), function(blk) { return(sum(blk.tab[1:blk])) })
    blk.op <- blk.ed - blk.tab + 1
    blk.rg <- cbind(blk.op, blk.ed)
    return(blk.rg + chr.op[chr] - 1)
  }))
  return(blk.rg)
}

#' Input genotype column by column if markers are dense
#'
#' Build date: Nov 14, 2018
#' Last update: Jul 30, 2019
#'
#' @author Dong Yin
#'
#' @param bigmtr total genotype matrix
#' @param mtr genotype matrix should be inputting
#' @param ed index of the last column in each process
#' @param mrk.dense whether markers are dense, it is TRUE when sequencing data
#'
#' @return none
#' @export
#'
#' @examples
#' bigmtr <- bigmemory::big.matrix(nrow = 1e4, ncol = 2e2, type = 'char')
#' bigmtr[1:5, 1:5]
#' options(bigmemory.typecast.warning=FALSE)
#' mtr <- matrix(0, 1e4, 1e2)
#' input.geno(bigmtr = bigmtr, mtr = mtr, ed = ncol(mtr), mrk.dense = FALSE)
#' bigmtr[1:5, 1:5]
input.geno <- function(bigmtr, mtr, ed, mrk.dense) {
  op <- ed + 1 - ncol(mtr)
  if (mrk.dense) {
    for (i in 1:ncol(mtr)) {
      bigmtr[, op+i-1] <- mtr[, i]
    }
  } else {
    bigmtr[, op:ed] <- mtr[]
  }
}


#' Get true position in the genotype matrix
#'
#' Build date: June 11, 2020
#' Last update: June 11, 2020
#'
#' @param index the position of individual in the genotype matrix
#' @param incols the column number of an individual in the input genotype matrix, it can be 1 or 2
#'
#' @return sub-genotype matrix
#' @export
#'
#' @examples
#' index <- c(1:2, 5:6)
#' gmt <- getgmt(index = index, incols = 2)
getgmt <- function(index, incols = 2) {
  if (incols == 2) {
    gmt.dam <- index * 2
    gmt.sir <- gmt.dam - 1
    gmt.comb <- c(gmt.sir, gmt.dam)
    gmt.comb[seq(1, length(gmt.comb), 2)] <- gmt.sir
    gmt.comb[seq(2, length(gmt.comb), 2)] <- gmt.dam
  } else if (incols == 1) {
    gmt.comb <- index
  } else {
    stop("Please input a correct incols!")
  }
  
  return(gmt.comb)
}
