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
#' @param rawgeno genotype matrix from outside
#' @param geno genotype matrix need dealing with
#' @param num.marker number of markers
#' @param num.ind number of individuals
#' @param prob weight of "0" and "1" in genotype matrix, the sum of element in vector equals 1
#' @param blk.rg represents the started and ended position blocks
#' @param recom.spot whether to consider recombination in every blocks
#' @param range.hot range of exchages in hot spot block
#' @param range.cold range of exchages in cold spot block
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
#' basepop.geno <- genotype(num.marker = num.marker, num.ind = num.ind, verbose = TRUE)
#' basepop.geno[1:5, 1:5]
#'
#' # get block information and recombination information
#' nmrk <- nrow(basepop.geno)
#' nind <- ncol(basepop.geno) / 2
#' pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
#' blk.rg <- cal.blk(pos.map)
#' recom.spot <- as.numeric(pos.map[blk.rg[, 1], 7])
#'
#' # genotype matrix after Exchange and Mutation
#' basepop.geno.em <-
#' genotype(geno = basepop.geno,
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
    logging.log("Input outer genotype matrix...\n", verbose = verbose)
    if (!is.matrix(rawgeno)) {
      rawgeno <- as.matrix(rawgeno)
    }
    outgeno <- rawgeno
    num.marker <- nrow(rawgeno)
    num.ind <- ncol(rawgeno) / 2

  } else if (!is.null(num.marker) && !is.null(num.ind)){
    logging.log("Establish genotype matrix of base-population...\n", verbose = verbose)
    outgeno <- matrix(sample(c(0, 1), num.marker*num.ind*2, prob = prob, replace = TRUE), num.marker, 2*num.ind)

  } else if (!is.null(geno) & !is.null(recom.spot)) {
    logging.log("Chromosome exchange and mutation on genotype matrix...\n", verbose = verbose)
    num.marker <- nrow(geno)
    num.ind <- ncol(geno) / 2
    # outgeno <- deepcopy(geno) # deepcopy() in bigmemory
    ind.swap <- sample(c(0, 1), num.ind, replace = TRUE)

    geno.swap <- function(ind) {
      if (ind.swap[ind] == 0) {
        t1 <- geno[, 2*ind-1]
        t2 <- geno[, 2*ind]
        return(cbind(t1, t2))
      }

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
    spot.total <- num.marker * 2 * num.ind
    num.mut <- ceiling(spot.total * rate.mut)
    row.mut <- sample(1:num.marker, num.mut)
    col.mut <- sample(1:(2*num.ind), num.mut)
    for (i in 1:length(row.mut)) {
      if (geno[row.mut[i], col.mut[i]] == 1) {
        outgeno[row.mut[i], col.mut[i]] <- 0
      } else {
        outgeno[row.mut[i], col.mut[i]] <- 1
      }
    }

  } else if (!is.null(geno) & is.null(recom.spot)) {
    logging.log("Mutation on genotype matrix...\n", verbose = verbose)
    num.marker <- nrow(geno)
    num.ind <- ncol(geno) / 2
    outgeno <- geno
    # outgeno <- deepcopy(geno)

    # mutation
    spot.total <- num.marker * 2 * num.ind
    num.mut <- ceiling(spot.total * rate.mut)
    row.mut <- sample(1:num.marker, num.mut)
    col.mut <- sample(1:(2*num.ind), num.mut)
    for (i in 1:length(row.mut)) {
      if (geno[row.mut[i], col.mut[i]] == 1) {
        outgeno[row.mut[i], col.mut[i]] <- 0
      } else {
        outgeno[row.mut[i], col.mut[i]] <- 1
      }
    }

  } else {
    stop("please input the correct genotype matrix!")
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
#' @param input.map map from outside
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
    stop("please input a map file!")
  if (num.marker != nrow(input.map))
    stop("the number of markers should be equal between genotype file and map file!")
  if (!is.matrix(input.map))
    input.map <- as.matrix(input.map)
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
#' @param mrk.dense whether markers are dense
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

# By Haohao Zhang
#' Remove big.matrix safely
#'
#' Build date: Aug 8, 2019
#' Last update: Aug 18, 2019
#'
#' @author Haohao Zhang and Dong Yin
#'
#' @param x filename of big.matrix
#' @param desc_suffix suffix of description file of big.matrix
#' @param bin_suffix suffix of binary file of big.matrix
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' library(bigmemory)
#' mat <- filebacked.big.matrix(
#'      nrow = 10,
#'      ncol = 10,
#'      init = 0,
#'      type = 'char',
#'      backingpath = ".",
#'      backingfile = 'genotype.geno.bin',
#'      descriptorfile = 'genotype.geno.desc')
#' options(bigmemory.typecast.warning=FALSE)
#'
#' remove_bigmatrix(x = "genotype")
remove_bigmatrix <- function(x, desc_suffix=".geno.desc", bin_suffix=".geno.bin") {
  name <- basename(x)
  path <- dirname(x)

  descfile <- paste0(name, desc_suffix)
  binfile  <- paste0(name, bin_suffix)

  remove_var <- function(binfile, envir) {
    for (v in ls(envir = envir)) {
      if (class(get(v, envir = envir)) == "big.matrix") {
        desc <- describe(get(v, envir = envir))@description
        if (desc$filename == binfile) { # TODO: Risky deletion
          rm(list = v, envir = envir)
          gc()
        }
      } else if (class(get(v, envir = envir)) == "list") {
        if (is.null(get(v, envir = envir)$geno)) next
        if (class(get(v, envir = envir)$geno) == "big.matrix") {
          desc <- describe(get(v, envir = envir)$geno)@description
          if (desc$filename == binfile) { # TODO: Risky deletion
            rm(list = v, envir = envir, inherits = TRUE)
            gc()
          }
        }
      }
    }
  }

  # Delete objects that occupy binfile in the global environment
  remove_var(binfile, as.environment(-1L))
  remove_var(binfile, globalenv())

  if (file.exists(descfile)) {
    file.remove(descfile)
  }
  if (file.exists(binfile)) {
    file.remove(binfile)
  }
}