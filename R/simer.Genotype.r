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


#' Genotype simulation
#'
#' Generating and editing genotype data.
#' 
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
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
#' \item{$geno$pop.geno}{the genotype data.}
#' \item{$geno$incols}{'1': one-column genotype represents an individual; '2': two-column genotype represents an individual.}
#' \item{$geno$pop.marker}{the number of markers.}
#' \item{$geno$pop.ind}{the number of individuals in the base population.}
#' \item{$geno$prob}{the genotype code probability.}
#' \item{$geno$rate.mut}{the mutation rate of the genotype data.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate genotype simulation parameters
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)
#' 
#' # Run genotype simulation
#' SP <- genotype(SP)
genotype <- function(SP = NULL, ncpus = 0, verbose = TRUE) {

### Start genotype simulation
  
  # genotype parameters
  if (is.matrix(SP$geno$pop.geno) | is.big.matrix(SP$geno$pop.geno)) {
    pop.geno <- SP$geno$pop.geno
    SP$geno$pop.geno <- list(1)
  } else if (is.data.frame(SP$geno$pop.geno)) {
    pop.geno <- as.matrix(SP$geno$pop.geno)
    SP$geno$pop.geno <- list(1)
  } else {
    pop.geno <- SP$geno$pop.geno[[length(SP$geno$pop.geno)]]
  }
  pop.map <- SP$map$pop.map
  incols <- SP$geno$incols
  pop.marker <- SP$geno$pop.marker
  pop.ind <- SP$geno$pop.ind
  prob <- SP$geno$prob
  rate.mut <- SP$geno$rate.mut
  
  if (is.null(SP)) {
    stop("'SP' should be specified!")
  }
  if (is.null(pop.geno) & is.null(pop.marker) & is.null(pop.ind)) {
    stop("Please input information of genotype!")
  }

  if (!is.null(pop.geno)) {
    logging.log(" Input outer genotype matrix...\n", verbose = verbose)
    if (is.big.matrix(pop.geno)) {
      bigmat <- pop.geno
    } else {
      bigmat <- big.matrix(
        nrow = nrow(pop.geno),
        ncol = ncol(pop.geno),
        init = 3,
        type = 'char')
      Mat2BigMat(bigmat@address, mat = pop.geno, threads = ncpus)
    }

  } else if (!is.null(pop.marker) & !is.null(pop.ind)) {
    logging.log(" Establish genotype matrix of base-population...\n", verbose = verbose)
    if (incols == 2) {
      codes <- c(0, 1)
      if (!is.null(prob) & length(prob) != 2) {
        stop("The length of prob should be 2!")
      }
    } else if (incols == 1) {
      codes <- c(0, 1, 2)
      if (!is.null(prob) & length(prob) != 3) {
        stop("The length of prob should be 3!")
      }
    } else {
      stop("'incols' should only be 1 or 2!")
    }
    pop.geno <- matrix(sample(codes, pop.marker*pop.ind*incols, prob = prob, replace = TRUE), pop.marker, incols*pop.ind)
    bigmat <- big.matrix(
      nrow = nrow(pop.geno),
      ncol = ncol(pop.geno),
      init = 3,
      type = 'char')
    Mat2BigMat(bigmat@address, mat = pop.geno, threads = ncpus)
    
  } else {
    stop("Please input the correct genotype matrix!")
  }
  rm(pop.geno); gc()
  
  SP$geno$pop.marker <- pop.marker <- nrow(bigmat)
  SP$geno$pop.ind <- pop.ind <- ncol(bigmat) / incols
  
  if (!is.null(pop.map)) {
    if (nrow(bigmat) != nrow(pop.map)) {
      stop("Marker number should be same in both 'pop.map' and 'pop.geno'!")
    }
    Recom <- pop.map$Recom
    if (!is.null(Recom) & incols == 2) {
      # logging.log(" Chromosome exchange on genotype matrix...\n", verbose = verbose)
      Recom <- which(Recom %% 2 == 1)
      ind.swap <- sample(c(0, 1), pop.ind, replace = TRUE)
      ind.swap <- which(ind.swap == 1)
      for (ind in ind.swap) {
        geno.swap <- bigmat[Recom, (2*ind)]
        bigmat[Recom, (2*ind)] <- bigmat[Recom, (2*ind-1)]
        bigmat[Recom, (2*ind-1)] <- geno.swap
      }
    }
  }
      
  if (!is.null(rate.mut)) {
    # logging.log(" Mutation on genotype matrix...\n", verbose = verbose)
    spot.total <- pop.marker * incols * pop.ind
    num.mut <- ceiling(spot.total * rate.mut)
    row.mut <- sample(1:pop.marker, num.mut)
    col.mut <- sample(1:(incols*pop.ind), num.mut)
    for (i in 1:length(row.mut)) {
      if (bigmat[row.mut[i], col.mut[i]] == 1) {
        bigmat[row.mut[i], col.mut[i]] <- 0
      } else {
        bigmat[row.mut[i], col.mut[i]] <- 1
      }
    }
  }
  
  if (is.null(SP$geno$pop.geno)) {
    SP$geno$pop.geno <- list(bigmat)
  } else {
    SP$geno$pop.geno[[length(SP$geno$pop.geno)]] <- bigmat
  }
  names(SP$geno$pop.geno)[length(SP$geno$pop.geno)] <- paste0("gen", length(SP$geno$pop.geno))
  return(SP)
}

#' Annotation simulation
#'
#' Generating a map with annotation information
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#' @param verbose whether to print detail.
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$map$pop.map}{the map data with annotation information.}
#' \item{$map$qtn.num}{the QTN number for (each group in) each trait.}
#' \item{$map$qtn.model}{the genetic model of QTN such as 'A + D'.}
#' \item{$map$qtn.dist}{the QTN distribution containing 'norm', 'geom', 'gamma' or 'beta'.}
#' \item{$map$qtn.sd}{the standard deviations for normal distribution.}
#' \item{$map$qtn.prob}{the probability of success for geometric distribution.}
#' \item{$map$qtn.shape}{the shape parameter for gamma distribution.}
#' \item{$map$qtn.scale}{the scale parameter for gamma distribution.}
#' \item{$map$qtn.shape1}{the shape1 parameter for beta distribution.}
#' \item{$map$qtn.shape2}{the shape2 parameter for beta distribution.}
#' \item{$map$qtn.ncp}{the ncp parameter for beta distribution.}
#' \item{$map$qtn.spot}{the QTN distribution probability in each block.}
#' \item{$map$len.block}{the block length.}
#' \item{$map$maf}{the maf threshold, markers less than this threshold will be exclude.}
#' \item{$map$recom.spot}{whether to generate recombination events.}
#' \item{$map$range.hot}{the recombination times range in the hot spot.}
#' \item{$map$range.cold}{the recombination times range in the cold spot.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate annotation simulation parameters
#' SP <- param.annot(qtn.num = list(tr1 = 10))
#' 
#' # Run annotation simulation
#' SP <- annotation(SP)
annotation <- function(SP, verbose = TRUE) {

### Start annotation simulation
  
  # annotation parameters
  pop.map <- SP$map$pop.map
  pop.geno <- SP$geno$pop.geno
  qtn.model <- SP$map$qtn.model
  qtn.index <- SP$map$qtn.index
  qtn.num <- SP$map$qtn.num
  qtn.dist <- SP$map$qtn.dist
  qtn.sd <- SP$map$qtn.sd
  qtn.prob <- SP$map$qtn.prob
  qtn.shape <- SP$map$qtn.shape
  qtn.scale <- SP$map$qtn.scale
  qtn.shape1 <- SP$map$qtn.shape1
  qtn.shape2 <- SP$map$qtn.shape2
  qtn.ncp <- SP$map$qtn.ncp
  qtn.spot <- SP$map$qtn.spot
  len.block <- SP$map$len.block
  maf <- SP$map$maf
  recom.spot <- SP$map$recom.spot
  range.hot <- SP$map$range.hot
  range.cold <- SP$map$range.cold
  
  if (is.null(pop.map)) {
    pop.map <- generate.map(pop.marker = 1e4)
  }
  
  if (!is.data.frame(pop.map)) {
    pop.map <- as.data.frame(pop.map)
  }
    
  if (!is.numeric(pop.map[, 3])) {
    pop.map[, 3] <- as.numeric(pop.map[, 3])
  }
  nc.map <- ncol(pop.map)
  
  if (nc.map < 5) {
    stop("A map should contain 'SNP', 'Chrom', 'BP', 'ALT', and 'REF'!")
  }
  
  if (is.null(pop.map$Block) & (recom.spot | qtn.spot)) {
    pop.map$Block <- pop.map[, 3] %/% len.block + 1
  }
  
  if (is.null(pop.map$Recom) & recom.spot) {
    chrs <- unique(pop.map[, 2])
    Recom <- rep(0, nrow(pop.map))
    for (i in 1:length(chrs)) {
      block.tab <- table(pop.map$Block[pop.map[, 2] == chrs[i]])
      nblock <- length(block.tab)
      sublock <- nblock %/% 3
      recom.chr <- rep(c(1, 0, 1), c(sublock, (nblock-2*sublock), sublock))
      for (j in 1:nblock) {
        recom.flag <- pop.map[, 2] == chrs[i] & pop.map$Block == j
        recom.sum <- sum(recom.flag)
        if (recom.sum != 0) {
          if (recom.chr[j] == 0) {
            recom.times <- sample(range.cold, 1)
          } else if (recom.chr[j] == 1) {
            recom.times <- sample(range.hot, 1)
          } else {
            stop("'recom.spot' only contains '0' and '1'!")
          } 
          recom.len <- recom.sum %/% 10 + 1
          recom.op <- sort(sample(1:(recom.sum - recom.len), recom.times, replace = TRUE))
          recom.ed <- recom.op + recom.len - 1
          recom.sub <- rep(0, recom.sum)
          for (k in 1:length(recom.op)) {
            recom.sub[recom.op[k]:recom.ed[k]] <- recom.sub[recom.op[k]:recom.ed[k]] + 1
          }
          Recom[recom.flag] <- recom.sub
        } # if (sum(recom.flag) != 0) {
      } # for (j in 1:nblock) {
    } # for (i in 1:length(chrs)) {
    pop.map$Recom <- Recom
  }
  
  if (is.null(pop.map$QTNProb) & qtn.spot) {
    nblock <- max(pop.map$Block)
    block.prob <- runif(nblock, 0, 1)
    pop.map$QTNProb <- block.prob[pop.map$Block]
  }
  
  if (!is.null(pop.map$QTNProb) & !is.null(maf)) {
    if (is.null(pop.geno)) {
      stop("MAF calculation need genotype data!")
    }
    if (nrow(pop.geno) != nrow(pop.map)) {
      stop("Marker number should be same in both 'pop.map' and 'pop.geno'!")
    }
    MAF <- rowSums(pop.geno[]) / ncol(pop.geno)
    MAF <- pmin(MAF, 1 - MAF)
    pop.map$QTNProb[MAF < maf] <- 0
    pop.map$MAF <- MAF
  }
  
  # select markers as QTNs
  if (is.null(qtn.index)) {
    nTrait <- length(qtn.num)
  } else {
    nTrait <- length(qtn.index)
  }
  qtn.trn <- rep(list(NULL), nTrait)
  names(qtn.trn) <- paste0('tr', 1:nTrait)
  
  # QTN number
  if (is.null(qtn.index)) {
    for (i in 1:nTrait) {
      if (sum(qtn.num[[i]]) <= 0) { next  }
      qtn.trn[[i]] <- sort(sample(1:nrow(pop.map), sum(qtn.num[[i]])))
      logging.log(" Number of selected markers of trait", i, ":", qtn.num[[i]], "\n", verbose = verbose)
    }
    SP$map$qtn.index <- qtn.trn
    
  } else {
    qtn.trn <- qtn.index
    SP$map$qtn.num <- qtn.num <- lapply(qtn.trn, length)
  }
  
  # QTN effect
  qtn.model <- toupper(qtn.model)
  qtn.model <- sort(unique(unlist(strsplit(qtn.model, split = "\\s\\+\\s"))))
  qtn.model.single <- qtn.model[nchar(qtn.model) == 1]
  nAD <- length(qtn.model.single)
  qtn.trn.eff <- as.data.frame(matrix(NA, nrow(pop.map), nAD * nTrait))
  qtn.eff.name <- expand.grid(qtn.model.single, paste0("QTN", 1:nTrait))
  names(qtn.trn.eff) <- paste0(qtn.eff.name[, 2], "_", qtn.eff.name[, 1])
  for (i in 1:nTrait) {
    for (j in 1:nAD) {
      qtn.trn.eff[qtn.trn[[i]], nAD*(i-1) + j] <- cal.eff(qtn.num[[i]], qtn.dist[[i]], qtn.sd[[i]], qtn.prob[[i]], qtn.shape[[i]], qtn.scale[[i]], qtn.shape1[[i]], qtn.shape2[[i]], qtn.ncp[[i]])
    }
  }
  pop.map <- cbind(pop.map, qtn.trn.eff)
  
  # QTN interaction effect
  pop.map.GxG <- NULL
  if (any(nchar(qtn.model) > 1)) {
    qtn.trn.inteff <- rep(list(NULL), nTrait)
    for (i in 1:nTrait) {
      GxG.tmp <- GxG.network(pop.map, qtn.trn[[i]], qtn.model)
      GxG.tmp.eff <- cal.eff(length(GxG.tmp), qtn.dist[[i]], qtn.sd[[i]], qtn.prob[[i]], qtn.shape[[i]], qtn.scale[[i]], qtn.shape1[[i]], qtn.shape2[[i]], qtn.ncp[[i]])
      GxG.tmp <- data.frame(GxG.tmp, GxG.tmp.eff)
      names(GxG.tmp) <- c("GxG_name", paste0("GxG_eff", i))
      qtn.trn.inteff[[i]] <- GxG.tmp
    }
    pop.map.GxG <- Reduce(function(x, y) merge(x, y, by = "GxG_name", all = TRUE), qtn.trn.inteff, accumulate = FALSE)
  }
  
  if (is.null(pop.map.GxG)) {
    SP$map <- c(list(pop.map = pop.map), SP$map)
  } else {
    SP$map <- c(list(pop.map = pop.map, pop.map.GxG = pop.map.GxG), SP$map)
  }
  return(SP)
}

#' QTN genetic effects
#' 
#' Calculate for genetic effects vector of selected markers.
#' 
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param qtn.num integer: the QTN number of single trait; vector: the multiple group QTN number of single trait; matrix: the QTN number of multiple traits.
#' @param qtn.dist the QTN distribution containing 'norm', 'geom', 'gamma' or 'beta'.
#' @param qtn.sd the standard deviations for normal distribution.
#' @param qtn.prob the probability of success for geometric distribution.
#' @param qtn.shape the shape parameter for gamma distribution.
#' @param qtn.scale the scale parameter for gamma distribution.
#' @param qtn.shape1 the shape1 parameter for beta distribution.
#' @param qtn.shape2 the shape2 parameter for beta distribution.
#' @param qtn.ncp the ncp parameter for beta distribution.
#' 
#' @return a vector of genetic effect.
#' 
#' @export
#'
#' @examples
#' eff <- cal.eff(qtn.num = 10)
#' str(eff)
cal.eff <- function(qtn.num = 10, qtn.dist = "norm", qtn.sd = 1, qtn.prob = 0.5, qtn.shape = 1, qtn.scale = 1, qtn.shape1 = 1, qtn.shape2 = 1, qtn.ncp = 0) {

  if (sum(qtn.num) == 0) return(0)
  
  qtn.eff <- NULL
  for (nq in 1:length(qtn.num)) {
    if (qtn.dist[nq] == "norm") {
      qtn.eff <- c(qtn.eff, rnorm(qtn.num[nq], 0, qtn.sd[nq]))
    } else if (qtn.dist[nq] == "geom") {
      qtn.eff <- c(qtn.eff, rgeom(qtn.num[nq], qtn.prob[nq]))
    } else if (qtn.dist[nq] == "gamma") {
      qtn.eff <- c(qtn.eff, rgamma(qtn.num[nq], qtn.shape[nq], qtn.scale[nq]))
    } else if (qtn.dist[nq] == "beta") {
      qtn.eff <- c(qtn.eff, rbeta(qtn.num[nq], qtn.shape1[nq], qtn.shape2[nq], qtn.ncp[nq]))
    } else {
      stop("Please input a right QTN effect!")
    }
  }
  
  return(qtn.eff)
}

#' Genetic interaction network
#'
#' Generate genetic interaction effect combination network.
#' 
#' Build date: Mar 19, 2022
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param pop.map the map data with annotation information.
#' @param qtn.pos the index of QTNs in the map data.
#' @param qtn.model the genetic model of QTN such as 'A:D'.
#'
#' @return a data frame of genetic interaction effect.
#' 
#' @export
#'
#' @examples
#' pop.map <- generate.map(pop.marker = 1e4)
#' GxG.net <- GxG.network(pop.map)
#' head(GxG.net)
GxG.network <- function(pop.map = NULL, qtn.pos = 1:10, qtn.model = "A:D") {
  
  if (is.null(pop.map)) {
    stop("'pop.map' is necessary!")
  }
  
  qtn.model <- sort(unique(unlist(strsplit(qtn.model, split = "\\s\\+\\s"))))
  qtn.model.GxG <- qtn.model[nchar(qtn.model) > 1]
  qtn.model.GxG <- strsplit(qtn.model.GxG, split = ":")
  max.lev <- max(unlist(lapply(qtn.model.GxG, length)))
  if (length(qtn.model.GxG) == 0) {
    return(NULL)
  }
  
  pop.map.GxG <- NULL
  SNP <- pop.map[, 1]
  paths <- rep(list(NULL), max.lev)
  nQTNs <- length(qtn.pos)
  g <- ba.game(n = nQTNs, m = 2, directed = FALSE)
  for (i in 1:(nQTNs-1)) {
    for (j in (i+1):nQTNs) {
      ps <- all_simple_paths(g, from = i, to = j)
      for (k in 1:length(ps)) {
        if (length(ps[[k]]) <= max.lev) {
          ps_str <- paste(SNP[qtn.pos[ps[[k]]]], collapse = "-")
          ps_str <- paste(ps_str, qtn.model[nchar(qtn.model) == length(ps[[k]]) * 2 - 1], sep = "_")
          paths[[length(ps[[k]])]] <- c(paths[[length(ps[[k]])]], ps_str)
        } # end if (length(ps[[k]] <= max.lev)) {
      } # end for (k in 1:length(ps)) {
    } # end for (j in (i+1):nQTNs) {
  } # end for (i in 1:(nQTNs-1)) {
  
  comb_int <- unlist(paths)
  return(comb_int)
}

#' Marker information
#' 
#' Generate map data with marker information.
#' 
#' Build date: Mar 19, 2022
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param pop.marker the number of markers.
#' @param num.chr the number of chromosomes.
#' @param len.chr the length of chromosomes.
#'
#' @return a data frame with marker information.
#' 
#' @export
#'
#' @examples
#' pop.map <- generate.map(pop.marker = 1e4)
#' str(pop.map)
generate.map <- function(pop.marker = NULL, num.chr = 18, len.chr = 1.5e8) {
  
  if(is.null(pop.marker)) {
    stop("Please specify the number of markers!")
  }
  
  num.every <- rep(pop.marker %/% num.chr, num.chr)
  num.every[num.chr] <- num.every[num.chr] + pop.marker %% num.chr
  
  SNP <- paste("M", 1:(pop.marker), sep = "")
  Chrom <- rep(1:num.chr, num.every)
  BP <- do.call('c', lapply(1:num.chr, function(chr) {
    return(sort(sample(1:len.chr, num.every[chr]))) 
  }))
  
  base <- c("A", "T", "C", "G")
  ALT <- sample(base, pop.marker, replace = TRUE)
  REF <- sample(base, pop.marker, replace = TRUE)
  ff <- ALT == REF
  while (sum(ff) > 0) {
    REF[ff] <- sample(base, sum(ff), replace = TRUE)
    ff <- ALT == REF
  }
  
  map <- data.frame(SNP, Chrom, BP, ALT, REF)
  return(map)
}

#' Genotype code convertor 1
#' 
#' Convert genotype matrix from (0, 1) to (0, 1, 2).
#'
#' Build date: Nov 14, 2018
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param pop.geno genotype matrix of (0, 1).
#'
#' @return genotype matrix of (0, 1, 2).
#' 
#' @export
#'
#' @examples
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2, incols = 2)
#' SP <- genotype(SP)
#' geno1 <- SP$geno$pop.geno$gen1
#' geno2 <- geno.cvt1(geno1)
#' geno1[1:6, 1:4]
#' geno2[1:6, 1:2]
geno.cvt1 <- function(pop.geno) {
  if (is.null(pop.geno)) return(NULL)
  num.ind <- ncol(pop.geno) / 2
  v.odd <- (1:num.ind) * 2 - 1
  v.even <- (1:num.ind) * 2
  geno <- pop.geno[, v.odd] + pop.geno[, v.even]
  return(geno)
}

#' Genotype code convertor 2
#' 
#' Convert genotype matrix from (0, 1, 2) to (0, 1).
#'
#' Build date: Jul 11, 2020
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param pop.geno genotype matrix of (0, 1, 2).
#' 
#' @return genotype matrix of (0, 1).
#' 
#' @export
#'
#' @examples
#' SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2, incols = 1)
#' SP <- genotype(SP)
#' geno1 <- SP$geno$pop.geno$gen1
#' geno2 <- geno.cvt2(geno1)
#' geno1[1:6, 1:2]
#' geno2[1:6, 1:4]
geno.cvt2 <- function(pop.geno) {
  if (is.null(pop.geno)) return(NULL)
  nind <- ncol(pop.geno)
  nmrk <- nrow(pop.geno)
  geno <- matrix(3, nmrk, 2*nind)
  
  v.odd <- (1:nind) * 2 - 1
  v.even <- (1:nind) * 2
  
  for (i in 1:nind) {
    tmp1 <- pop.geno[, i]
    tmp1[tmp1 == 2] <- 1
    geno[, v.even[i]] <- tmp1
    tmp2 <- pop.geno[, i] - 1
    tmp2[tmp2 == -1] <- 0
    geno[, v.odd[i]] <- tmp2
  }
  
  return(geno)
}
