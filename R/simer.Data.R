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


#' Simer.Data: To prepare data for Simer package
#'
#' Build date: May 26, 2021
#' Last update: Aug 23, 2021
#'
#' @author Dong Yin
#'
#' @param fileMVP Genotype in MVP format
#' @param fileBed Genotype in PLINK binary format
#' @param filePlinkPed Genotype in PLINK ped format
#' @param genoType type parameter in bigmemory, genotype data. The default is char, it is highly recommended *NOT* to modify this parameter.
#' @param filterGeno threshold of sample miss rate
#' @param filterHWE threshold of Hardy-Weinberg Test
#' @param filterMind threshold of variant miss rate 
#' @param filterMAF threshold of Minor Allele Frequency
#' @param filePed Pedigree, 3-columns or 15 columns pedigree
#' @param standardID whether kid id is 15-character standard
#' @param fileSir the file name of candidate sires
#' @param fileDam the file name of candidate dams
#' @param exclThres if conflict ratio is more than exclThres, exclude this parent
#' @param assignThres if conflict ratio is less than assignThres, assign this parent to the individual
#' @param pedSep the separator of pedigree file
#' @param filePhe Phenotype, the first column is taxa name, the subsequent columns are traits
#' @param planPhe breeding plans about phenotype
#' @param pheCols the column neeeding extracting
#' @param pheSep seperator for phenotype file.
#' @param missing the missing value
#' @param SNP.impute "Left", "Middle", "Right", or NULL for skip impute.
#' @param outpath the output path
#' @param out prefix of output file name
#' @param maxLine number of SNPs, only used for saving memory when calculate kinship matrix
#' @param priority "speed" or "memory"
#' @param ncpus The number of threads used, if NULL, (logical core number - 1) is automatically used
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return NULL
#'
#' @examples
#' options(simer.OutputLog2File = FALSE)
#' filePhe <- system.file("extdata", "phenotype.txt", package = "simer")
#' filePhe <- simer.Data.Pheno(filePhe = filePhe, out = tempfile("outfile"))
simer.Data <- function(fileMVP = NULL, fileBed = NULL, filePlinkPed = NULL, genoType = 'char', filterGeno=NULL, filterHWE=NULL, filterMind=NULL, filterMAF=NULL,
                       filePed = NULL, standardID = FALSE, fileSir=NULL, fileDam=NULL, exclThres=0.01, assignThres=0.005, pedSep='\t', 
                       filePhe = NULL, planPhe = NULL, pheCols = NULL, pheSep = '\t', missing = c(NA, 'NA', '-9', 9999),
                       SNP.impute = "Major",
                       outpath = getwd(), out = 'simer', maxLine = 10000, priority = "speed", ncpus = 0, verbose = TRUE) {
  
  logging.initialize("Simer.Data", outpath)
  
  if (!is.null(out)) {  out <- file.path(outpath, out) }
  
  if (!is.null(filePlinkPed)) {
    system(paste0("plink --file ", filePlinkPed, " --make-bed --out ", filePlinkPed))
    fileBed <- filePlinkPed
  }
  
  if (!is.null(fileBed)) {
    logging.log("******************** Data Format Convert ********************\n", verbose = verbose)
    if (is.null(out)) { out <- paste0(fileBed, ".qc") }
    MVP.Data.Bfile2MVP(
      bfile = fileBed, 
      out = out, 
      maxLine = maxLine, 
      priority = priority, 
      type.geno = genoType,
      verbose = verbose,
      threads = ncpus
    )
    fileMVP <- out
  }
  
  genoFileName <- NULL
  if (!is.null(fileMVP)) {
    logging.log("*************** Genotype Data Quality Control ***************\n", verbose = verbose)
    genoFileName <-
      simer.Data.Geno(
        fileMVP = fileMVP,
        filePed = filePed,
        out = out,
        genoType = genoType,
        filterGeno = filterGeno,
        filterHWE = filterHWE,
        filterMind = filterMind,
        filterMAF = filterMAF,
        ncpus = ncpus,
        verbose = verbose)
  }
  
  pedFileName <- NULL
  if (!is.null(filePed)) {
    logging.log("*************** Pedigree Data Quality Control ***************\n", verbose = verbose)
    pedFileName <- 
      simer.Data.Ped(
        filePed = filePed, 
        fileMVP = NULL,
        out = out, 
        standardID = standardID, 
        fileSir = fileSir, 
        fileDam = fileDam, 
        exclThres = exclThres, 
        assignThres = assignThres, 
        sep = pedSep, 
        ncpus = ncpus, 
        verbose = verbose
    )
  }
  
  if (!is.null(filePhe)) {
    logging.log("*************** Phenotype Data Quality Control **************\n", verbose = verbose)
    pheFileName <- 
      simer.Data.Pheno(
        filePhe = filePhe, 
        filePed = pedFileName,
        out = out, 
        planPhe = planPhe, 
        pheCols = pheCols, 
        sep = pheSep, 
        missing = missing, 
        verbose = verbose)
  }
   
}

#' simer.Data.MVP2MVP: Change the path of MVP files
#' Author: Dong Yin
#' Build date: May 26, 2021
#' Last update: May 26, 2021
#' 
#' @param fileMVP the prefix of MVP file
#' @param genoType type parameter in bigmemory, genotype data. The default is char, it is highly recommended *NOT* to modify this parameter.
#' @param out the name of output file
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' Output file:
#' <out>.geno.desc
#' <out>.geno.bin
#' <out>.geno.ind
#' <out>.geno.map
#' 
#' @examples
#' mvpPath <- system.file("extdata", "01bigmemory", "demo", package = "simer")
#' simer.Data.MVP2MVP(mvpPath, out = tempfile("outfile"))
simer.Data.MVP2MVP <- function(fileMVP, genoType='char', out='simer', verbose=TRUE) {
  t1 <- as.numeric(Sys.time())
  
  if (fileMVP != out) { remove_bigmatrix(out) }
  fileDesc <- normalizePath(paste0(fileMVP, '.geno.desc'), mustWork = TRUE)
  fileInd <- normalizePath(paste0(fileMVP, '.geno.ind'), mustWork = TRUE)
  fileMap <- normalizePath(paste0(fileMVP, '.geno.map'), mustWork = TRUE)
  
  backingfile <- paste0(basename(out), ".geno.bin")
  descriptorfile <- paste0(basename(out), ".geno.desc")
  
  bigmat <- attach.big.matrix(fileDesc)
  deepcopy(x = bigmat,
           type = genoType,
           backingfile = backingfile,
           backingpath =dirname(out),
           descriptorfile = descriptorfile,
  )
  file.copy(fileInd, paste0(out, ".geno.ind"))
  file.copy(fileMap, paste0(out, ".geno.map"))
  
  t2 <- as.numeric(Sys.time())
  logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(invisible(dim(bigmat)))
}

#' simer.Data.Impute: Impute the missing value
#' Author: Dong Yin
#' Build date: May 26, 2021
#' Last update: Sep 18, 2021
#' 
#' @param fileMVP Genotype in MVP format
#' @param fileBed Genotype in PLINK binary format
#' @param out the name of output file
#' @param maxLine number of SNPs, only used for saving memory when calculate kinship matrix
#' @param ncpus The number of threads used, if NULL, (logical core number - 1) is automatically used
#' @param verbose whether to print detail.
#'
#' @return genotype data after imputing
#' @export
#'
#' @examples
#' # need beagle
#' # fileMVP <- system.file("extdata", "01bigmemory", "demo", package = "simer")
#' # fileMVPimp <- simer.Data.Impute(fileMVP = fileMVP)
simer.Data.Impute <- function(fileMVP = NULL, fileBed = NULL, out = NULL, maxLine = 1e4, ncpus = 0, verbose = TRUE) {
  
  if (sum(is.null(fileMVP), is.null(fileBed)) != 1) {
    stop("Only a file type can be input!")
  }
  
  if (!is.null(fileMVP)) {
    descFile <- normalizePath(paste0(fileMVP, ".geno.desc"), mustWork = TRUE)
    bigmat <- attach.big.matrix(descFile)
    if (!hasNA(bigmat@address)) {
      message("No NA in genotype, imputation has been skipped.")
      return()
    } else {
      if (is.null(out)) { out <- paste0(fileMVP, ".imp") }
      mapFile <- normalizePath(paste0(fileMVP, ".geno.map"), mustWork = TRUE)
      map <- read.table(mapFile, header = TRUE)
      tmpout <- tempfile()
      MVP.Data.MVP2Bfile(bigmat = bigmat, map = map, out = tmpout)
      system(paste('plink --bfile', tmpout, "--recode vcf-iid --out", tmpout))
    }
  }
  
  if (!is.null(fileBed)) {
    famFile <- normalizePath(paste0(fileBed, '.fam'), mustWork = TRUE)
    bedFile <- normalizePath(paste0(fileBed, '.bed'), mustWork = TRUE)
    fam <- read.table(famFile, header = FALSE)
    n <- nrow(fam)
    hasNA <- hasNABed(bedFile, n, maxLine, ncpus, verbose)
    if (!hasNA) {
      message("No NA in genotype, imputation has been skipped.")
      return()
    } else {
      if (is.null(out)) { out <- paste0(fileBed, ".imp") }
      system(paste('plink --bfile', tmpout, "--recode vcf-iid --out", tmpout))
    }
  }
  
  system('if [ ! -f beagle.28Jun21.220.jar ]; then
              echo
              echo "Downloading beagle.28Jun21.220.jar"
              wget http://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar
              fi')
  
  system(paste0("java -jar beagle.28Jun21.220.jar gt=", tmpout, ".vcf out=", tmpout))
  
  system(paste0('plink --vcf ', tmpout, ".vcf.gz --make-bed --out ", tmpout))
  
  MVP.Data.Bfile2MVP(bfile = tmpout, out = out)
  
  return(out)
}

#' simer.Data.Geno: Data quality control of genotype data
#' Author: Dong Yin
#' Build date: June 3, 2021
#' Last update: Aug 22, 2021
#' 
#' @param fileMVP the prefix of MVP file
#' @param filePed the name of pedigree file
#' @param out the name of output file
#' @param genoType type parameter in bigmemory, genotype data. The default is char, it is highly recommended *NOT* to modify this parameter.
#' @param filterGeno threshold of sample miss rate
#' @param filterHWE threshold of Hardy-Weinberg Test
#' @param filterMind threshold of variant miss rate 
#' @param filterMAF threshold of Minor Allele Frequency
#' @param ncpus The number of threads used, if NULL, (logical core number - 1) is automatically used
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' Output file:
#' <out>.geno.desc
#' <out>.geno.bin
#' <out>.geno.ind
#' 
#' @examples
#' fileMVP <- system.file("extdata", "01bigmemory", "demo", package = "simer")
#' filePed <- system.file("extdata", "pedigree.txt", package = "simer")
#' simer.Data.Geno(fileMVP=fileMVP, filePed=filePed, genoType='char', 
#'   out=tempfile("outfile"), 
#'   filterGeno=0.1, filterHWE=0.001, filterMind=0.1, filterMAF=0.05)
simer.Data.Geno <- function(fileMVP, filePed=NULL, out='simer', genoType='char',
                            filterGeno=NULL, filterHWE=NULL, filterMind=NULL, filterMAF=NULL,
                            ncpus=0, verbose=TRUE) {
  
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Genotype Data.\n", verbose = verbose)
  
  if (length(fileMVP) == 0) { fileMVP <- NULL }
  if (length(filePed) == 0) { filePed <- NULL }
  if (is.null(out)) { out <- paste0(fileMVP, ".qc") }
  if (fileMVP != out) { remove_bigmatrix(out) }
  
  fileDesc <- normalizePath(paste0(fileMVP, '.geno.desc'), mustWork = TRUE)
  fileInd <- normalizePath(paste0(fileMVP, '.geno.ind'), mustWork = TRUE)
  fileMap <- normalizePath(paste0(fileMVP, '.geno.map'), mustWork = TRUE)
  
  genoInd <- read.table(fileInd, sep = '\t', header = FALSE)[, 1]
  genoMap <- read.table(fileMap, sep = '\t', header = TRUE)
  
  if (!is.null(filePed)) {
    ped <-  read.table(filePed, sep = '\t', header = TRUE)
    ped <- unique(unlist(ped))
    keepInds <- which(genoInd %in% ped)
  } else {
    keepInds <- NULL
  }
  
  backingfile <- paste0(basename(out), ".geno.bin")
  descriptorfile <- paste0(basename(out), ".geno.desc")
  
  bigmat <- attach.big.matrix(fileDesc)
  genoInfo <- GenoFilter(bigmat@address, keepInds, filterGeno, filterHWE, filterMind, filterMAF, ncpus, verbose)
  keepRows <- genoInfo$keepRows
  keepCols <- genoInfo$keepCols
  
  deepcopy(x = bigmat,
           cols = keepCols,
           rows = keepRows,
           type = genoType,
           backingfile = backingfile,
           backingpath =dirname(out),
           descriptorfile = descriptorfile,
  )
  
  genoInd <- genoInd[keepCols]
  genoMap <- genoMap[keepRows, ]
  write.table(genoInd, paste0(out, ".geno.ind"), quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
  write.table(genoMap, paste0(out, ".geno.map"), quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  t2 <- as.numeric(Sys.time())
  logging.log(" Preparation for GENOTYPE data for", basename(fileMVP), "is done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(out)
}

#' simer.Data.Ped: To check pedigree file
#' Author: LiLin Yin and Dong Yin
#' Build date: May 6, 2021
#' Last update: June 3, 2021
#' 
#' @param filePed the name of pedigree file need correcting
#' @param fileMVP Genotype in MVP format
#' @param out the name of output file
#' @param standardID whether kid id is 15-character standard
#' @param fileSir the file name of candidate sires
#' @param fileDam the file name of candidate dams
#' @param exclThres if conflict ratio is more than exclThres, exclude this parent
#' @param assignThres if conflict ratio is less than assignThres, assign this parent to the individual
#' @param header whether the file contains header
#' @param sep separator of the file
#' @param ncpus The number of threads used, if NULL, (logical core number - 1) is automatically used
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' Output file:
#' <out>.report.txt
#' <out>.error.txt
#' <out>.qc.txt
#' 
#' @examples
#' # simer.Data.Ped needs genotype data
#' out <- tempfile("outfile")
#' fileMVP <- system.file("extdata", "01bigmemory", "demo", package = "simer")
#' simer.Data.MVP2MVP(fileMVP, out = out)
#' 
#' filePed <- system.file("extdata", "pedigree.txt", package = "simer")
#' simer.Data.Ped(filePed = filePed, out = out, standardID = FALSE)
simer.Data.Ped <- function(filePed, fileMVP=NULL, out=NULL, standardID=FALSE, fileSir=NULL, fileDam=NULL, 
                           exclThres=0.01, assignThres=0.005, header=TRUE, sep='\t', ncpus=0, verbose=TRUE) {
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Pedigree Data.\n", verbose = verbose)
  # read data
  # if (!is.vector(filePed)) { filePed <- c(filePed) }
  if (length(filePed) == 0) { filePed <- NULL }
  
  # pedigree files
  if (is.data.frame(filePed)) {
    pedigree <- filePed
  } else {
    if (length(filePed) == 1) {
      pedigree <- read.table(filePed, sep = sep, header = header, stringsAsFactors = FALSE)
    } else {
      pedigree <- do.call(rbind, lapply(1:length(filePed), function(i) {
        return(read.table(filePed[i], sep = sep, header = header, stringsAsFactors = FALSE))
      }))
    }
  }
  
  if (length(fileSir) == 0) { fileSir <- NULL }
  if (length(fileDam) == 0) { fileDam <- NULL }
  if (!is.null(fileSir)) {
    candSir <- unlist(read.table(fileSir, stringsAsFactors = FALSE))
  } else {
    candSir <- NULL
  }
  if (!is.null(fileDam)) {
    candDam <- unlist(read.table(fileDam, stringsAsFactors = FALSE))
  } else {
    candDam <- NULL
  }

  # thanks for YinLL for sharing code
  # get 3-column pedigree
  if (!is.matrix(pedigree))	pedigree <- as.matrix(pedigree)
  pedigree <- apply(pedigree, 2, as.character)
  pedigree[pedigree == ""] <- "0"
  pedigree[is.na(pedigree)] <- "0"
  if (ncol(pedigree) != 3) {
    if (ncol(pedigree) != 15)
      stop("Please check your data! pegigree information in 15 columns which contain 3 generations' information are needed!")
    if (standardID) {
      id <- unique(c(pedigree))
      if (sum(nchar(id) != 15) > 0)
        stop(paste("The format of below individuals don't meet the requirements:","\n", id[nchar(id) != 15], sep = ""))
    }
    # Ind Sir  SS  SD SSS SSD SDS SDD Dam  DS  DD DSS DSD DDS DDD
    #   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    # pedx <- rbind(
    #   pedigree[, c(1, 2, 9)],
    #   pedigree[, c(2, 3, 4)],
    #   pedigree[, c(9, 10, 11)],
    #   pedigree[, c(3, 5, 6)],
    #   pedigree[, c(4, 7, 8)],
    #   pedigree[, c(10, 12, 13)],
    #   pedigree[, c(11, 14, 15)],
    #   cbind(pedigree[, 5], 0, 0),
    #   cbind(pedigree[, 6], 0, 0),
    #   cbind(pedigree[, 7], 0, 0),
    #   cbind(pedigree[, 8], 0, 0),
    #   cbind(pedigree[, 12], 0, 0),
    #   cbind(pedigree[, 13], 0, 0),
    #   cbind(pedigree[, 14], 0, 0),
    #   cbind(pedigree[, 15], 0, 0)
    # )
    # Ind Sir Dam  SS  SD  DS  DD SSS SSD SDS SDD DSS DSD DDS DDD
    #   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    pedx <- rbind(
      pedigree[, c(1, 2, 3)],
      pedigree[, c(2, 4, 5)],
      pedigree[, c(3, 6, 7)],
      pedigree[, c(4, 8, 9)],
      pedigree[, c(5, 10, 11)],
      pedigree[, c(6, 12, 13)],
      pedigree[, c(7, 14, 15)],
      cbind(pedigree[, 8], 0, 0),
      cbind(pedigree[, 9], 0, 0),
      cbind(pedigree[, 10], 0, 0),
      cbind(pedigree[, 11], 0, 0),
      cbind(pedigree[, 12], 0, 0),
      cbind(pedigree[, 13], 0, 0),
      cbind(pedigree[, 14], 0, 0),
      cbind(pedigree[, 15], 0, 0)
    )
    pedx <- pedx[!duplicated(pedx[, 1]), ]
  }else{
    pedx <- pedigree
    pedx0 <- setdiff(pedx[,c(2:3)], pedx[, 1])
    pedx0 <- pedx0[pedx0 != "0"]
    if(length(pedx0) > 0){
      pedx0 <- cbind(pedx0, "0", "0")
      colnames(pedx0) <- colnames(pedx)
      pedx <- rbind(pedx0, pedx)
    }
  }
  pedx <- pedx[pedx[, 1] != "0", ]
  pedError <- pedx[duplicated(pedx[, 1]), ]
  pedx <- pedx[!duplicated(pedx[, 1]), ]
  
  # read genotype data
  if (length(fileMVP) == 0) { fileMVP <- NULL }
  if (!is.null(fileMVP)) {
    hasGeno <- TRUE
    genoFile <- normalizePath(paste0(fileMVP, ".geno.desc"), mustWork = TRUE)
    genoIDFile <- normalizePath(paste0(fileMVP, ".geno.ind"), mustWork = TRUE)
    genoID <- as.character(unlist(read.table(genoIDFile, header = FALSE)))
    geno <- attach.big.matrix(genoFile)
  } else {
    hasGeno <- FALSE
  }
  
  if (standardID == TRUE) {
    pedx <- cbind(pedx, substr(pedx[, 1], 8, nchar(pedx[, 1])))
    index.born <- order(as.numeric(pedx[, 4]))
    pedx <- pedx[index.born, ]
    ped <- pedx[, 1:3]
    birthDate <- as.numeric(pedx[, 4])
    rm(pedx); gc()
    if (hasGeno) {
      ped <- PedigreeCorrector(geno@address, genoID, ped, candSir, candDam, exclThres, assignThres, birthDate, ncpus, verbose)
    }
    
  }else{
    if (hasGeno) {
      birthDate = NULL
      pedx <- PedigreeCorrector(geno@address, genoID, pedx, candSir, candDam, exclThres, assignThres, birthDate, ncpus, verbose)
    }
    
    #print("Making needed files")
    pedx1 <- pedx[pedx[, 2] == "0" & pedx[, 3] == "0", ]
    pedx2 <- pedx[!(pedx[, 2] == "0" & pedx[, 3] == "0"), ]
    go = TRUE
    while(go == TRUE) {
      Cpedx <- pedx1[, 1]
      index <- (pedx2[, 2] %in% Cpedx) & (pedx2[, 3] %in% Cpedx)
      if (sum(index) == 0) {
        index.sir <- pedx2[, 2] %in% Cpedx
        index.dam <- pedx2[, 3] %in% Cpedx
        if (sum(index.sir) != 0 | sum(index.dam) != 0) {
          # only one parent can be found
          pedError <- rbind(pedError, pedx2[index.sir | index.dam, 1:3])
          pedx2[index.sir, 3] <- "0"
          pedx2[index.dam, 2] <- "0"
          pedx1 <- rbind(pedx1, pedx2[index.sir | index.dam, ])
          pedx2 <- pedx2[!(index.sir | index.dam), ]
        } else {
          # no parent can be found
          pedx02 <- setdiff(pedx2[,c(2:3)], pedx2[, 1])
          pedx02 <- pedx02[pedx02 != "0"]
          if(length(pedx02) > 0){
            pedx02 <- cbind(pedx02, "0", "0")
            colnames(pedx02) <- colnames(pedx2)
            pedx1 <- rbind(pedx1, pedx02)
          } else {
            pedError <- rbind(pedError, pedx2[, 1:3])
            pedx2[, 2:3] <- "0"
            pedx1 <- rbind(pedx1, pedx2)
            pedx2 <- pedx2[-(1:nrow(pedx2)), ]
          }
        }
      } else {
        pedx1 <- rbind(pedx1, pedx2[index, ])
        pedx2 <- pedx2[!index, ]
      }
      if ("character" %in% class(pedx2)) pedx2 <- matrix(pedx2, 1)
      if (nrow(pedx2) == 0) go = FALSE
    }
    ped <- pedx1
    rm(pedx1);rm(pedx2);gc()
  }
  
  if (hasGeno) {
    pedError <- rbind(pedError, ped[ped$sirState=="NotFound" | ped$damState=="NotFound", 1:3])
    pedUse <- ped[, 1:3]
    pedUse$sir[ped$sirState == "NotFound"] <- "0"
    pedUse$dam[ped$damState == "NotFound"] <- "0"
  } else {
    pedError <- rbind(pedError, ped[ped[, 1]==ped[, 2] | ped[, 1] == ped[, 3], 1:3])
    pedUse <- ped[, 1:3]
    pedUse[pedUse[, 1] == pedUse[, 2], 2] <- "0"
    pedUse[pedUse[, 1] == pedUse[, 3], 3] <- "0"
  }
  
  if (is.null(out)) {
    out <- unlist(strsplit(filePed, split = '.', fixed = TRUE))[1]
  }
  write.table(ped, paste0(out, ".report.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  write.table(pedError, paste0(out, ".error.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  write.table(pedUse, paste0(out, ".qc.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  
  t2 <- as.numeric(Sys.time())
  logging.log(" Preparation for PEDIGREE data is done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(paste0(out, ".qc.txt"))
}

#' simer.Data.Pheno: Data quality control of phenotype data
#' Author: Haohao Zhang and Dong Yin
#' Build date: June 13, 2021
#' Last update: July 8, 2021
#'
#' @param filePhe the phenotype files, it can be a vector
#' @param filePed the pedigree files, it can be a vector
#' @param out the name prefix of output file
#' @param planPhe breeding plans about phenotype
#' @param pheCols the column needing extracting
#' @param header the header of file
#' @param sep the separator of file
#' @param missing the missing value
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' Output file:
#' <out>.qc.txt
#' 
#' @examples
#' filePhe <- system.file("extdata", "phenotype.txt", package = "simer")
#' simer.Data.Pheno(filePhe = filePhe, out = tempfile("outfile"))
simer.Data.Pheno <- function(filePhe, filePed=NULL, out=NULL, planPhe=NULL, pheCols=NULL, header=TRUE, sep='\t', missing=c(NA, 'NA', 'Na', '.', '-', 'NAN', 'nan', 'na', 'N/A', 'n/a', '<NA>', '', '-9', 9999), verbose=TRUE) {
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Phenotype Data.\n", verbose = verbose)
  # read data
  # if (!is.vector(filePhe)) { filePhe <- c(filePhe) }
  if (length(filePhe) == 0) { filePhe <- NULL  }
  if (length(filePed) == 0) { filePed <- NULL }
  
  if (!is.null(filePed)) {
    ped <-  read.table(filePed, sep = '\t', header = TRUE)
    ped <- unique(unlist(ped))
  }
  
  phenoQC <- function(filePhe, planPhe) {
    if (is.character(filePhe)) {
      pheno <- read.table(filePhe, sep = sep, header = header)
    } else {
      pheno <- filePhe
    }
    
    if (!is.null(filePed)) {
      pheno <- pheno[pheno[, 1] %in% ped, ]
    }
    
    if (is.null(planPhe)) {
      # auto select columns
      if(is.null(pheCols)) {
        pheCols <- c(1:ncol(pheno))
      }
      # check phenotype file
      if (length(pheCols) < 2) {
        stop("ERROR: At least 2 columns in the phenotype file should be specified, please check the parameter 'pheSep'. ")
      }
      hasRep <- FALSE
      pheList <- pheno[, pheCols]
      
    } else {
      logging.log(" JOB NAME:", planPhe$job_name, "\n", verbose = verbose)
      pheDef <- sapply(planPhe$job_traits, function(plan) {
        return(plan$definition)
      })
      pheName <- sapply(planPhe$job_traits, function(plan) {
        if (!is.null(plan$trait)) { return(plan$trait) }
        if (!is.null(plan$environment)) { return(plan$environment) }
      })
      envName <- sapply(planPhe$job_traits, function(plan) {
        return(plan$environment)
      })
      pheCond <- lapply(planPhe$job_traits, function(plan) {
        return(plan$range)
      })
      useFlag <- sapply(1:length(pheDef), function(i) {
        flag <- FALSE
        for (j in 1:ncol(pheno)) {
          flagt <- grepl(pattern = names(pheno)[j], pheDef[i])
          if (flagt) { flag = TRUE; break; }
        }
        return(flag)
      })
      usePheDef <- pheDef[useFlag]
      usePheName <- pheName[useFlag]
      useEnvName <- unlist(envName[useFlag])
      usePheCond <- pheCond[useFlag]
      
      # get new phenotype
      newPheDef <- setdiff(usePheDef, usePheName)
      if (length(newPheDef) == 0) {
        pheList <- pheno
      } else { 
        newPheName <- usePheName[match(newPheDef, usePheDef)]
        rmCol <- match(newPheName, names(pheno))
        rmCol <- rmCol[!is.na(rmCol)]
        if (length(rmCol) > 0) { pheno <- pheno[, -rmCol, drop = FALSE]  }
        newPheList <- do.call(cbind, lapply(1:length(newPheDef), function(i) {
          return(data.frame(with(pheno, eval(parse(text = newPheDef[i])))))
        }))
        names(newPheList) <- newPheName
        pheList <- cbind(pheno, newPheList)
      }
      
      # data filter & select & arrange
      if (length(planPhe$filter) > 0) {
        filterRow <- with(pheList, eval(parse(text = planPhe$filter[1])))
        filterRow[is.na(filterRow)] <- FALSE
        pheList <- pheList[filterRow, ] 
      }
      if (length(planPhe$select) > 0) {
        pheList <- pheList[, planPhe$select] 
      }
      if (length(planPhe$arrange) > 0) {
        if (length(planPhe$decreasing) > 0) {
          decreasing <- planPhe$decreasing
        } else {
          decreasing <- FALSE
        }
        orderCmd <- paste0("order(", paste0("pheList$", planPhe$arrange, collapse = ","), ",decreasing = decreasing)")
        pheList <- pheList[eval(parse(text = orderCmd)), ]
      }

      # remove abnormal values
      noRange <- FALSE
      for (i in 1:length(usePheName)) {
        pheList[[usePheName[i]]][is.infinite(pheList[[usePheName[i]]])] <- NA
        pheList[[usePheName[i]]][is.nan(pheList[[usePheName[i]]])] <- NA
        if (length(usePheCond[[i]]) == 0) {
          if (is.numeric(pheList[[usePheName[i]]])) {
            noRange <- TRUE
            mn <- mean(pheList[[usePheName[i]]], na.rm = TRUE)
            sd <- sd(pheList[[usePheName[i]]], na.rm = TRUE)
            f1 <- pheList[[usePheName[i]]] >= (mn - 3 * sd)
            f2 <- pheList[[usePheName[i]]] <= (mn + 3 * sd)
            f1[is.na(f1)] <- TRUE
            f2[is.na(f2)] <- TRUE
            pheList[[usePheName[i]]][!(f1 & f2)] <- NA
          }
          
        } else if (is.numeric(usePheCond[[i]])) {
          if (length(usePheCond[[i]]) != 2) {
            stop("Numeric parameter 'range' should only contain the minimum and maximum!")
          }
          
          f1 <- pheList[[usePheName[i]]] >= usePheCond[[i]][1]
          f2 <- pheList[[usePheName[i]]] <= usePheCond[[i]][2]
          f1[is.na(f1)] <- TRUE
          f2[is.na(f2)] <- TRUE
          pheList[[usePheName[i]]][!(f1 & f2)] <- NA
          
        } else if (!is.numeric(usePheCond[[i]])) {
          f1 <- pheList[[usePheName[i]]] %in% usePheCond[[i]]
          pheList[[usePheName[i]]][!f1] <- NA 
          
        } else {
          stop("This 'range' may be wrong!")
        }
      }
      
      if (noRange) {
        logging.log(" The traits without classification have been cleaned by [mean-3*sd, mean+3*sd]!\n", verbose = verbose)
      }
      
      # check non-repeat record trait
      hasRep <- planPhe$repeated_records
      if (!hasRep) {
        pheList <- pheList[!duplicated(pheList[, 1]), ]
      }
      
    } # if (is.null(planPhe)) {
    
    # remove spaces of elements
    finalPhe <- data.frame(lapply(pheList, function(x){ gsub("\\s+", "", x) }))
    # remove column ofs full of NAs
    drop <- c()
    for (i in 2:ncol(finalPhe)) {
      finalPhe[finalPhe[, i] %in% missing, i] <- NA
      if (all(is.na(finalPhe[, i]))) {
        drop <- c(drop, i)
      }
    }
    if (length(drop) > 0) {
      finalPhe <- finalPhe[, -drop]
    }

    # rename header
    if (!header)  {
      colnames(finalPhe)[1] <- 'Taxa'
      traits <- 2:ncol(finalPhe)
      colnames(finalPhe)[traits] <- paste0('t', traits - 1)
    }
    
    # Output
    if (is.null(out)) {
      out <- unlist(strsplit(filePhe, split = '.', fixed = TRUE))[1]
    }
    if (hasRep) {
      pheFileName <- paste0(out, ".repeat.qc.txt")
    } else {
      pheFileName <- paste0(out, ".qc.txt")
    }
    write.table(finalPhe, pheFileName, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    return(pheFileName)
  } # end function phenoQC
  
  pheFileName <- NULL
  for (i in 1:length(filePhe)) {
    pheFileName <- c(pheFileName, phenoQC(filePhe[i], planPhe[[i]]))
  }
  
  t2 <- as.numeric(Sys.time())
  logging.log(" Preparation for PHENOTYPE data for is Done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(pheFileName)
}

#' To find appropriate fixed effects, covariates, and random effects
#' Author: Dong Yin
#' Build date: July 17, 2021
#' Last update: Apr 21, 2022
#' 
#' @param jsonList the list of json parameters
#' @param header the header of file
#' @param sep the separator of file
#' @param ncpus the number of threads
#' @param verbose whether to print detail.
#'
#' @return the best effects for EBV model
#' @export
#'
#' @examples
#' jsonFile <- system.file("extdata", "demo2.json", package = "simer")
#' jsonList <- rjson::fromJSON(file = jsonFile)
#' # jsonList <- simer.Data.Env(jsonList = jsonList)
simer.Data.Env <- function(jsonList = NULL, header = TRUE, sep = '\t', ncpus = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  
  genoPath <- jsonList$genotype
  genoFiles <- list.files(genoPath)
  fileMVP <- grep(pattern = "geno.desc", genoFiles, value = TRUE)
  fileMVP <- file.path(genoPath, fileMVP)
  fileMVP <- substr(fileMVP, 1, nchar(fileMVP)-10)
  filePed <- jsonList$pedigree
  planPhe <- jsonList$analysis_plan
  
  for (i in 1:length(planPhe)) {
    logging.log(" JOB NAME:", planPhe[[i]]$job_name, "\n", verbose = verbose)
    filePhe <- planPhe[[i]]$sample_info
    pheno <- read.table(filePhe, header = header, sep = sep)
    multrait <- planPhe[[i]]$multi_trait
    randomRatio <- planPhe[[i]]$random_ratio
    planPhe[[i]]$random_ratio <- NULL
    nTrait <- length(planPhe[[i]]$job_traits)
    jsonListN <- jsonList
    if (multrait) {
      planPheN <- lapply(1:nTrait, function(j) { return(planPhe[[i]]) })
      for (j in 1:nTrait) {
        planPheN[[j]]$multi_trait <- FALSE
        planPheN[[j]]$job_traits <- planPheN[[j]]$job_traits[j]
      }
    } else {
      planPheN <- planPhe[i]
    }
    for (j in 1:nTrait) {
      traits <- planPheN[[j]]$job_traits[[1]]$traits
      covariates <- unlist(planPheN[[j]]$job_traits[[1]]$covariates)
      fixedEffects <- unlist(planPheN[[j]]$job_traits[[1]]$fixed_effects)
      randomEffects <- unlist(planPheN[[j]]$job_traits[[1]]$random_effects)
      covariates <- covariates[covariates %in% names(pheno)]
      fixedEffects <- fixedEffects[fixedEffects %in% names(pheno)]
      randomEffects <- randomEffects[randomEffects %in% names(pheno)]
      finalPhe <- pheno[, c(names(pheno)[1], traits, covariates, fixedEffects, randomEffects)]
      # remove all NAs
      finalPhe <- na.omit(finalPhe)
      lapply(covariates, function(col) {  mode(finalPhe[, col]) <<- "numeric" })
      lapply(fixedEffects, function(col) {  mode(finalPhe[, col]) <<- "character" })
      # remove column of one level or full levels
      finalPhe <- checkEnv(finalPhe, c(covariates, fixedEffects, randomEffects))
      covariates <- covariates[covariates %in% names(finalPhe)]
      fixedEffects <- fixedEffects[fixedEffects %in% names(finalPhe)]
      randomEffects <- randomEffects[randomEffects %in% names(finalPhe)]
      lmPhe <- lm(paste(paste0(traits, "~1"), 
                        paste(unlist(c(covariates, fixedEffects)), collapse = "+"), 
                        sep = "+"), data = finalPhe)
      # choose a model by BIC in a stepwise algorithm
      file <- NULL
      if (verbose) {
        try(file <- get("logging.file", envir = package.env), silent = TRUE)
        sink(file = file, append = TRUE, split = TRUE)
      }
      slmPhe <- step(lmPhe, k = log(nrow(finalPhe)))
      if (verbose) {
        sink()
      }
      
      envName <- names(slmPhe$model)[-1]
      covariates <- covariates[covariates %in% envName]
      fixedEffects <- fixedEffects[fixedEffects %in% envName]
      # reset covariates, fixed effects
      planPheN[[j]]$job_traits[[1]]$covariates <- covariates
      planPheN[[j]]$job_traits[[1]]$fixed_effects <- fixedEffects
      jsonListN$analysis_plan <- list(planPheN[[j]])
      # select random effect which ratio less than threshold
      gebv <- simer.Data.cHIBLUP(jsonList = jsonListN, ncpus = ncpus, verbose = verbose)
      vc <- gebv[[1]]$varList[[1]]
      randomEffectRatio <- vc[1:length(randomEffects)] / sum(vc)
      randomEffects <- randomEffects[randomEffectRatio > randomRatio]
      # reset covariates, fixed effects and random effects
      planPhe[[i]]$job_traits[[j]]$covariates <- covariates
      planPhe[[i]]$job_traits[[j]]$fixed_effects <- fixedEffects
      planPhe[[i]]$job_traits[[j]]$random_effects <- randomEffects
      envFormula <- c(paste0(covariates, "(C)"), paste0(fixedEffects, "(F)"), paste0(randomEffects, "(R)"))
      envFormula <- envFormula[nchar(envFormula) > 3]
      logging.log(" *********************************************************\n",
                  "Model optimized by BIC and random variance ratio is:\n", 
                    paste(c(paste0(traits, "~1"), envFormula), collapse = '+'), "\n",
                  "*********************************************************\n", verbose = verbose)
    }
  }
  
  jsonList$analysis_plan <- planPhe
  t2 <- as.numeric(Sys.time())
  logging.log(" Model optimization is Done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(jsonList)
}

#' simer.Data.cHIBLUP: The function of calling HIBLUP software of C version
#' Author: Dong Yin
#' Build date: June 28, 2021
#' Last update: Apr 21, 2022
#'
#' @param SP a list of all simulation parameters
#' @param jsonList the list of json parameters
#' @param mode 'A' or 'AD', Additive effect model or Additive and Dominance mode
#' @param vc.method default is 'AI', the method of calucating variance components in HIBLUP
#' @param ncpus the number of threads
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return a list of gebv
#' 
#' @examples
#' jsonFile <- system.file("extdata", "demo2.json", package = "simer")
#' jsonList <- rjson::fromJSON(file = jsonFile)
#' # gebvs <- simer.Data.cHIBLUP(SP = SP)
#' # gebvs <- simer.Data.cHIBLUP(jsonList = jsonList)
simer.Data.cHIBLUP <- function(SP = NULL, jsonList = NULL, mode='A', vc.method = "AI", ncpus = 10, verbose=TRUE) {
  t1 <- as.numeric(Sys.time())
  
  if (!is.null(SP)) {
    replication <- SP$global$replication
    out <- SP$global$out
    outpath <- SP$global$outpath
    out.format <- SP$global$out.format
    incols <- SP$geno$incols
    pop.inds <- sapply(1:SP$reprod$pop.gen, function(i) {
      return(ncol(SP$geno$pop.geno[[i]]) / incols)
    })
    pop.ind <- sum(pop.inds)
    outpath = paste0(outpath, .Platform$file.sep, pop.ind, "_Simer_Data_plink")
    directory.rep <- paste0(outpath, .Platform$file.sep, "replication", replication)
    fileMVP <- file.path(directory.rep, out)
    filePed <- paste0(fileMVP, ".ped")
    nTrait <- length(SP$pheno$model)
    job_traits <- lapply(1:nTrait, function(i) {
      
    })
    planPhe <- list(list(
      job_name = "EBV Model Demo",
      sample_info = paste0(fileMVP, ".phe"),
      repeated_records = ifelse(SP$pheno$pop.rep == 1, FALSE, TRUE),
      multi_trait = ifelse(nTrait == 1, FALSE, TRUE),
      random_ratio = 0.05,
      
    ))
    planPhe <- jsonList$analysis_plan
  }
  
  if (!is.null(jsonList)) {
    genoPath <- jsonList$genotype
    genoFiles <- list.files(genoPath)
    fileMVP <- grep(pattern = "geno.desc", genoFiles, value = TRUE)
    fileMVP <- file.path(genoPath, fileMVP)
    fileMVP <- substr(fileMVP, 1, nchar(fileMVP)-10)
    filePed <- jsonList$pedigree
    planPhe <- jsonList$analysis_plan
  }
  
  # convert bigmemory to PLINK
  if (!is.null(fileMVP)) {
    filePhe <- planPhe[[1]]$sample_info
    bigmat <- attach.big.matrix(paste0(fileMVP, ".geno.desc"))
    map <- read.table(paste0(fileMVP, ".geno.map"), header = TRUE)
    pheno <- read.table(filePhe, header = TRUE)
    traits <- sapply(planPhe[[1]]$job_traits, function(x) return(x$trait))
    MVP.Data.MVP2Bfile(
      bigmat = bigmat, 
      map = map, 
      pheno = pheno[, c(1, match(traits[1], names(pheno)))],
      out = fileMVP,
      verbose = verbose
    )
  }
  
  gebvs <- NULL
  for (i in 1:length(planPhe)) {
    logging.log(" JOB NAME:", planPhe[[i]]$job_name, "\n", verbose = verbose)
    filePhe <- planPhe[[i]]$sample_info
    pheno <- read.table(filePhe, header = TRUE)
    traits <- sapply(planPhe[[i]]$job_traits, function(x) return(x$trait))
    covariates <- unique(c(unlist(sapply(planPhe[[i]]$job_traits, function(x) {
      return(unlist(c(x$covariates)))
    }))))
    fixedEffects <- unique(c(unlist(sapply(planPhe[[i]]$job_traits, function(x) {
      return(unlist(c(x$fixed_effects)))
    }))))
    randomEffects <- unique(c(unlist(sapply(planPhe[[i]]$job_traits, function(x) {
      return(unlist(c(x$random_effects)))
    }))))
    covariates <- covariates[covariates %in% names(pheno)]
    fixedEffects <- fixedEffects[fixedEffects %in% names(pheno)]
    randomEffects <- randomEffects[randomEffects %in% names(pheno)]
    # get useful phenotype data
    finalPhe <- pheno[, c(names(pheno)[1], traits, covariates, fixedEffects, randomEffects)]
    # remove all NAs
    finalPhe <- na.omit(finalPhe)
    lapply(covariates, function(col) {  mode(finalPhe[, col]) <<- "numeric" })
    lapply(fixedEffects, function(col) {  mode(finalPhe[, col]) <<- "character" })
    # remove column of one level or full levels
    finalPhe <- checkEnv(finalPhe, c(covariates, fixedEffects, randomEffects))
    # show data information
    logging.log(" Data Summary:\n", verbose = verbose)
    logging.print(summary(finalPhe), verbose = verbose)
    
    covariatesCmd <- lapply(planPhe[[i]]$job_traits, function(x) {
        x$covariates <- x$covariates[x$covariates %in% names(finalPhe)]
        env <- paste0(match(x$covariates, names(pheno)), collapse=',')
        if (env == "") return(0)
        return(env)
    })
    fixedEffectsCmd <- lapply(planPhe[[i]]$job_traits, function(x) {
        x$fixed_effects <- x$fixed_effects[x$fixed_effects %in% names(finalPhe)]
        env <- paste0(match(x$fixed_effects, names(pheno)), collapse=',')
        if (env == "") return(0)
        return(env)
    })
    randomEffectsCmd <- lapply(planPhe[[i]]$job_traits, function(x) {
        x$random_effects <- x$random_effects[x$random_effects %in% names(finalPhe)]
        env <- paste0(match(x$random_effects, names(pheno)), collapse=',')
        if (env == "") return(0)
        return(env)
    })
    phenoCmd <- match(traits, names(pheno))
    if (!planPhe[[i]]$multi_trait) {
      covariatesCmd <- covariatesCmd[[1]]
      fixedEffectsCmd <- fixedEffectsCmd[[1]]
      randomEffectsCmd <- randomEffectsCmd[[1]]
      phenoCmd <- phenoCmd[1]
      traits <- traits[1]
    }
    covariatesCmd <- paste(covariatesCmd, collapse=' ')
    fixedEffectsCmd <- paste(fixedEffectsCmd, collapse=' ')
    randomEffectsCmd <- paste(randomEffectsCmd, collapse=' ')
    phenoCmd <- paste(phenoCmd, collapse=' ')
    nTraitCmd <- ifelse(planPhe[[i]]$multi_trait, "--multi-trait", "--single-trait")
    out <- paste(traits, collapse = '_')

    completeCmd <- 
      paste("hiblup", nTraitCmd,
        paste("--pheno", filePhe),
        paste("--pheno-pos", phenoCmd),
        paste("--dcovar", fixedEffectsCmd),
        paste("--qcovar", covariatesCmd),
        paste("--rand", randomEffectsCmd),
        ifelse(is.null(filePed), "", paste("--pedigree", filePed)),
        ifelse(is.null(fileMVP), "", paste("--bfile", fileMVP)),
        ifelse(mode == "A", "--add", paste("--add", "--dom")),
        paste("--vc-method", vc.method),
        paste("--threads", ncpus),
        paste("--out", out)
      )
    
    system(completeCmd)
    
    gebv <- NULL

    varFile <- paste0(out, ".vars")
    vars <- read.table(varFile, header = TRUE)
    varList <- lapply(traits, function(trait) {
      if (planPhe[[i]]$multi_trait) {
        return(vars[grep(pattern = trait, vars[, 1]), 2])
      } else {
        return(vars[, 2])
      } 
    })
    names(varList) <- traits
    gebv$varList <- varList
    logging.log(" The variance components are: (R, G, E)\n", verbose = verbose)
    logging.print(gebv$varList, verbose = verbose)

    if (planPhe[[i]]$multi_trait) {
      covarFile <- paste0(out, ".covars")
      covars <- read.table(covarFile, header = TRUE)
      covA <- corA <- matrix(0, length(traits), length(traits))
      for (j in 1:length(traits)) {
        for (k in 1:j) {
          if (j == k) {
            covA[j, k] <- varList[[j]][length(varList[[j]])-1]
            corA[j, k] <- 1
          } else {
            covA[j, k] <- covA[k, j] <- covars[j+k-2, 2]
            corA[j, k] <- corA[k, j] <- covars[j+k-2, 4]
          }
        }
      }
      dimnames(covA) <- dimnames(corA) <- list(paste0(" ", traits), paste0(" ", traits))
      gebv$covA <- covA
      gebv$corA <- corA
      logging.log(" The genetic covariance are:\n", verbose = verbose)
      logging.print(gebv$covA, verbose = verbose)
      logging.log(" The genetic correlation are:\n", verbose = verbose)
      logging.print(gebv$corA, verbose = verbose)
    }
    
    gebvs[[i]] <- gebv
  }
  
  t2 <- as.numeric(Sys.time())
  logging.log(" Genetic assessment for ", basename(filePhe), "is Done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(gebvs)
}

#' simer.SELIND: The function of General Selection Index
#' Author: Dong Yin
#' Build date: Aug 26, 2021
#' Last update: Apr 21, 2022
#'
#' @param jsonList the list of json parameters
#' @param ncpus the number of threads
#' @param verbose whether to print detail.
#'
#' @return phenotype economic weight
#' @export
#' @references Y. S. Chen, Z. L. Sheng (1988) The Theory of General Selection Index. Genetic Report, 15(3): P185-P190
#'
#' @examples
#' jsonFile <- system.file("extdata", "demo2.json", package = "simer")
#' jsonList <- rjson::fromJSON(file = jsonFile)
#' # jsonList <- simer.Data.SELIND(jsonList = jsonList)
simer.Data.SELIND <- function(jsonList = NULL, ncpus = 10, verbose=TRUE) {
  t1 <- as.numeric(Sys.time())
  
  BVIndex <- jsonList$breeding_value_index
  genoPath <- jsonList$genotype
  genoFiles <- list.files(genoPath)
  fileMVP <- grep(pattern = "geno.desc", genoFiles, value = TRUE)
  fileMVP <- file.path(genoPath, fileMVP)
  fileMVP <- substr(fileMVP, 1, nchar(fileMVP)-10)
  filePed <- jsonList$pedigree
  planPhe <- jsonList$analysis_plan
  
  str1 <- unlist(strsplit(BVIndex, split = c("\\+|\\*")))
  str1 <- gsub("^\\s+|\\s+$", "", str1)
  strIsNA <- is.na(suppressWarnings(as.numeric(str1)))
  BVWeight <- as.numeric(str1[!strIsNA])
  names(BVWeight) <- str1[strIsNA]
  
  gebvs <- simer.Data.cHIBLUP(jsonList = jsonList, ncpus = ncpus, verbose = verbose)
  
  covPList <- NULL
  covAList <- NULL
  pheNames <- NULL
  for (i in 1:length(planPhe)) {
    # prepare phenotype data
    filePhe <- planPhe[[i]]$sample_info
    pheno <- read.table(filePhe, header = TRUE)
    pheName <- sapply(planPhe[[i]]$job_traits, function(x) {
      return(x$trait)
    })
    pheNames <- c(pheNames, pheName)
    if (!all(pheName %in% names(BVWeight))) {
      stop(pheName[!(pheName %in% names(BVWeight))], " are not in the 'BVIndex'!")
    }
    if (planPhe[[i]]$repeated_records) {
      simer.mean <- function(x) { return(mean(x, na.rm = TRUE)) }
      usePhe <- sapply(pheName, function(name) {
        return(tapply(pheno[, name], as.factor(pheno[, 1]), FUN = simer.mean))
      })
      covP <- var(usePhe, na.rm = TRUE)
    } else {
      covP <- var(pheno[, pheName, drop = FALSE], na.rm = TRUE)
    }
    covPList[[i]] <- covP
    if (planPhe[[i]]$multi_trait) {
      covAList[[i]] <- gebvs[[i]]$covA
    } else {
      covAList[[i]] <- gebvs[[i]]$varList[[1]][1]
    }
  }
  
  P <- as.matrix(Matrix::bdiag(covPList))
  A <- as.matrix(Matrix::bdiag(covAList))
  iP <- try(solve(P), silent = TRUE)
  if (inherits(iP, "try-error")) {
    iP <- MASS::ginv(P)
  }

  # selection index
  if (any(sort(pheNames) != sort(names(BVWeight)))) {
    stop("Trait names should be consistent between planPhe and BVWeight!")
  }
  BVWeight <- BVWeight[match(pheNames, names(BVWeight))]
  b <- iP %*% A %*% BVWeight
  b <- round(as.vector(b), digits = 2)
  selIndex <- paste(paste(b, names(BVWeight), sep = "*"), collapse = " + ")

  logging.log(" *********************************************************\n",
                "General Selection Index is:\n", 
                selIndex, "\n",
                "*********************************************************\n", verbose = verbose)

  jsonList$selection_index <- selIndex
  jsonList$breeding_value_index <- NULL
  t2 <- as.numeric(Sys.time())
  logging.log(" Selection Index Construction is Done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(jsonList)
}

#' Check the levels of environmental factors
#'
#' Build date: Sep 10, 2021
#' Last update: Sep 10, 2021
#'
#' @author Dong Yin
#'
#' @param data data needing check
#' @param envName the environmental factor name in the data
#' 
#' @return data without environmental factors of wrong level
#' @export
#'
#' @examples
#' data <- data.frame(a = rep(1, 3), b = 1:3, c = c(1, 1, 5))
#' envName <- c("a", "b", "c")
#' # data <- checkEnv(data = data, envName = envName)
checkEnv <- function(data, envName) {
  if (is.numeric(envName)) {
    envName <- names(data)[envName]
  }
  
  # remove column of one level or full levels
  drop <- c()
  for (i in 1:ncol(data)) {
    if (names(data)[i] %in% envName) {
      numUni <- length(unique(data[[i]]))
      if (numUni == 1 || numUni == length(data[[i]])) {
        drop <- c(drop, i)
      }
    }
  }
  if (length(drop) > 0) {
    warning(paste(names(data)[drop], collapse=', '), ' has been remove because of its one or full levels!')
    data <- data[, -drop, drop = FALSE]
  }
  return(data)
}
