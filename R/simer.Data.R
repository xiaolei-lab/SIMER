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


#' Data handling
#' 
#' Make data quality control for genotype, phenotype, and pedigree.
#' 
#' Build date: May 26, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#'
#' @param jsonList a list of data quality control parameters.
#' @param out the prefix of output files.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$genotype}{the path of genotype data.}
#' \item{$pedigree}{the filename of pedigree data.}
#' \item{$selection_index}{the selection index for all traits.}
#' \item{$breeding_value_index}{the breeding value index for all traits.}
#' \item{$quality_control_plan}{a list of parameters for data quality control.}
#' \item{$analysis_plan}{a list of parameters for genetic evaluation.}
#' }
#' 
#' @examples
#' # Read JSON file
#' jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
#' jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)
#' 
#' \dontrun{
#' # It needs 'plink' and 'hiblup' software
#' jsonList <- simer.Data(jsonList = jsonList)
#' }
simer.Data <- function(jsonList = NULL, out = 'simer.qc', ncpus = 0, verbose = TRUE) {
  
  # global parameters
  outpath <- dirname(out)
  if (!dir.exists(outpath)) { dir.create(outpath) }
  maxLine <- 10000
  priority <- "speed"
  
  genoPath <- unlist(jsonList$genotype)
  fileMVP <- fileBed <- filePlinkPed <-  NULL
  if (length(genoPath) != 0) {
    genoFiles <- list.files(genoPath)
    fileMVP <- grep(pattern = "geno.desc", genoFiles, value = TRUE)
    fileMVP <- file.path(genoPath, fileMVP)
    if (length(fileMVP) > 0) {
      fileMVP <- substr(fileMVP, 1, nchar(fileMVP)-10)
    } else {
      fileMVP <- NULL
    }
    fileBed <- grep(pattern = "bed", genoFiles, value = TRUE)
    fileBed <- file.path(genoPath, fileBed)
    if (length(fileBed) > 0) {
      fileBed <- substr(fileBed, 1, nchar(fileBed)-4)
    } else {
      fileBed <- NULL
    }
    filePlinkPed <- grep(pattern = ".ped$", genoFiles, value = TRUE)
    if (length(filePlinkPed) > 0) {
      filePlinkPed <- substr(filePlinkPed, 1, nchar(filePlinkPed)-4)
      filePlinkPed <- file.path(genoPath, filePlinkPed)
    } else {
      filePlinkPed <- NULL
    }
  }
  genoType <- "char"
  filter_geno <- unlist(jsonList$quality_control_plan$genotype_quality_control$filter)
  filterGeno <- unlist(jsonList$quality_control_plan$genotype_quality_control$filter_geno)
  filterHWE <- unlist(jsonList$quality_control_plan$genotype_quality_control$filter_hwe)
  filterMind <- unlist(jsonList$quality_control_plan$genotype_quality_control$filter_mind)
  filterMAF <- unlist(jsonList$quality_control_plan$genotype_quality_control$filter_maf)
  
  filePed <- unlist(jsonList$pedigree)
  standardID <- unlist(jsonList$quality_control_plan$pedigree_quality_control$standard_ID)
  fileSir <- unlist(jsonList$quality_control_plan$pedigree_quality_control$candidate_sire_file)
  fileDam <- unlist(jsonList$quality_control_plan$pedigree_quality_control$candidate_dam_file)
  exclThres <- unlist(jsonList$quality_control_plan$pedigree_quality_control$exclude_threshold)
  assignThres <- unlist(jsonList$quality_control_plan$pedigree_quality_control$assign_threshold)
  pedSep <- "\t"
  
  filePhe <- unlist(sapply(jsonList$quality_control_plan$phenotype_quality_control, function(x) return(x$sample_info)))
  planPhe <- jsonList$quality_control_plan$phenotype_quality_control
  pheCols <- NULL
  pheSep <- "\t"
  missing = c(NA, 'NA', 'Na', '.', '-', 'NAN', 'nan', 'na', 'N/A', 'n/a', '<NA>', '', '-9', 9999)
  
  genoFileName <- pedFileName <- pheFileName <- NULL
  
  logging.initialize("Simer.Data", outpath)
  
  if (length(fileBed) != 0) {
    logging.log("*************** Genotype Data Quality Control ***************\n", verbose = verbose)
    genoFileName <-
      simer.Data.Geno(
        fileMVP = fileMVP,
        fileBed = fileBed, 
        filePlinkPed = filePlinkPed,
        filePed = filePed,
        filePhe = filePhe,
        out = out,
        genoType = genoType,
        filter = filter_geno,
        filterGeno = filterGeno,
        filterHWE = filterHWE,
        filterMind = filterMind,
        filterMAF = filterMAF,
        ncpus = ncpus,
        verbose = verbose)
    
    jsonList$genotype <- dirname(genoFileName)
  }
  
  if (length(filePed) != 0) {
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
    
    jsonList$pedigree <- pedFileName
  }
  
  if (length(filePhe) != 0) {
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
    
    for (i in 1:length(filePhe)) {
      jsonList$analysis_plan[[i]]$sample_info <- pheFileName[i]
    }
  }
  
  return(jsonList)
}

#' Genotype data conversion
#' 
#' Convert genotype data from MVP format to MVP format.
#' 
#' Build date: May 26, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param fileMVP the prefix of MVP file.
#' @param genoType type parameter in bigmemory data. The default is 'char', it is highly recommended *NOT* to modify this parameter.
#' @param out the prefix of output files.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return
#' the function returns files
#' \describe{
#' \item{<out>.geno.desc}{the description file of genotype data.}
#' \item{<out>.geno.bin}{the binary file of genotype data.}
#' \item{<out>.geno.ind}{the genotyped individual file.}
#' \item{<out>.geno.map}{the marker information data file.}
#' }
#' 
#' @examples
#' # Get the prefix of genotype data
#' fileMVP <- system.file("extdata", "01bigmemory", "demo", package = "simer")
#' 
#' # Convert genotype data from MVP to MVP
#' simer.Data.MVP2MVP(fileMVP, out = tempfile("outfile"))
simer.Data.MVP2MVP <- function(fileMVP, genoType = 'char', out = 'simer', verbose = TRUE) {
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

#' Genotype data imputation
#' 
#' Impute the missing value within genotype data.
#' 
#' Build date: May 26, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param fileMVP genotype in MVP format.
#' @param fileBed genotype in PLINK binary format.
#' @param out the name of output file.
#' @param maxLine number of SNPs, only used for saving memory when calculate kinship matrix.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @return 
#' the function returns files
#' \describe{
#' \item{<out>.geno.desc}{the description file of genotype data.}
#' \item{<out>.geno.bin}{the binary file of genotype data.}
#' \item{<out>.geno.ind}{the genotyped individual file.}
#' \item{<out>.geno.map}{the marker information data file.}
#' }
#' 
#' @export
#'
#' @examples
#' # Get the prefix of genotype data
#' fileMVP <- system.file("extdata", "02plinkb", "demo", package = "simer")
#' 
#' \dontrun{
#' # It needs 'beagle' software
#' fileMVPimp <- simer.Data.Impute(fileBed = fileBed)
#' }
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
      simer.Data.MVP2Bfile(bigmat = bigmat, map = map, out = tmpout, threads = ncpus, verbose = verbose)
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
  
  simer.Data.Bfile2MVP(bfile = tmpout, out = out, threads = ncpus, verbose = verbose)
  
  return(out)
}

#' Genotype data quality control
#' 
#' Data quality control for genotype data in MVP format and PLINK format.
#' 
#' Build date: May 26, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param fileMVP genotype in MVP format.
#' @param fileBed genotype in PLINK binary format.
#' @param filePlinkPed genotype in PLINK numeric format.
#' @param filePed the filename of pedigree data.
#' @param filePhe the filename of phenotype data, it can be a vector.
#' @param out the prefix of output files.
#' @param genoType type parameter in bigmemory, genotype data. The default is char, it is highly recommended *NOT* to modify this parameter.
#' @param filter filter of genotyped individual.
#' @param filterGeno threshold of sample miss rate.
#' @param filterHWE threshold of Hardy-Weinberg Test.
#' @param filterMind threshold of variant miss rate.
#' @param filterMAF threshold of Minor Allele Frequency.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' the function returns files
#' \describe{
#' \item{<out>.bed}{the .bed file of PLINK binary format.}
#' \item{<out>.bim}{the .bim file of PLINK binary format.}
#' \item{<out>.fam}{the .fam file of PLINK binary format.}
#' }
#' 
#' @examples
#' # Get the prefix of genotype data
#' fileBed <- system.file("extdata", "02plinkb", "demo", package = "simer")
#' 
#' \dontrun{
#' # It needs 'plink' software
#' simer.Data.Geno(fileBed=fileBed)
#' }
simer.Data.Geno <- function(fileMVP = NULL, fileBed = NULL, filePlinkPed = NULL, filePed = NULL, filePhe = NULL, out = 'simer.qc', genoType = 'char',
                            filter = NULL, filterGeno = NULL, filterHWE = NULL, filterMind = NULL, filterMAF = NULL,
                            ncpus = 0, verbose = TRUE) {
  
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Genotype Data.\n", verbose = verbose)
  
  if (length(filePed) != 0) {
    ped <-  read.table(filePed, sep = '\t', header = TRUE)
    keepInds <- unique(unlist(ped))
  } else {
    keepInds <- NULL
  }
  
  if (length(filePhe) != 0) {
    if (length(filter) > 0) {
      pheList <- read.table(filePhe[1], header = TRUE)
      filterRow <- pheList[with(pheList, eval(parse(text = filter))), 1]
      if (is.null(keepInds)) {
        keepInds <- filterRow
      } else {
        keepInds <- intersect(keepInds, filterRow)
      }
    }
  }
  
  if (FALSE) {
    if (length(fileMVP) == 0) { fileMVP <- NULL }
    if (length(filePed) == 0) { filePed <- NULL }
    if (is.null(out)) { out <- paste0(fileMVP, ".qc") }
    if (fileMVP != out) { remove_bigmatrix(out) }
    
    fileDesc <- normalizePath(paste0(fileMVP, '.geno.desc'), mustWork = TRUE)
    fileInd <- normalizePath(paste0(fileMVP, '.geno.ind'), mustWork = TRUE)
    fileMap <- normalizePath(paste0(fileMVP, '.geno.map'), mustWork = TRUE)
    
    genoInd <- read.table(fileInd, sep = '\t', header = FALSE)[, 1]
    genoMap <- read.table(fileMap, sep = '\t', header = TRUE)
    
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
  }
  
  if (!is.null(fileBed) | !is.null(filePlinkPed)) {
    if (!is.null(keepInds)) {
      keepInds <- cbind(keepInds, keepInds)
      write.table(keepInds, "simer.geno.ind", quote = FALSE, sep = ' ', row.names = FALSE, col.names = FALSE)
    }
    completeCmd <- 
      paste("plink", ifelse(is.null(fileBed), " --file", "--bfile"), fileBed,
            ifelse(length(keepInds) == 0, "", paste("--keep simer.geno.ind")),
            ifelse(length(filterGeno) == 0, "", paste("--geno", filterGeno)),
            ifelse(length(filterHWE) == 0, "", paste("--hwe", filterHWE)),
            ifelse(length(filterMind) == 0, "", paste("--mind", filterMind)),
            ifelse(length(filterMAF) == 0, "", paste("--maf", filterMAF)),
            "--make-bed --out", out)
    
    system(completeCmd)
  }
  
  t2 <- as.numeric(Sys.time())
  logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(out)
}

#' Pedigree data quality control
#' 
#' Data quality control for pedigree data.
#' 
#' Build date: May 6, 2021
#' Last update: Apr 28, 2022
#'
#' @author Lilin Yin and Dong Yin
#' 
#' @param filePed the filename of pedigree need correcting.
#' @param fileMVP genotype in MVP format.
#' @param out the prefix of output file.
#' @param standardID whether kid id is 15-character standard.
#' @param fileSir the filename of candidate sires.
#' @param fileDam the filename of candidate dams.
#' @param exclThres if conflict ratio is more than exclThres, exclude this parent.
#' @param assignThres if conflict ratio is less than assignThres, assign this parent to the individual.
#' @param header whether the file contains header.
#' @param sep separator of the file.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' the function returns files
#' \describe{
#' \item{<out>.report.ped}{the report file containing correction condition.}
#' \item{<out>.error.ped}{the file containing pedigree error.}
#' \item{<out>.ped}{the pedigree file after correction.}
#' }
#' 
#' @examples
#' # Get the filename of pedigree data
#' filePed <- system.file("extdata", "05others", "pedigree.txt", package = "simer")
#' 
#' # Run pedigree correction
#' simer.Data.Ped(filePed = filePed, out = tempfile("outfile"))
simer.Data.Ped <- function(filePed, fileMVP = NULL, out = NULL, standardID = FALSE, fileSir = NULL, fileDam = NULL, 
                           exclThres = 0.01, assignThres = 0.005, header = TRUE, sep = '\t', ncpus = 0, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Pedigree Data.\n", verbose = verbose)
  
  # if (!is.vector(filePed)) { filePed <- c(filePed) }
  if (length(filePed) == 0) { filePed <- NULL }
  
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
    
    # print("Making needed files")
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
  write.table(ped, paste0(out, ".report.ped"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  write.table(pedError, paste0(out, ".error.ped"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  write.table(pedUse, paste0(out, ".ped"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  
  t2 <- as.numeric(Sys.time())
  logging.log(" Preparation for PEDIGREE data is done within", format_time(t2 - t1), "\n\n", verbose = verbose)
  return(paste0(out, ".ped"))
}

#' Phenotype data quality control
#' 
#' Data quality control for phenotype data.
#' 
#' Build date: June 13, 2021
#' Last update: Apr 28, 2022
#'
#' @author Haohao Zhang and Dong Yin
#' 
#' @param filePhe the phenotype files, it can be a vector.
#' @param filePed the pedigree files, it can be a vector.
#' @param out the prefix of output file.
#' @param planPhe the plans for phenotype quality control.
#' @param pheCols the column needing extracting.
#' @param header the header of file.
#' @param sep the separator of file.
#' @param missing the missing value.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return 
#' the function returns files
#' \describe{
#' \item{<out>.phe}{the phenotype file after correction.}
#' }
#' 
#' @examples
#' # Get the filename of phenotype data
#' filePhe <- system.file("extdata", "05others", "phenotype.txt", package = "simer")
#' 
#' # Run phenotype correction
#' simer.Data.Pheno(filePhe = filePhe, out = tempfile("outfile"))
simer.Data.Pheno <- function(filePhe = NULL, filePed = NULL, out = NULL, planPhe = NULL, pheCols = NULL, header = TRUE, sep = '\t', missing = c(NA, 'NA', 'Na', '.', '-', 'NAN', 'nan', 'na', 'N/A', 'n/a', '<NA>', '', '-9', 9999), verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  logging.log(" Start Checking Phenotype Data.\n", verbose = verbose)
  
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
      # logging.log(" JOB NAME:", planPhe$job_name, "\n", verbose = verbose)
      pheDef <- sapply(planPhe$job_traits, function(plan) {
        return(unlist(plan$definition))
      })
      pheName <- sapply(planPhe$job_traits, function(plan) {
        if (!is.null(unlist(plan$trait))) { return(unlist(plan$trait)) }
        if (!is.null(unlist(plan$environment))) { return(unlist(plan$environment)) }
      })
      envName <- sapply(planPhe$job_traits, function(plan) {
        return(unlist(plan$environment))
      })
      pheCond <- lapply(planPhe$job_traits, function(plan) {
        return(unlist(plan$range))
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
      
      # filter & select & arrange
      if (length(unlist(planPhe$filter)) > 0) {
        filterRow <- with(pheList, eval(parse(text = unlist(planPhe$filter))))
        filterRow[is.na(filterRow)] <- FALSE
        pheList <- pheList[filterRow, ] 
      }
      if (length(unlist(planPhe$select)) > 0) {
        pheList <- pheList[, unlist(planPhe$select)] 
      }
      if (length(unlist(planPhe$arrange)) > 0) {
        if (length(unlist(planPhe$decreasing)) > 0) {
          decreasing <- unlist(planPhe$decreasing)
        } else {
          decreasing <- FALSE
        }
        orderCmd <- paste0("order(", paste0("pheList$", unlist(planPhe$arrange), collapse = ","), ",decreasing = decreasing)")
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
      hasRep <- unlist(planPhe$repeated_records)
      if (!hasRep) {
        pheList <- pheList[!duplicated(pheList[, 1]), ]
      }
      
    }
    
    # remove spaces in the phenotype data
    finalPhe <- data.frame(lapply(pheList, function(x){ gsub("\\s+", "", x) }))
    # remove column with full of NAs
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
    
    if (is.null(out)) {
      out <- unlist(strsplit(filePhe, split = '.', fixed = TRUE))[1]
    }
    if (hasRep) {
      pheFileName <- paste0(out, ".repeat.phe")
    } else {
      pheFileName <- paste0(out, ".phe")
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

#' Environmental factor selection
#' 
#' To find appropriate fixed effects, covariates, and random effects.
#' 
#' Build date: July 17, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param jsonList the list of environmental factor selection parameters.
#' @param hiblupPath the path of HIBLUP software.
#' @param header the header of file.
#' @param sep the separator of file.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$genotype}{the path of genotype data.}
#' \item{$pedigree}{the filename of pedigree data.}
#' \item{$selection_index}{the selection index for all traits.}
#' \item{$breeding_value_index}{the breeding value index for all traits.}
#' \item{$quality_control_plan}{a list of parameters for data quality control.}
#' \item{$analysis_plan}{a list of parameters for genetic evaluation.}
#' }
#' 
#' @export
#'
#' @examples
#' # Read JSON file
#' jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
#' jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)
#' 
#' \dontrun{
#' # It needs 'hiblup' solfware
#' jsonList <- simer.Data.Env(jsonList = jsonList)
#' }
simer.Data.Env <- function(jsonList = NULL, hiblupPath = '', header = TRUE, sep = '\t', ncpus = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  
  planPhe <- jsonList$analysis_plan
  auto_optim <- unlist(jsonList$auto_optimization)
  
  for (i in 1:length(planPhe)) {
    # logging.log(" JOB NAME:", planPhe[[i]]$job_name, "\n", verbose = verbose)
    filePhe <- unlist(planPhe[[i]]$sample_info)
    pheno <- read.table(filePhe, header = header, sep = sep)
    multrait <- unlist(planPhe[[i]]$multi_trait)
    randomRatio <- unlist(planPhe[[i]]$random_ratio)
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
      traits <- unlist(planPheN[[j]]$job_traits[[1]]$traits)
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
      finalPhe <- checkEnv(finalPhe, c(covariates, fixedEffects, randomEffects), verbose = verbose)
      covariates <- covariates[covariates %in% names(finalPhe)]
      fixedEffects <- fixedEffects[fixedEffects %in% names(finalPhe)]
      randomEffects <- randomEffects[randomEffects %in% names(finalPhe)]
      
      if (auto_optim) {
        if (length(c(covariates, fixedEffects)) != 0) {
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
          # reset covariates and fixed effects
          planPheN[[j]]$job_traits[[1]]$covariates <- covariates
          planPheN[[j]]$job_traits[[1]]$fixed_effects <- fixedEffects
        }
        
        jsonListN$analysis_plan <- list(planPheN[[j]])
        # select random effect which ratio less than threshold
        gebv <- simer.Data.cHIBLUP(jsonList = jsonListN, hiblupPath = hiblupPath, ncpus = ncpus, verbose = verbose)
        vc <- gebv[[1]]$varList[[1]]
        randomEffectRatio <- vc[1:length(randomEffects)] / sum(vc)
        randomEffects <- randomEffects[randomEffectRatio > randomRatio]
        # reset covariates, fixed effects, and random effects
        planPhe[[i]]$job_traits[[j]]$covariates <- covariates
        planPhe[[i]]$job_traits[[j]]$fixed_effects <- fixedEffects
        planPhe[[i]]$job_traits[[j]]$random_effects <- randomEffects
      }
      
      envFormula <- c(
        paste0(planPhe[[i]]$job_traits[[j]]$fixed_effects, "(F)"), 
        paste0(planPhe[[i]]$job_traits[[j]]$covariates, "(C)"), 
        paste0(planPhe[[i]]$job_traits[[j]]$random_effects, "(R)"))
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

#' Genetic evaluation
#' 
#' The function of calling HIBLUP software of C version.
#' 
#' Build date: June 28, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param jsonList the list of genetic evaluation parameters.
#' @param hiblupPath the path of HIBLUP software.
#' @param mode 'A' or 'AD', Additive effect model or Additive and Dominance model.
#' @param vc.method default is 'AI', the method of calculating variance components in HIBLUP software.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#' 
#' @export
#' 
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$randList}{a list of estimated random effects.}
#' \item{$varList}{a list of variance components.}
#' \item{$covA}{the genetic covariance matrix for all traits.}
#' \item{$corA}{the genetic correlation matrix for all traits.}
#' }
#' 
#' @examples
#' # Read JSON file
#' jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
#' jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)
#' 
#' \dontrun{
#' # It needs 'hiblup' software
#' gebvs <- simer.Data.cHIBLUP(jsonList = jsonList)
#' }
simer.Data.cHIBLUP <- function(jsonList = NULL, hiblupPath = '', mode = "A", vc.method = "AI", ncpus = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  
  genoPath <- unlist(jsonList$genotype)
  if (is.null(genoPath)) {  genoPath <- "" }
  genoFiles <- list.files(genoPath)
  fileBed <- grep(pattern = "bed", genoFiles, value = TRUE)
  fileBed <- file.path(genoPath, fileBed)
  fileBed <- substr(fileBed, 1, nchar(fileBed) - 4)
  filePed <- unlist(jsonList$pedigree)
  planPhe <- jsonList$analysis_plan
  
  gebvs <- NULL
  for (i in 1:length(planPhe)) {
    # logging.log(" JOB NAME:", planPhe[[i]]$job_name, "\n", verbose = verbose)
    filePhe <- unlist(planPhe[[i]]$sample_info)
    pheno <- read.table(filePhe, header = TRUE)
    traits <- sapply(planPhe[[i]]$job_traits, function(x) return(unlist(x$trait)))
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
    finalPhe <- checkEnv(finalPhe, c(covariates, fixedEffects, randomEffects), verbose = verbose)
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
        if (unlist(planPhe[[i]]$repeated_records)) {
          env <- paste0(c(1, env), collapse = ',')
        }
        if (env == "") return(0)
        return(env)
    })
    phenoCmd <- match(traits, names(pheno))
    if (!unlist(planPhe[[i]]$multi_trait)) {
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
    nTraitCmd <- ifelse(unlist(planPhe[[i]]$multi_trait), "--multi-trait", "--single-trait")
    out <- paste(traits, collapse = '_')

    completeCmd <- 
      paste(paste0(hiblupPath, "hiblup"), nTraitCmd,
        paste("--pheno", filePhe),
        paste("--pheno-pos", phenoCmd),
        paste("--dcovar", fixedEffectsCmd),
        paste("--qcovar", covariatesCmd),
        paste("--rand", randomEffectsCmd),
        ifelse(length(filePed) == 0, "", paste("--pedigree", filePed)),
        ifelse(length(fileBed) == 0, "", paste("--bfile", fileBed)),
        ifelse(mode == "A", "--add", paste("--add", "--dom")),
        paste("--vc-method", vc.method),
        paste("--threads", ncpus),
        paste("--out", out)
      )
    
    system(completeCmd)
    
    gebv <- NULL
    
    randList <- lapply(traits, function(trait) {
      if (unlist(planPhe[[i]]$multi_trait)) {
        randFile <- paste0(out, ".", trait, ".rand")
      } else {
        randFile <- paste0(trait, ".rand")
      }
      rand <- read.table(randFile, header = TRUE)
      return(rand[, c(1, ncol(rand) - 1)])
    })
    names(randList) <- traits
    gebv$randList <- randList
    
    varFile <- paste0(out, ".vars")
    vars <- read.table(varFile, header = TRUE)
    varList <- lapply(traits, function(trait) {
      if (unlist(planPhe[[i]]$multi_trait)) {
        return(vars[grep(pattern = trait, vars[, 1]), 2])
      } else {
        return(vars[, 2])
      } 
    })
    names(varList) <- traits
    gebv$varList <- varList
    logging.log(" The variance components are: (R, G, E)\n", verbose = verbose)
    logging.print(gebv$varList, verbose = verbose)

    if (unlist(planPhe[[i]]$multi_trait)) {
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

#' Selection index construction
#' 
#' The function of General Selection Index.
#' 
#' Build date: Aug 26, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param jsonList the list of selection index construction parameters.
#' @param hiblupPath the path of HIBLUP software.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return 
#' the function returns a list containing
#' \describe{
#' \item{$genotype}{the path of genotype data.}
#' \item{$pedigree}{the filename of pedigree data.}
#' \item{$selection_index}{the selection index for all traits.}
#' \item{$breeding_value_index}{the breeding value index for all traits.}
#' \item{$quality_control_plan}{a list of parameters for data quality control.}
#' \item{$analysis_plan}{a list of parameters for genetic evaluation.}
#' }
#' 
#' @export
#' 
#' @references Y. S. Chen, Z. L. Sheng (1988) The Theory of General Selection Index. Genetic Report, 15(3): P185-P190
#'
#' @examples
#' # Read JSON file
#' jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
#' jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)
#' 
#' \dontrun{
#' # It needs 'hiblup' software
#' jsonList <- simer.Data.SELIND(jsonList = jsonList)
#' }
simer.Data.SELIND <- function(jsonList = NULL, hiblupPath = '', ncpus = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  
  BVIndex <- unlist(jsonList$breeding_value_index)
  planPhe <- jsonList$analysis_plan
  auto_optim <- unlist(jsonList$auto_optimization)
  
  str1 <- unlist(strsplit(BVIndex, split = c("\\+|\\*")))
  str1 <- gsub("^\\s+|\\s+$", "", str1)
  strIsNA <- is.na(suppressWarnings(as.numeric(str1)))
  BVWeight <- as.numeric(str1[!strIsNA])
  names(BVWeight) <- str1[strIsNA]
  
  covPList <- NULL
  covAList <- NULL
  pheNames <- NULL
  usePhes <- NULL
  for (i in 1:length(planPhe)) {
    filePhe <- unlist(planPhe[[i]]$sample_info)
    pheno <- read.table(filePhe, header = TRUE)
    if (is.null(pheno$gen)) {
      pheno$gen <- 1
    }
    pheno <- pheno[pheno$gen == max(pheno$gen), ]
    pheName <- sapply(planPhe[[i]]$job_traits, function(x) {
      return(unlist(x$traits))
    })
    pheNames <- c(pheNames, pheName)
    if (!all(pheName %in% names(BVWeight))) {
      stop(pheName[!(pheName %in% names(BVWeight))], " are not in the 'BVIndex'!")
    }
    if (unlist(planPhe[[i]]$repeated_records)) {
      simer.mean <- function(x) { return(mean(x, na.rm = TRUE)) }
      usePhe <- sapply(pheName, function(name) {
        return(tapply(pheno[, name], as.factor(pheno[, 1]), FUN = simer.mean))
      })
      usePhe <- data.frame(rownames(usePhe), usePhe)
      names(usePhe)[1] <- names(pheno)[1]
    } else {
      usePhe <- pheno[, c(names(pheno)[1], pheName), drop = FALSE]
    }
    if (is.null(usePhes)) {
      usePhes <- usePhe
    } else {
      usePhes <- merge(x = usePhes, y = usePhe, by = names(pheno)[1], all = TRUE)
    }
    usePhe <- usePhe[, -1, drop = FALSE]
    covP <- var(usePhe, na.rm = TRUE)
    covPList[[i]] <- covP
  }
  
  usePhes <- usePhes[, -1]
  usePhes[is.na(usePhes)] <- 0
  
  if (auto_optim) {
    gebvs <- simer.Data.cHIBLUP(jsonList = jsonList, hiblupPath = hiblupPath, ncpus = ncpus, verbose = verbose)
    
    for (i in 1:length(planPhe)) {
      if (unlist(planPhe[[i]]$multi_trait)) {
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
  
  } else {
    if (is.null(unlist(jsonList$selection_index))) {
      stop("A selection index is necessary!")
    }
    str1 <- unlist(strsplit(unlist(jsonList$selection_index), split = c("\\+|\\*")))
    str1 <- gsub("^\\s+|\\s+$", "", str1)
    strIsNA <- is.na(suppressWarnings(as.numeric(str1)))
    b <- as.numeric(str1[!strIsNA])
    names(b) <- str1[strIsNA]
  } 
  
  selIndex <- paste(paste(b, names(BVWeight), sep = "*"), collapse = " + ")
  
  # genetic progress
  scores <- sort(as.matrix(usePhes) %*% b, decreasing = TRUE)
  geneticProgress <- round(mean(scores[1:(0.1*length(scores))]) - mean(scores), digits = 2)
  selIndex <- paste0("100 + ", selIndex)
  
  logging.log(" *********************************************************\n",
                "General Selection Index ~ Genetic Progress is:\n", 
                paste0(selIndex, " ~ ", geneticProgress), "\n",
                "*********************************************************\n", verbose = verbose)

  jsonList$selection_index <- selIndex
  jsonList$genetic_progress <- geneticProgress
  jsonList$breeding_value_index <- NULL
  t2 <- as.numeric(Sys.time())
  logging.log(" Selection Index Construction is Done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(jsonList)
}

#' Data quality control
#'
#' Make data quality control by JSON file.
#' 
#' Build date: Oct 19, 2020
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param jsonFile the path of JSON file.
#' @param hiblupPath the path of HIBLUP software.
#' @param out the prefix of output files.
#' @param dataQC whether to make data quality control.
#' @param buildModel whether to build EBV model.
#' @param buildIndex whether to build Selection Index.
#' @param ncpus the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print detail.
#'
#' @return 
#' the function returns a list containing
#' \describe{
#' \item{$genotype}{the path of genotype data.}
#' \item{$pedigree}{the filename of pedigree data.}
#' \item{$selection_index}{the selection index for all traits.}
#' \item{$breeding_value_index}{the breeding value index for all traits.}
#' \item{$quality_control_plan}{a list of parameters for data quality control.}
#' \item{$analysis_plan}{a list of parameters for genetic evaluation.}
#' }
#' 
#' @export
#' 
#' @examples
#' # Get JSON file
#' jsonFile <- system.file("extdata", "04breeding_plan", "plan1.json", package = "simer")
#' 
#' \dontrun{
#' # It needs 'plink' and 'hiblup' software
#' jsonList <- simer.Data.Json(jsonFile = jsonFile)
#' }
simer.Data.Json <- function(jsonFile, hiblupPath = '', out = "simer.qc", dataQC = TRUE, buildModel = TRUE, buildIndex = TRUE, ncpus = 10, verbose = TRUE) {
  
  jsonList <- jsonlite::fromJSON(txt = jsonFile, simplifyVector = FALSE)
  if (length(jsonList$threads) != 0) { ncpus <- jsonList$threads  }
  
  ## step 1. data quality control
  if (dataQC) {
    jsonList <- simer.Data(jsonList = jsonList, out = out, ncpus = ncpus, verbose = verbose)
  }
  
  ## step 2. find the best environmental factors for EBV model
  if (buildModel) {
    jsonList <- simer.Data.Env(jsonList = jsonList, hiblupPath = hiblupPath, ncpus = ncpus, verbose = verbose)
  }
  newJsonFile <- paste0(out, ".model.json")
  newJson <- jsonlite::toJSON(jsonList, pretty = TRUE, auto_unbox = TRUE)
  
  ## step 3. construct selection index
  if (buildIndex) {
    jsonList <- simer.Data.SELIND(jsonList = jsonList, hiblupPath = hiblupPath, ncpus = ncpus, verbose = verbose)
  }
  newJson <- jsonlite::toJSON(jsonList, pretty = TRUE, auto_unbox = TRUE)
  
  if (verbose) {
    cat(newJson, file = newJsonFile)
  }
  
  return(jsonList)
}

#' Environmental factor checking
#'
#' Check the levels of environmental factors.
#' 
#' Build date: Sep 10, 2021
#' Last update: Apr 28, 2022
#'
#' @author Dong Yin
#' 
#' @param data data needing check.
#' @param envName the environmental factor name within the data.
#' @param verbose whether to print detail.
#' 
#' @return data without environmental factors of wrong level.
#' 
#' @export
#'
#' @examples
#' data <- data.frame(a = c(1, 1, 2), b = c(2, 2, 3), c = c(3, 3, 4))
#' envName <- c("a", "b", "c")
#' data <- checkEnv(data = data, envName = envName)
checkEnv <- function(data, envName, verbose = TRUE) {
  if (is.numeric(envName)) {
    envName <- names(data)[envName]
  }
  
  # remove column(s) of one level or full levels
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
    logging.log(paste(names(data)[drop], collapse=', '), ' has been remove because of its one or full levels!\n', verbose = verbose)
    data <- data[, -drop, drop = FALSE]
  }
  return(data)
}

#' simer.Data.MVP2Bfile: To transform MVP data to binary format
#' 
#' transforming MVP data to binary format.
#' 
#' Build date: Sep 12, 2018
#' Last update: July 20, 2022
#'
#' @author Haohao Zhang and Dong Yin
#' 
#' @param bigmat Genotype in bigmatrix format (0,1,2).
#' @param map the map file.
#' @param pheno the phenotype file.
#' @param out the name of output file.
#' @param threads the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' @param verbose whether to print the reminder.
#'
#' @return NULL
#' Output files:
#' .bed, .bim, .fam
#' 
#' @export
#' 
#' @examples
#' # Generate bigmat and map
#' bigmat <- as.big.matrix(matrix(1:6, 3, 2))
#' map <- generate.map(pop.marker = 3)
#' 
#' # Data converting
#' simer.Data.MVP2Bfile(bigmat, map, out=tempfile("outfile"))
simer.Data.MVP2Bfile <- function(bigmat, map, pheno = NULL, out = 'simer', threads = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  
  logging.log(paste0("inds: ", ncol(bigmat), "\tmarkers:", nrow(bigmat), '\n'), verbose = verbose)
  
  # write bed file
  write_bfile(bigmat@address, out, threads = threads, verbose = verbose)
  
  # write fam
  #  1. Family ID ('FID')
  #  2. Within-family ID ('IID'; cannot be '0')
  #  3. Within-family ID of father ('0' if father isn't in dataset)
  #  4. Within-family ID of mother ('0' if mother isn't in dataset)
  #  5. Sex code ('1' = male, '2' = female, '0' = unknown)
  #  6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
  
  if (is.null(pheno)) {
    ind <- paste0("ind", 1:ncol(bigmat))
    sir <- rep(0, ncol(bigmat))
    dam <- rep(0, ncol(bigmat))
    sex <- rep(0, ncol(bigmat))
    pheno <- rep(-9, ncol(bigmat))
    message("pheno is NULL, automatically named individuals.")
    
  } else if (ncol(pheno) == 1) {
    ind <- pheno[, 1]
    sir <- rep(0, ncol(bigmat))
    dam <- rep(0, ncol(bigmat))
    sex <- rep(0, ncol(bigmat))
    pheno <- rep(-9, ncol(bigmat))
    
  } else {
    if (ncol(pheno) > 2) { 
      message("Only the first phenotype is written to the fam file, and the remaining phenotypes are ignored.")
    }
    ind <- pheno[, 1]
    sir <- pheno[, 5]
    dam <- pheno[, 6]
    sex <- pheno[, 7]
    pheno <- pheno[, 8]
  }
  
  fam <- cbind(ind, ind, sir, dam, sex, pheno)
  write.table(fam, paste0(out, '.fam'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
  
  # write bim
  #  1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  #  2. Variant identifier
  #  3. Position in morgans or centimorgans (safe to use dummy value of '0')
  #  4. Base-pair coordinate (normally 1-based, but 0 ok; limited to 231-2)
  #  5. Allele 1 (corresponding to clear bits in .bed; usually minor)
  #  6. Allele 2 (corresponding to set bits in .bed; usually major)
  bim <- cbind(map[, 2], map[, 1], 0, map[, 3], map[, 4], map[, 5])
  write.table(bim, paste0(out, '.bim'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
  t2 <- as.numeric(Sys.time())
  logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
}

#' simer.Data.Bfile2MVP: To transform plink binary data to MVP package
#' 
#' transforming plink binary data to MVP package.
#' 
#' Build date: Sep 12, 2018
#' Last update: July 25, 2022
#'
#' @author Haohao Zhang and Dong Yin
#' 
#' @param bfile Genotype in binary format (.bed, .bim, .fam).
#' @param out the name of output file.
#' @param maxLine the max number of line to write to big matrix for each loop.
#' @param priority 'memory' or 'speed'.
#' @param type.geno the type of genotype elements.
#' @param threads number of thread for transforming.
#' @param verbose whether to print the reminder.
#'
#' @return number of individuals and markers.
#' Output files:
#' genotype.desc, genotype.bin: genotype file in bigmemory format
#' phenotype.phe: ordered phenotype file, same taxa order with genotype file
#' map.map: SNP information
#' 
#' @export
#' 
#' @examples
#' # Get bfile path
#' bfilePath <- file.path(system.file("extdata", "02plinkb", package = "simer"), "demo")
#' 
#' # Data converting
#' simer.Data.Bfile2MVP(bfilePath, tempfile("outfile"))
simer.Data.Bfile2MVP <- function(bfile, out = 'simer', maxLine = 1e4, priority = 'speed', type.geno = 'char', threads = 10, verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  bim_file <- normalizePath(paste0(bfile, '.bim'), mustWork = TRUE)
  fam_file <- normalizePath(paste0(bfile, '.fam'), mustWork = TRUE)
  bed_file <- normalizePath(paste0(bfile, '.bed'), mustWork = TRUE)
  # check old file
  backingfile <- paste0(basename(out), ".geno.bin")
  descriptorfile <- paste0(basename(out), ".geno.desc")
  remove_bigmatrix(out)
  
  # parser map
  logging.log("Reading file...\n", verbose = verbose)
  m <- simer.Data.Map(bim_file, out = out, cols = c(2, 1, 4, 6, 5), header = FALSE)
  
  # parser phenotype, ind file
  fam <- read.table(fam_file, header = FALSE)
  n <- nrow(fam)
  write.table(fam[, 2], paste0( out, '.geno.ind'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  logging.log(paste0("inds: ", n, "\tmarkers:", m, '\n'), verbose = verbose)
  
  # parse genotype
  bigmat <- filebacked.big.matrix(
    nrow = m,
    ncol = n,
    type = type.geno,
    backingfile = backingfile,
    backingpath = dirname(out),
    descriptorfile = descriptorfile,
    dimnames = c(NULL, NULL)
  )
  
  if (priority == "speed") { maxLine <- -1 }
  read_bfile(bed_file = bed_file, pBigMat = bigmat@address, maxLine = maxLine, threads = threads, verbose = verbose)
  t2 <- as.numeric(Sys.time())
  logging.log("Preparation for GENOTYPE data is done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(invisible(c(m, n)))
}

#' simer.Data.Map: To check map file
#' 
#' checking map file.
#' 
#' Build date: Sep 12, 2018
#' Last update: July 25, 2022
#'
#' @author Haohao Zhang and Dong Yin
#' 
#' @param map the name of map file or map object(data.frame or matrix)
#' @param out the name of output file
#' @param cols selected columns
#' @param header whether the file contains header
#' @param sep seperator of the file
#' @param verbose whether to print detail.
#' 
#' @return 
#' Output file:
#' <out>.map
#' 
#' @export
#' 
#' @examples
#' # Get map path
#' mapPath <- system.file("extdata", "01bigmemory", "demo.geno.map", package = "simer")
#' 
#' # Check map data
#' simer.Data.Map(mapPath, tempfile("outfile"))
simer.Data.Map <- function(map, out = 'simer', cols = 1:5, header = TRUE, sep = '\t', verbose = TRUE) {
  t1 <- as.numeric(Sys.time())
  if (is.character(map) && !is.data.frame(map)) {
    map <- read.table(map, header = header, stringsAsFactors = FALSE)
  }
  map <- map[, cols]
  colnames(map) <- c("SNP", "CHROM", "POS", "REF", "ALT")
  if (length(unique(map[, 1])) != nrow(map)) {
    warning("WARNING: SNP is not unique and has been automatically renamed.")
    map[, 1] <- paste(map[, 2], map[, 3], sep = "-")
  }
  allels <- map[, 4:5]
  allels[allels == 0] <- '.'
  map[, 4:5] <- allels
  
  write.table(map, paste0(out, ".geno.map"), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  t2 <- as.numeric(Sys.time())
  logging.log("Preparation for MAP data is done within", format_time(t2 - t1), "\n", verbose = verbose)
  return(nrow(map))
}
