#' The main function of Data Quality Control
#'
#' Build date: Oct 19, 2020
#' Last update: May 26, 2021
#'
#' @author Dong Yin
#'
#' @param jsonFile the path of json file
#' @param verbose whether to print detail.
#'
#' @return report list
#' @export
#' @examples
#' jsonFile <- system.file("extdata", "demo1.json", package = "simer")
#' # JsonQC() needs set work directory at data directory 
#' # aa <- JsonQC(jsonFile = jsonFile)
JsonQC <- function(jsonFile, verbose=TRUE) {
  
  jsonList <- rjson::fromJSON(file = jsonFile)
  outpath <- getwd()
  out <- NULL
  
  # check genotype parameters
  genoPath <- jsonList$genotype
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
  filterGeno <- jsonList$quality_control_plan$genotype_quality_control$filter_geno
  filterHWE <- jsonList$quality_control_plan$genotype_quality_control$filter_hwe
  filterMind <- jsonList$quality_control_plan$genotype_quality_control$filter_mind
  filterMAF <- jsonList$quality_control_plan$genotype_quality_control$filter_maf
  
  # check pedigree parameters
  filePed <- jsonList$pedigree
  standardID <- jsonList$quality_control_plan$pedigree_quality_control$standard_ID
  fileSir <- jsonList$quality_control_plan$pedigree_quality_control$candidate_sire_file
  fileDam <- jsonList$quality_control_plan$pedigree_quality_control$candidate_dam_file
  exclThres <- jsonList$quality_control_plan$pedigree_quality_control$exclude_threshold
  assignThres <- jsonList$quality_control_plan$pedigree_quality_control$assign_threshold
  
  # check phenotype parameters
  filePhe <- sapply(jsonList$quality_control_plan$phenotype_quality_control, function(x) return(x$sample_info))
  planPhe <- jsonList$quality_control_plan$phenotype_quality_control

  ## step 1. data quality control
  simer.Data(fileMVP = fileMVP, fileBed = fileBed, filePlinkPed = filePlinkPed, genoType = 'char', filterGeno = filterGeno, filterHWE = filterHWE, filterMind = filterMind, filterMAF = filterMAF,
             filePhe = filePhe, planPhe = planPhe, pheCols = NULL, pheSep = '\t', missing=c(NA, 'NA', 'Na', '.', '-', 'NAN', 'nan', 'na', 'N/A', 'n/a', '<NA>', '', '-9', 9999),
             filePed = filePed, standardID = standardID, fileSir = fileSir, fileDam = fileDam, exclThres = exclThres, assignThres = assignThres, pedSep = '\t', 
             SNP.impute = "Major",
             outpath = outpath, out = out, maxLine = 10000, priority = "speed", ncpus = NULL, verbose = verbose)
  
}

#' The main function of Constructing Model and Index
#'
#' Build date: Aug 16, 2021
#' Last update: Aug 16, 2021
#'
#' @author Dong Yin
#'
#' @param jsonFile the path of json file
#' @param buildModel whether to build EBV model
#' @param buildIndex whether to build Selection Index
#' @param ncpus the number of threads
#' @param verbose whether to print detail.
#'
#' @return report list
#' @export
#' @examples
#' jsonFile <- system.file("extdata", "demo2.json", package = "simer")
#' # aa <- JsonModel(jsonFile = jsonFile)
JsonModel <- function(jsonFile, buildModel = TRUE, buildIndex = TRUE, ncpus = 10, verbose = TRUE) {
  
  jsonList <- rjson::fromJSON(file = jsonFile)
  outpath <- getwd()
  out <- NULL
  logging.initialize("Simer.Data", outpath)
  
  # check genotype parameters
  genoPath <- jsonList$genotype
  fileMVP <- NULL
  if (length(genoPath) != 0) {
    genoFiles <- list.files(genoPath)
    fileMVP <- fileBed <- fileNum <- NULL
    fileMVP <- grep(pattern = "geno.desc", genoFiles, value = TRUE)
    fileMVP <- file.path(genoPath, fileMVP)
    if (length(fileMVP) > 0) {
      fileMVP <- substr(fileMVP, 1, nchar(fileMVP)-10)
    } else {
      fileMVP <- NULL
    }
  }
  
  # check pedigree parameters
  filePed <- jsonList$pedigree

  # check phenotype parameters
  filePhe <- sapply(jsonList$analysis_plan, function(x) return(x$sample_info))
  planPhe <- jsonList$analysis_plan
  
  ## step 2. find the best effects for EBV model
  if (buildModel) {
    planPhe <- simer.Data.Env(planPhe = planPhe, fileMVP, filePed, ncpus = ncpus, verbose = verbose)
    jsonList$analysis_plan <- planPhe
  }
  
  newJsonFile <- paste0(substr(jsonFile, 1, nchar(jsonFile)-5), ".model.json")
  newJson <- rjson::toJSON(jsonList)
  cat(newJson, file = newJsonFile)

  ## step 3. construct selection index
  if (buildIndex) {
    BVIndex <- jsonList$breeding_value_index
    selIndex <- simer.Data.SELIND(BVIndex, planPhe, fileMVP, filePed, ncpus = ncpus, verbose = verbose)
    jsonList$selection_index <- selIndex
    jsonList$breeding_value_index <- NULL
  }
  
  newJson <- rjson::toJSON(jsonList)
  cat(newJson, file = newJsonFile)
  
}
