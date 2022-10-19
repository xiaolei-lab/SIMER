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


#' Simer
#' 
#' Main function of Simer.
#'
#' Build date: Jan 7, 2019
#' Last update: Apr 29, 2022
#'
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#'
#' @param SP a list of all simulation parameters.
#' 
#' @return 
#' the function returns a list containing
#' \describe{
#' \item{$global}{a list of global parameters.}
#' \item{$map}{a list of marker information parameters.}
#' \item{$geno}{a list of genotype simulation parameters.}
#' \item{$pheno}{a list of phenotype simulation parameters.}
#' \item{$sel}{a list of selection parameters.}
#' \item{$reprod}{a list of reproduction parameters.}
#' }
#' 
#' @export
#'
#' @examples
#' # Generate all simulation parameters
#' SP <- param.simer(out = "simer")
#' 
#' # Run Simer
#' SP <- simer(SP)
simer <- function(SP) {

### Start simer

# TODO: how to generate inbreeding sirs and uninbreeding dams
# TODO: optcontri.sel  
# TODO: add superior limit of homo
# TODO: add inbreeding coefficient
# TODO: add true block distribution
# TODO: inbreeding change in every generations
# TODO: breeding literature from ZhenST
# TODO: remove one column genotype
# TODO: data converter of simer
# TODO: core population
# TODO: pEBVs, gEBVs, and ssEBVs
  
  # global parameters
  replication <- SP$global$replication
  seed.sim <- SP$global$seed.sim
  out <- SP$global$out
  outpath <- SP$global$outpath
  ncpus <- SP$global$ncpus
  verbose <- SP$global$verbose
  
  # initialize logging
  if (!is.null(outpath)) {
    if (!dir.exists(outpath)) stop(paste0("Please check your output path: ", outpath))
    if (verbose) {
      logging.initialize("Simer", outpath = outpath)
    }
  }
  
  # welcome to simer
  simer.Version(width = 70, verbose = verbose) 
  
	################### MAIN_FUNCTION_SETTING ###################
  logging.log("--------------------------- replication ", replication, "---------------------------\n", verbose = verbose)
  op <- Sys.time()
  logging.log(" SIMER BEGIN AT", as.character(op), "\n", verbose = verbose)
  logging.log(" Random seed is", floor(seed.sim), "\n", verbose = verbose)
  set.seed(floor(seed.sim))
  
  ################### DATA SIMULATION ###################
  SP <- annotation(SP = SP, verbose = verbose)
  SP <- genotype(SP = SP, ncpus = ncpus, verbose = verbose)
  SP <- phenotype(SP = SP, verbose = verbose)
  SP <- selects(SP = SP, verbose = verbose)
  SP <- reproduces(SP = SP, ncpus = ncpus, verbose = verbose)
  
  ################### DATA WRITING ###################
  SP <- write.file(SP)
  
  print_accomplished(width = 70, verbose = verbose)
  ed <- Sys.time()
  logging.log(" SIMER DONE WITHIN TOTAL RUN TIME:", format_time(as.numeric(ed)-as.numeric(op)), "\n", verbose = verbose)
  
  return(SP)
}
