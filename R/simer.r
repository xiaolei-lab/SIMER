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
#' Last update: Apr 17, 2022
#'
#' @author Dong Yin, Lilin Yin, Haohao Zhang and Xiaolei Liu
#'
#' @param SP a list of all simulation parameters
#' 
#' @return a list with simulated population
#' @export
#'
#' @examples
#' SP <- param.simer(out = "simer")
#' SP <- simer(SP)
simer <- function(SP) {

# Start simer

# TODO: How to generate inbreeding sirs and uninbreeding dams
# TODO: optcontri.sel  
# TODO: add superior limit of homo
# TODO: add summary() to population information
# TODO: add inbreeding coefficient
# TODO: update index selection
# TODO: add true block distribution  
# TODO: inbreeding change in every generations
# TODO: breeding literature from ZhenST
# TODO: remove one column genotype
  
  # unfold global parameters
  replication <- SP$global$replication
  seed.sim <- SP$global$seed.sim
  incols <- SP$global$incols
  outcols <- SP$global$outcols
  out <- SP$global$out
  outpath <- SP$global$outpath
  selPath <- SP$global$selPath
  out.format <- SP$global$out.format
  out.geno.gen <- SP$global$out.geno.gen
  out.pheno.gen <- SP$global$out.pheno.gen
  ncpus <- SP$global$ncpus
  verbose <- SP$global$verbose
  SP$geno$incols <- incols
  
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
  set.seed(seed.sim)
  if (incols == 1) outcols <- 1
  
  ################### DATA SIMULATION ###################
  SP <- annotation(SP = SP, verbose = verbose)
  SP <- genotype(SP = SP, ncpus = ncpus, verbose = verbose)
  SP <- phenotype(SP = SP, verbose = verbose)
  SP <- selects(SP = SP, verbose = verbose)
  SP <- reproduces(SP = SP, ncpus = ncpus, verbose = verbose)
  
  ################### DATA WRITING ###################
  write.file(SP)
  
  print_accomplished(width = 70, verbose = verbose)
  # Return the last directory
  ed <- Sys.time()
  logging.log(" SIMER DONE WITHIN TOTAL RUN TIME:", format_time(as.numeric(ed)-as.numeric(op)), "\n", verbose = verbose)
  return(SP)
}
