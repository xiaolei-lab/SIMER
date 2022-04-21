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


logging.initialize <- function(module, outpath) {
  file <- NULL
  if (options("simer.OutputLog2File") == TRUE) {
    now <- Sys.time()
    file <- paste(module, format(now, "%Y%m%d_%H%M%S"), "log", sep = ".")
    file <- file.path(outpath, file)
  }
  
  assign("logging.file", file, envir = package.env)
}

#' Print or write log
#'
#' @param ... R objects
#' @param file A connection, or a character string naming the file to print to. If "" (the default), cat prints to the standard output connection, the console unless redirected by sink. If it is "|cmd", the output is piped to the command given by ‘cmd’, by opening a pipe connection.
#' @param sep a character vector of strings to append after each element.
#' @param fill a logical or (positive) numeric controlling how the output is broken into successive lines.
#' @param labels character vector of labels for the lines printed. Ignored if fill is FALSE.
#' @param verbose whether to print detail
#'
#' @return print information
#' @export
#'
#' @examples
#' logging.log('simer')
logging.log <- function(..., file = NULL, sep = " ", fill = FALSE, labels = NULL, verbose = TRUE) {
  if (verbose) {
    cat(..., sep = sep, fill = fill, labels = labels)
  }
  
  if (is.null(file)) {
    try(file <- get("logging.file", envir = package.env), silent = TRUE)
  }
  
  if (!is.null(file)) {
    cat(..., file = file, sep = sep, fill = fill, labels = labels, append = TRUE)
  }
}

#' Print things into file
#'
#' Build date: Feb 7, 2020
#' Last update: Feb 7, 2020
#' by using base::print
#'
#' @author Dong Yin
#'
#' @param x a matrix or a list
#' @param file output file name
#' @param append ogical. If TRUE, output will be appended to file; otherwise, it will overwrite the contents of file
#' @param verbose whether to print details
#'
#' @return print in the screen
#' @export
#'
#' @examples
#' x <- list(a = "a", b = "b")
#' simer.print(x)
logging.print <- function(x, file = NULL, append = TRUE, verbose = TRUE) {
  if (verbose) {
    if (is.numeric(x) & is.vector(x)) {
      if (length(x) <= 10) {
        cat("", x, "\n")
      } else {
        cat("", x[1:10], "...(more IDs in the logging file)\n")
      }
    } else {
      print(x)
    }
    
    if (is.null(file)) {
      try(file <- get("logging.file", envir = package.env), silent = TRUE)
    }
    if (!is.null(file)) {
      sink(file = file, append = append)
      if (is.numeric(x) & is.vector(x)) {
        cat("", x, "\n")
      } else {
        print(x)
      }
      sink()
    }
  }
}
