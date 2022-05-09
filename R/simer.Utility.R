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


#' Simer version
#' 
#' Print simer version.
#'
#' Build date: Aug 30, 2017
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param width the width of the message.
#' @param verbose whether to print detail.
#' 
#' @return version number.
#' 
#' @export
#'
#' @examples
#' simer.Version()
simer.Version <- function(width = 60, verbose = TRUE) {
  welcome <- "Welcome to SIMER"
  title   <- "Data Simulation for Life Science and Breeding"
  authors <- c("Designed and Maintained by Dong Yin, Xuanning Zhang, Lilin Yin, Haohao Zhang, and Xiaolei Liu", 
               "Contributors: Zhenshuang Tang, Jingya Xu, Xinyun Li, Mengjin Zhu, Xiaohui Yuan, and Shuhong Zhao")
  contact <- "Contact: xiaoleiliu@mail.hzau.edu.cn"
  logo_s  <- c(" ____ ___ __  __ _____ ____  ", 
               "/ ___|_ _|  \\/  | ____|  _ \\ ", 
               "\\___ \\| || |\\/| |  _| | |_) |", 
               " ___) | || |  | | |___|  _ < ", 
               "|____/___|_|  |_|_____|_| \\_\\")
                              
  version <- print_info(welcome = welcome, title = title, logo = logo_s, authors = authors, contact = contact, linechar = '=', width = width, verbose = verbose)
  return(invisible(version))
}

#' Progress bar
#' 
#' Print progress bar.
#' 
#' Build date: Aug 30, 2017
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param i the current loop number.
#' @param n the max loop number.
#' @param type type1 for "for" function.
#' @param symbol the symbol for the rate of progress.
#' @param tmp.file the opened file of "fifo" function.
#' @param symbol.head the head for the bar.
#' @param symbol.tail the tail for the bar.
#' @param fixed.points whether use the setted points which will be printed.
#' @param points the setted points which will be printed.
#' @param symbol.len the total length of progress bar.
#' @param verbose whether to print detail.
#'
#' @keywords internal
#' 
#' @return none.
print_bar <- function(i,
                      n,
                      type = c("type1", "type3"),
                      symbol = "-",
                      tmp.file = NULL,
                      symbol.head = ">>>",
                      symbol.tail = ">" ,
                      fixed.points = TRUE,
                      points = seq(0, 100, 1),
                      symbol.len = 48,
                      verbose = TRUE
) {
  switch(
    match.arg(type), 
    "type1"={
      if(fixed.points){
        point.index <- points
        point.index <- point.index[point.index > floor(100*(i-1)/n)]
        if(floor(100*i/n) %in% point.index){
          if(floor(100*i/n) != max(point.index)){
            print.len <- floor(symbol.len*i/n)
            logging.log(
              paste("\r", 
                    paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                    paste(rep(" ", symbol.len-print.len), collapse=""),
                    sprintf("%.2f%%", 100*i/n)
                    , sep=""),
              verbose = verbose
            )
          }else{
            print.len <- floor(symbol.len*i/n)
            logging.log(
              paste("\r", 
                paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                sprintf("%.2f%%", 100*i/n), "\n"
                , sep=""),
              verbose = verbose
            )
          }
        }
      }else{
        if(i < n){
          print.len <- floor(symbol.len*i/n)
          logging.log(
            paste("\r", 
                  paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                  paste(rep(" ", symbol.len-print.len), collapse=""),
                  sprintf("%.2f%%", 100*i/n)
                  , sep=""),
            verbose = verbose
          )
        }else{
          print.len <- floor(symbol.len*i/n)
          logging.log(
            paste("\r", 
                  paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                  sprintf("%.2f%%", 100*i/n), "\n"
                  , sep=""),
            verbose = verbose
          )
        }
      }
    },
    # "type2"={
    #    if(inherits(parallel:::mcfork(), "masterProcess")) {
    #      progress <- 0.0
    #      while(progress < n && !isIncomplete(tmp.file)){
    #        msg <- readBin(tmp.file, "double")
    #        progress <- progress + as.numeric(msg)
    #        print.len <- round(symbol.len * progress / n)
    #        if(fixed.points){
    #          if(progress %in% round(points * n / 100)){
    #              logging.log(paste("\r", 
    #                        paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
    #                        paste(rep(" ", symbol.len-print.len), collapse=""),
    #                        sprintf("%.2f%%", progress * 100 / n), sep=""))
    #          }
    #        }else{
    #          logging.log(paste("\r", 
    #                    paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
    #                    paste(rep(" ", symbol.len-print.len), collapse=""),
    #                    sprintf("%.2f%%", progress * 100 / n), sep=""))
    #        }
    #      }
    #      parallel:::mcexit()
    #    }
    # },
    "type3"={
      progress <- readBin(tmp.file, "double") + 1
      writeBin(progress, tmp.file)
      print.len <- round(symbol.len * progress / n)
      if(fixed.points){
        if(progress %in% round(points * n / 100)){
          logging.log(
            paste("\r", 
                  paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                  paste(rep(" ", symbol.len-print.len), collapse=""),
                  sprintf("%.2f%%", progress * 100 / n)
                  , sep=""),
            verbose = verbose
          )
        }
      }else{
        logging.log(
          paste("\r", 
                paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                paste(rep(" ", symbol.len-print.len), collapse=""),
                sprintf("%.2f%%", progress * 100 / n)
                , sep=""),
          verbose = verbose
        )
      }
    }
  )
}

#' Accomplishment
#' 
#' Print accomplishment information.
#' 
#' Build date: Aug 30, 2017
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param width the width of the message.
#' @param verbose whether to print detail.
#'
#' @keywords internal
#' 
#' @return none.
print_accomplished <- function(width = 60, verbose = TRUE) {
  logging.log(make_line("SIMER ACCOMPLISHED", width = width, linechar = '='), "\n", verbose = verbose)
}

#' Simer information
#' 
#' Print R Package information, include title, short_title, logo, version, authors, contact.
#'
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin and Haohao Zhang
#' 
#' @param welcome welcome text, for example: "Welcom to <Packagename>".
#' @param title long text to introduct package.
#' @param short_title short label, top-left of logo.
#' @param logo logo.
#' @param version short label, bottom-right of logo.
#' @param authors authors of software.
#' @param contact email or website.
#' @param line 1, 2, or char.
#' @param width banner width.
#' @param verbose whether to print detail.
#'
#' @export
#' 
#' @return welcome information.
#' 
#' @keywords internal
#'
#' @examples
#' welcome <- "Welcome to SIMER"
#' title   <- "Data Simulation for Life Science and Breeding"
#' authors <- c("Designed and Maintained by Dong Yin, Xuanning Zhang,
#'               Lilin Yin, Haohao Zhang, and Xiaolei Liu", 
#'              "Contributors: Zhenshuang Tang, Jingya Xu, Xinyun Li, 
#'               Mengjin Zhu, Xiaohui Yuan, Shuhong Zhao")
#' contact <- "Contact: xiaoleiliu@mail.hzau.edu.cn"
#' logo_s  <- c(" ____ ___ __  __ _____ ____  ", 
#'              "/ ___|_ _|  \\/  | ____|  _ \\ ", 
#'              "\\___ \\| || |\\/| |  _| | |_) |", 
#'              " ___) | || |  | | |___|  _ < ", 
#'              "|____/___|_|  |_|_____|_| \\_\\")
#' print_info(welcome = welcome, title = title, logo = logo_s, authors = authors, 
#'            contact = contact, linechar = '=', width = 70)
print_info <- function(welcome = NULL, title = NULL, short_title = NULL, logo = NULL, version = NULL, authors = NULL, contact = NULL, linechar = '=', width = NULL, verbose = TRUE) {
  msg <- c()
  # width
  if (is.null(width)) { width <- getOption('width') }
  # version
  if (is.null(version)) {
    if (getPackageName() == ".GlobalEnv") {
      version <- "devel"
    } else {
      version <- as.character(packageVersion(getPackageName()))
    }
  }
  # welcome
  if (is.null(welcome)) { 
    if (getPackageName() == ".GlobalEnv") {
      welcome <- ""
    } else {
      welcome <- paste0("Welcome to ", getPackageName())
    }
  }
  msg <- c(msg, make_line(welcome, linechar = linechar, width = width))
  # title
  if (!is.null(title)) {
    msg <- c(msg, rule_wrap(string = title, width = width, align = "center"))
  }
  
  # align logo
  logo_width <- max(sapply(logo, nchar))
  for (i in 1:length(logo)) {
    l <- paste0(logo[i], paste(rep(" ", logo_width - nchar(logo[i])), collapse = ""))
    l <- make_line(l, width)
    msg <- c(msg, l)
  }
  
  # paste short_title label to logo top-left
  if (!is.null(short_title)) {
    i <- length(msg) - length(logo) + 1
    msg[i] <- paste_label(msg[i], paste0(short_title), side = "left")
  }
  
  # paste version label to logo bottom-right
  msg[length(msg)] <- paste_label(msg[length(msg)], paste0("Version: ", version), side = "right")
  
  # authors
  if (!is.null(authors)) {
    msg <- c(msg, rule_wrap(string = authors, align = "left", linechar = " ", width = width))
  }
  # contact
  if (!is.null(contact)) {
    msg <- c(msg, rule_wrap(string = contact, align = "left", linechar = " ", width = width))
  }
  # bottom line
  msg <- c(msg, paste0(rep(linechar, width), collapse = ''))
  
  logging.log(msg, sep = "\n", verbose = verbose)
  
  return(version)
}

#' Line making
#' 
#' Add a line to the screen.
#' 
#' Build date: Dec 12, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param string a string.
#' @param width the width of the message.
#' @param linechar char in every line.
#' @param align the position of string.
#' @param margin the margin information, default 2.
#'
#' @keywords internal
#' 
#' @return none.
make_line <- function(string, width, linechar = " ", align = "center", margin = 1) {
  string <- paste0(paste0(rep(" ", margin), collapse = ""),
                   string,
                   paste0(rep(" ", margin), collapse = ""))
  
  if (align == "center") {
    if (width > nchar(string)) {
      left_width <- (width - nchar(string)) %/% 2
      right_width <- width - nchar(string) - left_width
      string <-
        paste0(paste0(rep(linechar, left_width), collapse = ""),
               string,
               paste0(rep(linechar, right_width), collapse = ""))
    }
  } else if (align == "left") {
    if (width > nchar(string)) {
      string <-
        paste0(linechar,
               string,
               paste0(rep(linechar, width - nchar(string) - 1), collapse = ""))
    }
  }
  return(string)
}

#' Accomplishment
#' 
#' Print accomplishment information.
#' 
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param string a string.
#' @param width the width of the message.
#' @param align the position of string.
#' @param linechar char in every line.
#'
#' @keywords internal
#' 
#' @return none.
rule_wrap <- function(string, width, align = "center", linechar = " ") {
  # define
  msg <- c()
  lines <- strwrap(string, width = width - 4)
  
  # wrap
  for (i in 1:length(lines)) {
    l <- make_line(lines[i], width = width, linechar = linechar, align = align)
    msg <- c(msg, l)
  }
  return(msg)
}

#' Pasting label
#' 
#' Paste label to a line.
#' 
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param line long text.
#' @param label short label.
#' @param side "right" or "left".
#' @param margin the margin information, default 2.
#' 
#' @keywords internal
#' 
#' @return none.
paste_label <- function(line, label, side = "right", margin = 2) {
  if (side == "right") {
    end   <- nchar(line) - margin
    start <- end - (nchar(label) - 1)
  } else {
    start <- 1 + margin
    end   <- start + (nchar(label) - 1)
  }
  substr(line, start, end) <- label
  return(line)
}

#' Time formating
#' 
#' Format the time.
#' 
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param x the total seconds.
#' 
#' @keywords internal
#' 
#' @export
#' 
#' @return running time.
#'
#' @examples
#' format_time(x = 7200)
format_time <- function(x) {
  h <- x %/% 3600
  m <- (x %% 3600) %/% 60
  s <- ((x %% 3600) %% 60)
  index <- which(c(h, m, s) != 0)
  num <- c(h, m, s)[index]
  num <- round(num, 0)
  char <- c("h", "m", "s")[index]
  return(paste0(num, char, collapse = ""))
}


#' Installation checking
#' 
#' Check if the software is installed.
#' 
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param package the package name.
#' 
#' @keywords internal
#' 
#' @return none.
load_if_installed <- function(package) {
  if (!identical(system.file(package = package), "")) {
    do.call('library', list(package))
    return(TRUE)
  } else {
    return(FALSE) 
  }
}


#' MKL environment
#' 
#' Run code in the MKL environment.
#' 
#' Build date: Oct 22, 2018
#' Last update: Apr 30, 2022
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param exprs the expression.
#' @param threads the number of threads used, if NULL, (logical core number - 1) is automatically used.
#' 
#' @keywords internal
#' 
#' @return none.
mkl_env <- function(exprs, threads = 1) {
  if (load_if_installed("RevoUtilsMath")) {
    math.cores <- eval(parse(text = "getMKLthreads()"))
    eval(parse(text = "setMKLthreads(threads)"))
  }
  result <- exprs
  if (load_if_installed("RevoUtilsMath")) {
    eval(parse(text = "setMKLthreads(math.cores)"))
  }
  return(result)
}

#' Big.matrix removing
#' 
#' Remove big.matrix safely.
#'
#' Build date: Aug 8, 2019
#' Last update: Apr 30, 2022
#'
#' @author Haohao Zhang and Dong Yin
#'
#' @param x the filename of big.matrix.
#' @param desc_suffix the suffix of description file of big.matrix.
#' @param bin_suffix the suffix of binary file of big.matrix.
#'
#' @return TRUE or FALSE
#'
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
#'      backingfile = 'simer.geno.bin',
#'      descriptorfile = 'simer.geno.desc')
#'
#' remove_bigmatrix(x = "simer")
remove_bigmatrix <- function(x, desc_suffix = ".geno.desc", bin_suffix = ".geno.bin") {
  name <- basename(x)
  path <- dirname(x)
  
  descfile <- paste0(x, desc_suffix)
  binfile  <- paste0(x, bin_suffix)
  
  remove_var <- function(binfile, envir) {
    for (v in ls(envir = envir)) {
      if (any(class(get(v, envir = envir)) == "big.matrix")) {
        desc <- describe(get(v, envir = envir))@description
        if (desc$filename == binfile) {
          rm(list = v, envir = envir)
          gc()
        }
      }
    }
  }
  
  # delete objects that occupy binfile in the global environment
  remove_var(binfile, as.environment(-1L))
  remove_var(binfile, globalenv())
  gc()
  
  if (file.exists(descfile)) {
   file.remove(descfile)
  }
  if (file.exists(binfile)) {
    file.remove(binfile)
  }
}

#' File writing
#' 
#' Write files of Simer.
#'
#' Build date: Jan 7, 2019
#' Last update: Apr 30, 2022
#'
#' @author Dong Yin
#'
#' @param SP a list of all simulation parameters.
#'
#' @return none.
#' 
#' @export
#'
#' @examples
#' outpath <- tempdir()
#' SP <- param.simer(out = "simer")
#' SP <- simer(SP)
#' SP$global$outpath <- outpath
#' write.file(SP)
#' unlink(file.path(outpath, "180_Simer_Data_numeric"), recursive = TRUE)
write.file <- function(SP) {
  
  # global parameters
  replication <- SP$global$replication
  out <- SP$global$out
  outpath <- SP$global$outpath
  out.format <- SP$global$out.format
  out.geno.gen <- SP$global$out.geno.gen
  out.pheno.gen <- SP$global$out.pheno.gen
  verbose <- SP$global$verbose
  incols <- SP$geno$incols
  
  if (is.null(outpath)) return(SP)
  
  pop.marker <- nrow(SP$geno$pop.geno[[1]])
  pop.inds <- sapply(out.geno.gen, function(i) {
    return(ncol(SP$geno$pop.geno[[i]]) / incols)
  })[out.geno.gen]
  pop.ind <- sum(pop.inds)
  
  if (max(out.geno.gen) > SP$reprod$pop.gen) {
    stop("'out.geno.gen' should be not more than 'pop.gen'!")
  }
  if (max(out.pheno.gen) > SP$reprod$pop.gen) {
    stop("'out.pheno.gen' should be not more than 'pop.gen'!")
  }
  
  if (out.format == "numeric") {
    outpath = paste0(outpath, .Platform$file.sep, pop.ind, "_Simer_Data_numeric")
  } else if (out.format == "plink"){
    outpath = paste0(outpath, .Platform$file.sep, pop.ind, "_Simer_Data_plink")
  } else {
    stop("out.format should be 'numeric' or 'plink'!")
  }
  if (!dir.exists(outpath)) dir.create(outpath)
  
  directory.rep <- paste0(outpath, .Platform$file.sep, "replication", replication)
  if (dir.exists(directory.rep)) {
    remove_bigmatrix(file.path(directory.rep, out))
    unlink(directory.rep, recursive = TRUE)
  }
  dir.create(directory.rep)
  
  geno.back <- paste0(out, ".geno.bin")
  geno.desc <- paste0(out, ".geno.desc")
  geno.total <- filebacked.big.matrix(
    nrow = pop.marker,
    ncol = pop.ind,
    init = 3,
    type = 'char',
    backingpath = directory.rep,
    backingfile = geno.back,
    descriptorfile = geno.desc
  )
  
  pop.inds <- Reduce("+", pop.inds, accumulate = TRUE)
  pop.inds <- c(0, pop.inds[-length(pop.inds)]) + 1
  for (i in 1:length(out.geno.gen)) {
    if (incols == 1) {
      BigMat2BigMat(geno.total@address, SP$geno$pop.geno[[out.geno.gen[i]]]@address, op = pop.inds[i])
    } else {
      Mat2BigMat(geno.total@address, geno.cvt1(SP$geno$pop.geno[[out.geno.gen[i]]][]), op = pop.inds[i])
    }
  }
  
  SP$global$useAllGeno <- TRUE
  SP <- phenotype(SP)
  pheno.geno <- NULL
  pheno.total <- do.call(rbind, lapply(out.pheno.gen, function(i) {
    return(SP$pheno$pop[[i]])
  }))
  
  if (out.format == "numeric") {
    write.table(pheno.total[, 1], file = file.path(directory.rep, paste0(out, ".geno.id")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(SP$map$pop.map, file = file.path(directory.rep, paste0(out, ".geno.map")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    if (!is.null(SP$map$pop.map.GxG)) {
      write.table(SP$map$pop.map.GxG, file = file.path(directory.rep, paste0(out, ".GxG.geno.map")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
  
  } else if (out.format == "plink") {
    pheno.geno <- do.call(rbind, lapply(out.geno.gen, function(i) {
      return(SP$pheno$pop[[i]])
    }))
    MVP.Data.MVP2Bfile(bigmat = geno.total, map = SP$map$pop.map, pheno = pheno.geno[, 1, drop = FALSE], out = file.path(directory.rep, out), verbose = verbose)
    remove_bigmatrix(file.path(directory.rep, out))
    geno.total <- 0
  }
  write.table(pheno.total[, c(1, 5, 6)], file = file.path(directory.rep, paste0(out, ".ped")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(pheno.total, file = file.path(directory.rep, paste0(out, ".phe")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  logging.log(" All files have been saved successfully!\n", verbose = verbose)
  
  rm(geno.total); rm(pheno.total); rm(pheno.geno); gc()
  
  return(SP)
}
