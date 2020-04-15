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


#' Print simer Banner
#'
#' Build date: Aug 30, 2017
#' Last update: Oct 21, 2019
#' 
#' @author Dong Yin, Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param width the width of the message
#' @param verbose whether to print detail.
#' 
#' @return version number.
#' @export
#'
#' @examples
#' simer.Version()
simer.Version <- function(width=60, verbose=TRUE) {
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


#' Print progress bar
#'
#' @param i the current loop number
#' @param n max loop number
#' @param type type1 for "for" function
#' @param symbol the symbol for the rate of progress
#' @param tmp.file the opened file of "fifo" function
#' @param symbol.head the head for the bar
#' @param symbol.tail the tail for the bar
#' @param fixed.points whether use the setted points which will be printed
#' @param points the setted points which will be printed
#' @param symbol.len the total length of progress bar
#'
#' @keywords internal
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
        #     if(inherits(parallel:::mcfork(), "masterProcess")) {
        #         progress <- 0.0
        #         while(progress < n && !isIncomplete(tmp.file)){
        #             msg <- readBin(tmp.file, "double")
        #             progress <- progress + as.numeric(msg)
        #             print.len <- round(symbol.len * progress / n)
        #             if(fixed.points){
        #                 if(progress %in% round(points * n / 100)){
        #                     logging.log(paste("\r", 
        #                               paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
        #                               paste(rep(" ", symbol.len-print.len), collapse=""),
        #                               sprintf("%.2f%%", progress * 100 / n), sep=""))
        #                 }
        #             }else{
        #                 logging.log(paste("\r", 
        #                           paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
        #                           paste(rep(" ", symbol.len-print.len), collapse=""),
        #                           sprintf("%.2f%%", progress * 100 / n), sep=""))
        #             }
        #         }
        #         parallel:::mcexit()
        #     }
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


print_accomplished <- function(width = 60, verbose = TRUE) {
    logging.log(make_line("SIMER ACCOMPLISHED", width = width, linechar = '='), "\n", verbose = verbose)
}

#' Print R Package information, include title, short_title, logo, version, authors, contact
#'
#' Build date: Oct 22, 2018
#' Last update: Oct 21, 2019
#' 
#' @keywords internal
#' @author Dong Yin and Haohao Zhang
#' 
#' @param welcome welcome text, for example: "Welcom to <Packagename>"
#' @param title long text to introduct package
#' @param short_title short label, top-left of logo
#' @param logo logo
#' @param version short label, bottom-right of logo
#' @param authors 
#' @param contact email or website
#' @param line 1, 2, or char
#' @param width banner width
#'
#' @export
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
print_info <- function(welcome=NULL, title=NULL, short_title=NULL, logo=NULL, version=NULL, authors=NULL, contact=NULL, linechar = '=', width=NULL, verbose=TRUE) {
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

#' make line
#' 
#' Build date: Dec 12, 2018
#' Last update: Dec 12, 2018
#' 
#' @keywords internal
#' @author Haohao Zhang
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

#' wrap text to multiple line, align left, right or center.
#' 
#' Build date: Oct 22, 2018
#' Last update: Dec 12, 2018
#' by using base::strwarp.
#' 
#' @keywords internal
#' @author Haohao Zhang
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

#' Paste label to a line
#' 
#' Build date: Oct 22, 2018
#' Last update: Oct 22, 2018
#' 
#' @param line long text
#' @param label short label
#' @param side "right" or "left"
#' @param margin default 2
#' 
#' @keywords internal
#' @author Haohao Zhang
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

#' format time
#' 
#' @param x seconds
#' 
#' @keywords internal
#' @export
#' @author Lilin Yin
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


load_if_installed <- function(package) {
    if (!identical(system.file(package = package), "")) {
        do.call('library', list(package))
        return(TRUE)
    } else {
        return(FALSE) 
    }
}


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


#' Remove big.matrix safely
#'
#' Build date: Aug 8, 2019
#' Last update: Feb 14, 2019
#'
#' @author Haohao Zhang and Dong Yin
#'
#' @param x filename of big.matrix
#' @param desc_suffix suffix of description file of big.matrix
#' @param bin_suffix suffix of binary file of big.matrix
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
#' options(bigmemory.typecast.warning=FALSE)
#'
#' remove_bigmatrix(x = "simer")
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
                }
            } else if (class(get(v, envir = envir)) == "list") {
                if (is.null(get(v, envir = envir)$geno)) next
                if (class(get(v, envir = envir)$geno) == "big.matrix") {
                    desc <- describe(get(v, envir = envir)$geno)@description
                    if (desc$filename == binfile) { # TODO: Risky deletion
                        rm(list = v, envir = envir)
                    }
                }
            }
        }
    }
    
    # Delete objects that occupy binfile in the global environment
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


#' Write files of simer
#'
#' Build date: Jan 7, 2019
#' Last update: Oct 16, 2019
#'
#' @author Dong Yin
#'
#' @param pop population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye
#' @param geno genotype matrix of population
#' @param map map information of markers
#' @param out.geno.index indice of individuals outputting genotype
#' @param out.pheno.index indice of individuals outputting phenotype
#' @param seed.map random seed of map file
#' @param out prefix of output file name
#' @param outpath path of output files
#' @param out.format format of output, "numeric" or "plink"
#' @param verbose whether to print detail
#'
#' @return None
#' @export
#'
#' @examples
#' \donttest{
#' data(simdata)
#' nmrk <- nrow(input.map)
#' pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
#' basepop <- getpop(nind = 100, from = 1, ratio = 0.1)
#' basepop.geno <- rawgeno
#' basepop.geno <- as.big.matrix(basepop.geno)
#' effs <-
#'     cal.effs(pop.geno = basepop.geno,
#'              cal.model = "A",
#'              num.qtn.tr1 = c(2, 6, 10),
#'              sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
#'              dist.qtn.tr1 = rep("normal", 6),
#'              prob.tr1 = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
#'              shape.tr1 = c(1, 1, 1, 1, 1, 1),
#'              scale.tr1 = c(1, 1, 1, 1, 1, 1),
#'              multrait = FALSE,
#'              num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
#'              sd.trn = diag(c(1, 0.5)),
#'              qtn.spot = rep(0.1, 10),
#'              maf = 0, 
#'              verbose = TRUE)
#' str(basepop)
#' pop.pheno <-
#'     phenotype(effs = effs,
#'               pop = basepop,
#'               pop.geno = basepop.geno,
#'               pos.map = NULL,
#'               h2.tr1 = c(0.3, 0.1, 0.05, 0.05, 0.05, 0.01),
#'               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
#'               h2.trn = c(0.3, 0.5),  
#'               sel.crit = "pheno", 
#'               pop.total = basepop, 
#'               sel.on = TRUE, 
#'               inner.env = NULL, 
#'               verbose = TRUE)
#' basepop <- pop.pheno$pop
#' pop.pheno$pop <- NULL           
#' idx <- basepop$index
#' seed.map <- 888888
#' # convert (0, 1) to (0, 1, 2)
#' basepop.geno <- geno.cvt(basepop.geno)
#' basepop.geno <- as.big.matrix(basepop.geno)
#' write.file(pop = basepop, geno = basepop.geno, map = pos.map, 
#'     out.geno.index = idx, out.pheno.index = idx, seed.map = seed.map, 
#'     outpath = tempdir(), out.format = "numeric", verbose = TRUE)
#' file.remove(file.path(tempdir(), "simer.geno.id"))
#' file.remove(file.path(tempdir(), "simer.map"))
#' file.remove(file.path(tempdir(), "simer.ped"))
#' file.remove(file.path(tempdir(), "simer.phe"))
#' }
write.file <- function(pop, geno, map, out.geno.index, out.pheno.index, seed.map, out = "simer", outpath, out.format, verbose) {
    if (is.null(outpath)) return(invisible())
    
    if (out.format == "numeric") {
        write.table(pop[out.geno.index, 2], file = file.path(outpath, paste0(out, ".geno.id")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        logging.log(" Generate genoid successfully!\n", verbose = verbose)
        write.table(map, file = file.path(outpath, paste0(out, ".map")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        logging.log(" Generate map successfully!\n", verbose = verbose)
        write.table(pop[out.pheno.index, c(2, 5, 6)], file = file.path(outpath, paste0(out, ".ped")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        logging.log(" Generate pedigree successfully!\n", verbose = verbose)
        write.table(pop[out.pheno.index, c(1, 2, 5, 7, 8:ncol(pop))], file = file.path(outpath, paste0(out, ".phe")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        logging.log(" Generate phenotype successfully!\n", verbose = verbose)
        
    } else if (out.format == "plink") {
        f1 <- grep(pattern = "TBV|TGV|pheno|ebv|u1", x = names(pop), value = TRUE)
        pheno <- subset(pop, select = c("index", f1))
        pheno <- pheno[out.geno.index, ]
        MVP.Data.MVP2Bfile(bigmat = geno, map = map, pheno = pheno, out = file.path(outpath, "mvp.plink"), verbose = verbose)
    }  
}

#' Build design matrix according to covariate
#'
#' Build date: Jan 20, 2020
#' Last update: Jan 20, 2020
#'
#' @author Dong Yin
#'
#' @param cv.name covariate names
#' @param data a data frame
#'
#' @return design matrix
#' @export
#'
#' @examples
#' dat <- data.frame(
#' f1 = sample(LETTERS[1:3], 20, TRUE),
#' f2 = sample(LETTERS[4:5], 20, TRUE),
#' row.names = paste0("id_", 1:20))
#' cv.name <- c("f1", "f2")
#' CV <- build.CV(cv.name = cv.name, data = dat)
build.CV <- function(cv.name, data) {
    design <- do.call(cbind, lapply(cv.name, function(cv){
        apply(outer(data[[cv]], unique(data[[cv]]), FUN = "=="), 2, as.integer)[, -1]
    }))
    rownames(design) <- rownames(data)
    colnames(design) <- unlist(sapply(cv.name, function(cv) unique(data[[cv]])[-1]))
    if (ncol(design) > 1)
        design <- design[, !duplicated(colnames(design))] # duplicated colnames happen sometimes
    return(design)
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
simer.print <- function(x, file = NULL, append = TRUE, verbose = TRUE) {
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

#' Show file content
#'
#' Build date: Feb 11, 2020
#' Last update: Feb 11, 2020
#'
#' @author Dong Yin
#'
#' @param filename filename of a file
#' @param verbose whether to print details
#'
#' @return none
#' @export
#'
#' @examples
#' selPath <- system.file("extdata", "01breeding_plan", package = "simer")
#' filename <- file.path(selPath, "breeding_plan01.txt")
#' simer.show.file(filename = filename, verbose = TRUE)
simer.show.file <- function(filename = NULL, verbose = TRUE) {
    if (is.null(filename)) stop("Please input a filename!")
    fileImage <- file(description=filename, open="r")
    inFile <- TRUE
    while (inFile) {
        tt <- readLines(fileImage, n=1)
        if (length(tt) == 0) {
            inFile = FALSE
            next
        }
        logging.log(tt, "\n", verbose = verbose)
    }
    close.connection(fileImage)
    return(invisible())
}