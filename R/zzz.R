package.env <- NULL

.onLoad <- function(libname, pkgname) {
    # Limit number of threads in veclib (MacOS MRO)
    if (Sys.info()["sysname"] == "Darwin") {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    } 
    
    # set option
    op <- options()
    op.simer <- list(
        simer.OutputLog2File = TRUE
    )
    toset <- !(names(op.simer) %in% names(op))
    if (any(toset)) { 
        options(op.simer[toset])
    }
    
    # package level environment
    package.env <<- new.env()
    
    return(invisible())
}

.onUnload <- function(libpath) {
    options(simer.OutputLog2File = NULL)
}

.onAttach <- function(...){
    packageStartupMessage("Full description, Bug report, Suggestion and the latest version:")
    packageStartupMessage("http://git.soyomics.com:9000/yindong/simer")
}