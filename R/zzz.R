package.env <- NULL

allparam <- c(
  "replication"  ,    "seed.sim"     ,    "out"          ,    "outpath"      ,
  "out.format"   ,    "pop.gen"      ,    "out.geno.gen" ,    "out.pheno.gen",
  "useAllGeno"   ,    "missing.geno" ,    "missing.phe"  ,    "ncpus"        ,
  "verbose"      ,    "pop.map"      ,    "species"      ,    "pop.marker"   ,
  "num.chr"      ,    "len.chr"      ,    "qtn.model"    ,    "qtn.index"    ,
  "qtn.num"      ,    "qtn.dist"     ,    "qtn.var"      ,    "qtn.prob"     ,
  "qtn.shape"    ,    "qtn.scale"    ,    "qtn.shape1"   ,    "qtn.shape2"   ,
  "qtn.ncp"      ,    "qtn.spot"     ,    "len.block"    ,    "maf"          ,
  "recom.spot"   ,    "range.hot"    ,    "range.cold"   ,    "pop.geno"     ,
  "incols"       ,    "pop.marker"   ,    "pop.ind"      ,    "prob"         ,
  "rate.mut"     ,    "cld"          ,    "pop"          ,    "pop.ind"      ,
  "pop.rep"      ,    "pop.rep.bal"  ,    "pop.env"      ,    "phe.type"     ,
  "phe.model"    ,    "phe.h2A"      ,    "phe.h2D"      ,    "phe.h2GxG"    ,
  "phe.h2GxE"    ,    "phe.h2PE"     ,    "phe.var"      ,    "phe.corA"     ,
  "phe.corD"     ,    "phe.corGxG"   ,    "phe.corPE"    ,    "phe.corE"     ,
  "pop.sel"      ,    "ps"           ,    "decr"         ,    "sel.crit"     ,
  "sel.single"   ,    "sel.multi"    ,    "index.wt"     ,    "index.tdm"    ,
  "goal.perc"    ,    "pass.perc"    ,    "pop.gen"      ,    "reprod.way"   ,
  "sex.rate"     ,    "prog"         
)

.onLoad <- function(libname, pkgname) {
  # limit number of threads in veclib (MacOS MRO)
  if (Sys.info()["sysname"] == "Darwin") {
    Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
  } 
  
  # set option
  op <- options()
  op.simer <- list(
    simer.OutputLog2File = TRUE,
    simer.show.warning = TRUE
  )
  toset <- !(names(op.simer) %in% names(op))
  if (any(toset)) { 
    options(op.simer[toset])
  }
  options(bigmemory.typecast.warning = FALSE)
  
  # package level environment
  package.env <<- new.env()
  
  return(invisible())
}

.onUnload <- function(libpath) {
  options(simer.OutputLog2File = NULL)
}

.onAttach <- function(...){
  packageStartupMessage("Full description, Bug report, Suggestion and the latest version:")
  packageStartupMessage("https://github.com/xiaolei-lab/SIMER")
}