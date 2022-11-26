library(simer)

# get genotype
prefix <- "pig_small"
path <- "../PICData/"
filename <- paste0(file.path(path, prefix), ".geno.desc")
geno <- attach.big.matrix(filename)
nmrk <- nrow(geno)
nind <- ncol(geno)

rep <- 100
nQTN <- 1000
seeds <- 1:rep

pop.map <- generate.map(pop.marker = nmrk)
qtn.num <- list(tr1 = nQTN)

fileHead <- c("progress1", "progress2", "progress3")
cat(paste0(paste(fileHead, collapse=","), "\n"), file="res_Simer_breeding_program.csv", append = FALSE)
for (i in 1:rep) {
    set.seed(seeds[i])
    cat("###### This is the", i, "Replication ######\n")
    cat("###### The random seed is", seeds[i], "######\n")
    
    # set Simulation Parameters
    SP <- param.annot(pop.map = pop.map, qtn.num = qtn.num, qtn.model = "A")
    SP <- param.geno(SP = SP, pop.geno = geno)
    SP <- param.pheno(
      SP = SP,
      pop.ind = nind,
      phe.model = list(tr1 = "T1 = A + E"),
      phe.var = list(tr1 = 100),
      phe.h2A = list(tr1 = 0.8)
    )

    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)

    # get phenotype
    pheno <- SP$pheno$pop$gen1
    pheno[, 2] <- cut(1:nrow(pheno), 10, labels = FALSE)

    # write phenotype
    write.table(pheno, paste0("Simer_breeding_program.phe"), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)
    
    # plan 1
    jsonFile <- "plan1.json"
    jsonList <- simer.Data.Json(jsonFile = jsonFile, hiblupPath = '../HIBLUP/', out = "./plan1/Simer_breeding_program_plan1")
    progress1 <- as.numeric(unlist(strsplit(jsonList$selection_index, split = "\\s+=\\s+"))[2])
    
    # plan 2
    jsonFile <- "plan2.json"
    jsonList <- simer.Data.Json(jsonFile = jsonFile, hiblupPath = '../HIBLUP/', out = "./plan2/Simer_breeding_program_plan2")
    progress2 <- as.numeric(unlist(strsplit(jsonList$selection_index, split = "\\s+=\\s+"))[2])

    # plan 3
    jsonFile <- "plan3.json"
    jsonList <- simer.Data.Json(jsonFile = jsonFile, hiblupPath = '../HIBLUP/', out = "./plan3/Simer_breeding_program_plan3")
    progress3 <- as.numeric(unlist(strsplit(jsonList$selection_index, split = "\\s+=\\s+"))[2])
    
    res <- c(progress1, progress2, progress3)

    cat(paste0(paste(res, collapse=","), "\n"), file="res_Simer_breeding_program.csv", append=TRUE)

}

res <- read.csv('res_Simer_breeding_program.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_Simer_breeding_program.csv', quote = FALSE)
