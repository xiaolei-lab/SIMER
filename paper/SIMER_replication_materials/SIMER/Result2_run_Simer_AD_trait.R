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

pheno_col <- 2:4
pop.map <- generate.map(pop.marker = nmrk)
qtn.num <- nQTN

fileHead <- c("h2A_tr1", "h2D_tr1", "h2A_tr2", "h2D_tr2", "h2A_tr3", "h2D_tr3")
cat(paste0(paste(fileHead, collapse=","), "\n"), file="res_Simer_AD_trait.csv", append = FALSE)
for (i in 1:rep) {
    set.seed(seeds[i])
    cat("###### This is the", i, "Replication ######\n")
    cat("###### The random seed is", seeds[i], "######\n")
    
    # T1
    SP <- param.annot(pop.map = pop.map, qtn.num = qtn.num, qtn.model = "A + D")
    SP <- param.geno(SP = SP, pop.geno = geno)
    SP <- param.pheno(
      SP = SP,
      pop.ind = nind,
      phe.model = list(
        tr1 = "T1 = A + D + E" # "T1" (Trait 1) consists of Additive effect, Dominant effect, and Residual effect
      ),
      phe.var = list(tr1 = 100),
      phe.h2A = list(tr1 = 0.2 * 0.9),
      phe.h2D = list(tr1 = 0.2 * 0.1)
    )
    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)
    pheno1 <- SP$pheno$pop$gen1[, 8, drop = FALSE]

    # T2
    SP <- param.annot(pop.map = pop.map, qtn.num = qtn.num, qtn.model = "A + D")
    SP <- param.geno(SP = SP, pop.geno = geno)
    SP <- param.pheno(
      SP = SP,
      pop.ind = nind,
      phe.model = list(
        tr1 = "T1 = A + D + E" # "T1" (Trait 1) consists of Additive effect, Dominant effect, and Residual effect
      ),
      phe.var = list(tr1 = 100),
      phe.h2A = list(tr1 = 0.5 * 0.9),
      phe.h2D = list(tr1 = 0.5 * 0.1)
    )
    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)
    pheno2 <- SP$pheno$pop$gen1[, 8, drop = FALSE]
    
    # T3
    SP <- param.annot(pop.map = pop.map, qtn.num = qtn.num, qtn.model = "A + D")
    SP <- param.geno(SP = SP, pop.geno = geno)
    SP <- param.pheno(
      SP = SP,
      pop.ind = nind,
      phe.model = list(
        tr1 = "T1 = A + D + E" # "T1" (Trait 1) consists of Additive effect, Dominant effect, and Residual effect
      ),
      phe.var = list(tr1 = 100),
      phe.h2A = list(tr1 = 0.8 * 0.9),
      phe.h2D = list(tr1 = 0.8 * 0.1)
    )
    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)
    pheno3 <- SP$pheno$pop$gen1[, 8, drop = FALSE]

    pheno  <- data.frame(index = 1:nind, T1 = pheno1[, 1], T2 = pheno2[, 1], T3 = pheno3[, 1])
    
    # write phenotype
    write.table(pheno, paste0("Simer_AD_trait.phe.rep", i), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)
    
    res <- NULL
    
    # genetic evaluation for T1
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("Simer_AD_trait.phe.rep", i)),
        paste("--pheno-pos 2"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("Simer_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("Simer_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    # genetic evaluation for T2
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("Simer_AD_trait.phe.rep", i)),
        paste("--pheno-pos 3"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("Simer_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("Simer_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    # genetic evaluation for T3
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("Simer_AD_trait.phe.rep", i)),
        paste("--pheno-pos 4"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("Simer_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("Simer_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    cat(paste0(paste(res, collapse=","), "\n"), file="res_Simer_AD_trait.csv", append=TRUE)
    
}

res <- read.csv('res_Simer_AD_trait.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_Simer_AD_trait.csv', quote = FALSE)
