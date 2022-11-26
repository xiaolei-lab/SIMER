library(simer)
# install.packages("AlphaSimR")
library(AlphaSimR)

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

# create genetic map for one chromosomes, each 100 Morgan long
genMap <- list(seq(0, 100, length.out = nmrk))

# create haplotypes for outbred individuals
geno1 <- geno.cvt2(geno[])
haplotypes <- list(chr1 = t(geno1))

# creating Founder Haplotypes
founderPop <- newMapPop(genMap = genMap, haplotypes = haplotypes)

fileHead <- c("h2A_tr1", "h2D_tr1", "h2A_tr2", "h2D_tr2", "h2A_tr3", "h2D_tr3")
cat(paste0(paste(fileHead, collapse=","), "\n"), file="res_AlphaSimR_AD_trait.csv", append=FALSE)
for (i in 1:rep) {
    set.seed(seeds[i])
    cat("###### This is the", i, "Replication ######\n")
    cat("###### The random seed is", seeds[i], "######\n")
    
    # T1
    SP = SimParam$new(founderPop)
    SP$addTraitAD(nQtlPerChr = nQTN, mean = 0, var = 20, meanDD = 1, useVarA = FALSE)
    SP$setVarE(H2 = 0.2)
    pop <- newPop(founderPop, simParam = SP)
    pheno1 <- pheno(pop)

    # T2
    SP = SimParam$new(founderPop)
    SP$addTraitAD(nQtlPerChr = nQTN, mean = 0, var = 50, meanDD = 1, useVarA = FALSE)
    SP$setVarE(H2 = 0.5)
    pop <- newPop(founderPop, simParam = SP)
    pheno2 <- pheno(pop)

    # T3
    SP = SimParam$new(founderPop)
    SP$addTraitAD(nQtlPerChr = nQTN, mean = 0, var = 80, meanDD = 1, useVarA = FALSE)
    SP$setVarE(H2 = 0.8)
    pop <- newPop(founderPop, simParam = SP)
    pheno3 <- pheno(pop)
    
    pheno  <- data.frame(index = 1:nind, T1 = pheno1[, 1], T2 = pheno2[, 1], T3 = pheno3[, 1])    
    
    # write phenotype
    write.table(pheno, paste0("AlphaSimR_AD_trait.phe.rep", i), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)
    
    res <- NULL
    
    # genetic evaluation for T1
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("AlphaSimR_AD_trait.phe.rep", i)),
        paste("--pheno-pos 2"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("AlphaSimR_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("AlphaSimR_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    # genetic evaluation for T2
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("AlphaSimR_AD_trait.phe.rep", i)),
        paste("--pheno-pos 3"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("AlphaSimR_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("AlphaSimR_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    # genetic evaluation for T3
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("AlphaSimR_AD_trait.phe.rep", i)),
        paste("--pheno-pos 4"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add --dom", "--threads 10"),
        paste("--out", paste0("AlphaSimR_AD_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("AlphaSimR_AD_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    
    cat(paste0(paste(res, collapse=","), "\n"), file="res_AlphaSimR_AD_trait.csv", append=TRUE)
    
}

res <- read.csv('res_AlphaSimR_AD_trait.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_AlphaSimR_AD_trait.csv', quote = FALSE)
