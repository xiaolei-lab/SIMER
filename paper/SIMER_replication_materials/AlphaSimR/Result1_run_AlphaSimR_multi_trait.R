library(simer)
# install.packages("AlphaSimR")
library(AlphaSimR)

# get genotype
prefix <- "cad1k_small"
path <- "../HumanData/"
filename <- paste0(file.path(path, prefix), ".geno.desc")
geno <- attach.big.matrix(filename)
nmrk <- nrow(geno)
nind <- ncol(geno)

rep <- 100
nQTN <- 1000
seeds <- 1:rep

pheno_col <- 2:4
corA <- matrix(c(1.0, 0.2, 0.5,
                 0.2, 1.0, 0.8,
                 0.5, 0.8, 1.0), 3, 3, byrow = TRUE)

# create genetic map for one chromosomes, 100 Morgan long
genMap = list(seq(0, 100, length.out = nmrk))

# create haplotypes for outbred individuals
geno1 <- geno.cvt2(geno[])
haplotypes <- list(chr1 = t(geno1))

# create Founder Haplotypes
founderPop <- newMapPop(genMap = genMap, haplotypes = haplotypes)

fileHead <- c("h2_tr1", "h2_tr2", "h2_tr3", "cor_tr12", "cor_tr13", "cor_tr23")
cat(paste0(paste(fileHead, collapse = ","), "\n"), file = "res_AlphaSimR_multi_trait.csv", append = FALSE)
for (i in 1:rep) {
    set.seed(seeds[i])
    cat("###### This is the", i, "Replication ######\n")
    cat("###### The random seed is", seeds[i], "######\n")
    
    # set Simulation Parameters
    SP <- SimParam$new(founderPop)
    SP$addTraitA(nQtlPerChr = nQTN, mean = rep(0, 3), var = 100 * c(0.2, 0.5, 0.8), corA = corA)
    SP$setVarE(h2 = c(0.2, 0.5, 0.8))
    
    # model the Breeding Program
    pop <- newPop(founderPop)
    
    # get phenotype
    pheno0 <- pheno(pop)
    pheno  <- data.frame(index = 1:nind, T1 = pheno0[, 1], T2 = pheno0[, 2], T3 = pheno0[, 3])
    
    # write phenotype
    write.table(pheno, paste0("AlphaSimR_multi_trait.phe.rep", i), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)

    # genetic evaluation
    completeCmd <- 
      paste("../HIBLUP/hiblup --multi-trait",
        paste("--pheno", paste0("AlphaSimR_multi_trait.phe.rep", i)),
        paste("--pheno-pos", paste(pheno_col, collapse = " ")),
        paste("--bfile", file.path(path, prefix)),
        paste("--add", "--threads 10"),
        paste("--out", paste0("AlphaSimR_multi_trait_rep", i))
      )
    system(completeCmd)
    
    res <- NULL
    vars <- read.table(paste0("AlphaSimR_multi_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[c(1, 3, 5), 4])
    covars <- read.table(paste0("AlphaSimR_multi_trait_rep", i, ".covars"), header = TRUE)
    res <- c(res, covars[1:3, 4])
    
    cat(paste0(paste(res, collapse=","), "\n"), file = "res_AlphaSimR_multi_trait.csv", append = TRUE)

}

res <- read.csv('res_AlphaSimR_multi_trait.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_AlphaSimR_multi_trait.csv', quote = FALSE)
