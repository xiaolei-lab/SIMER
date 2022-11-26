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

pheno_col <- 5
pop.map <- generate.map(pop.marker = nmrk)
qtn.num <- nQTN

# prepare environmental factor list
pop.env <- list(
  F1 = list( # fixed effect 1
    level = c("1", "2"),
    effect = list(tr1 = c(50, 20))
  ), 
  F2 = list( # fixed effect 2
    level = c("d1", "d2", "d3"),
    effect = list(tr1 = c(10, 20, 30))
  ),
  C1 = list( # covariate 1
    level = c(70, 80, 90),
    slope = list(tr1 = 1.5)
  ),
  R1 = list( # random effect 1
    level = paste0("l", 1:100),
    ratio = list(tr1 = 0.1)
  )
)

fileHead <- c("R1", "h2", "pop_mean", "mu", "F1A", "F1B", "F2A", "F2B", "F2C", "C1")
cat(paste0(paste(fileHead, collapse=","), "\n"), file="res_Simer_LMM_trait.csv", append = FALSE)
for (i in 1:rep) {
    set.seed(seeds[i])
    cat("###### This is the", i, "Replication ######\n")
    cat("###### The random seed is", seeds[i], "######\n")
    
    # T1
    SP <- param.annot(pop.map = pop.map, qtn.num = qtn.num, qtn.model = "A")
    SP <- param.geno(SP = SP, pop.geno = geno)
    SP <- param.pheno(
      SP = SP,
      pop.ind = nind,
      pop.env = pop.env,
      phe.model = list(
        tr1 = "T1 = A + F1 + F2 + C1 + R1 + E" # "T1" (Trait 1) consists of Additive effect, F1, F2, C1, R1, and Residual effect
      ),
      phe.var = list(tr1 = 100),
      phe.h2A = list(tr1 = 0.8)
    )
    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)
    
    pheno <- SP$pheno$pop$gen1[, c(1, 8:12), drop = FALSE]
    
    # write phenotype
    write.table(pheno, paste0("Simer_LMM_trait.phe.rep", i), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)
    res <- NULL

    # genetic evaluation for T1
    completeCmd <- 
      paste("../HIBLUP/hiblup --single-trait",
        paste("--pheno", paste0("Simer_LMM_trait.phe.rep", i)),
        paste("--pheno-pos 6"),
        paste("--qcovar 4 --dcovar 2,3 --rand 5"),
        paste("--bfile", file.path(path, prefix)),
        paste("--add", "--threads 10"),
        paste("--out", paste0("Simer_LMM_trait_rep", i))
      )
    system(completeCmd)
    vars <- read.table(paste0("Simer_LMM_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[1:2, 4])
    res <- c(res, mean(pheno[, 6]))
    beta <- read.table(paste0("Simer_LMM_trait_rep", i, ".beta"), header = TRUE)
    res <- c(res, beta[, 2])

    cat(paste0(paste(res, collapse=","), "\n"), file="res_Simer_LMM_trait.csv", append = TRUE)
    
}

res <- read.csv('res_Simer_LMM_trait.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_Simer_LMM_trait.csv', quote = FALSE)
