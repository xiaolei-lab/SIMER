library(simer)

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
pop.map <- generate.map(pop.marker = nmrk)
qtn.num <- list(tr1 = nQTN, tr2 = nQTN, tr3 = nQTN)
phe.corA <- matrix(c(1.0, 0.2, 0.5,
                     0.2, 1.0, 0.8,
                     0.5, 0.8, 1.0), 3, 3, byrow = TRUE)

fileHead <- c("h2_tr1", "h2_tr2", "h2_tr3", "cor_tr12", "cor_tr13", "cor_tr23")
cat(paste0(paste(fileHead, collapse=","), "\n"), file = "res_Simer_multi_trait.csv", append = FALSE)
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
      phe.model = list(
        tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
        tr2 = "T2 = A + E", # "T2" (Trait 2) consists of Additive effect and Residual effect
        tr3 = "T3 = A + E"  # "T3" (Trait 3) consists of Additive effect and Residual effect
      ),
      phe.var = list(tr1 = 100, tr2 = 100, tr3 = 100),
      phe.h2A = list(tr1 = 0.2, tr2 = 0.5, tr3 = 0.8),
      phe.corA = phe.corA
    )

    SP <- annotation(SP)
    SP <- genotype(SP)
    SP <- phenotype(SP)

    # get phenotype
    pheno0 <- SP$pheno$pop$gen1[, 8:10]
    pheno  <- data.frame(index = 1:nind, T1 = pheno0[, 1], T2 = pheno0[, 2], T3 = pheno0[, 3])
    
    # write phenotype
    write.table(pheno, paste0("Simer_multi_trait.phe.rep", i), sep = '\t', quote = FALSE, row.names = FALSE, col.name = TRUE)

    # genetic evaluation
    completeCmd <- 
      paste("../HIBLUP/hiblup --multi-trait",
        paste("--pheno", paste0("Simer_multi_trait.phe.rep", i)),
        paste("--pheno-pos", paste(pheno_col, collapse = " ")),
        paste("--bfile", file.path(path, prefix)),
        paste("--add", "--threads 10"),
        paste("--out", paste0("Simer_multi_trait_rep", i))
      )
    system(completeCmd)
    
    res <- NULL
    vars <- read.table(paste0("Simer_multi_trait_rep", i, ".vars"), header = TRUE)
    res <- c(res, vars[c(1, 3, 5), 4])
    covars <- read.table(paste0("Simer_multi_trait_rep", i, ".covars"), header = TRUE)
    res <- c(res, covars[1:3, 4])
    
    cat(paste0(paste(res, collapse=","), "\n"), file = "res_Simer_multi_trait.csv", append = TRUE)

}

res <- read.csv('res_Simer_multi_trait.csv')
Mean <- round(colMeans(res), digits = 4)
S.D. <- round(sqrt(diag(var(res))) / 10, digits = 4)
final_res <- rbind(Mean, S.D.)
write.csv(final_res, 'final_res_Simer_multi_trait.csv', quote = FALSE)

