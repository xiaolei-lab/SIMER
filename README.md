# SIMER [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/XiaoleiLiuBio/rMVP/issues/new) 

## Data Simulation for Life Science and Breeding

### Authors:
Designed and Maintained by [Dong Yin](https://github.com/Foredawnbio), [Xuanning Zhang](https://github.com/zxnyaonuli), [Lilin Yin](https://github.com/YinLiLin), [Haohao Zhang](https://github.com/hyacz), and [**Xiaolei Liu**](https://github.com/XiaoleiLiuBio).  
Contributors: Zhenshuang Tang, Jingya Xu, Xinyun Li, Mengjin Zhu, Xiaohui Yuan, and Shuhong Zhao

Questions, suggestions, and bug reports are welcome and appreciated: [xiaoleiliu@mail.hzau.edu.cn](mailto:xiaoleiliu@mail.hzau.edu.cn)

### Contents
<!-- TOC updateOnSave:false -->

- [Installation](#installation)
- [Data Preparation](#data-preparation)
    - [Genotype](#genotype)
    - [Genotypic map](#genotypic-map)
    - [Pedigree](#pedigree)<img src="https://raw.githubusercontent.com/xiaolei-lab/SIMER/master/inst/extdata/00simer_logo/simer_logo.png" height="250" align="right" />
- [Data Input](#data-input)
    - [Basic](#basic)
    - [Optional](#optional)
- [Quick Start](#quick-start)
- [Genotype Simulation](#genotype-simulation)
    - [Gallery of genotype simulation input parameters](#gallery-of-genotype-simulation-input-parameters)
    - [Generate genotype matrix of base population](#generate-genotype-matrix-of-base-population)
    - [Set block information and recombination information](#set-block-information-and-recombination-information)
    - [Add chromosome crossovers and mutaions to genotype matrix](#add-chromosome-crossovers-and-mutaions-to-genotype-matrix)
- [Phenotype Simulation](#phenotype-simulation)  
    - [Gallery of phenotype simulation input parameters](#gallery-of-phenotype-simulation-input-parameters)  
    - [Generate base population information](#generate-base-population-information)  
    - [Generate phenotype of single trait by A model](#generate-phenotype-of-single-trait-by-A-model)  
    - [Generate phenotype of single trait by AD model](#generate-phenotype-of-single-trait-by-AD-model)  
    - [Generate phenotype of single trait by ADI model](#generate-phenotype-of-single-trait-by-ADI-model)  
    - [Generate phenotype of multiple traits](#generate-phenotype-of-multiple-traits)  
    - [Different QTN effect distributions](#different-QTN-effect-distributions)  
    - [Different selection criteria](#different-selection-criteria)  
    - [Multiple groups QTN effects](#multiple-groups-QTN-effects)  
- [Selection](#selection)
    - [Gallery of selection input parameters](#gallery-of-selection-input-parameters)
    - [Individual selection on single trait](#individual-selection-on-single-trait)
    - [Family selection on single trait](#family-selection-on-single-trait)
    - [Within family selection on single trait](#within-family-selection-on-single-trait)
    - [Combined selection on single trait](#combined-selection-on-single-trait)
    - [Tandem selection on multiple traits](#tandem-selection-on-multiple-traits)
    - [Independent culling selection on multiple traits](#independent-culling-selection-on-multiple-traits)
    - [Index selection on multiple traits](#index-selection-on-multiple-traits)
- [Reproduction](#reproduction) 
    - [Gallery of reproduction input parameters](#gallery-of-reproduction-input-parameters)
    - [Clone](#clone)
    - [Double haploid](#double-haploid)
    - [Self pollination](#self-pollination)
    - [Random mating](#random-mating)
    - [Random mating without self pollination](#random-mating-without-self-pollination)
    - [User designed pedigree mating](#user-designed-pedigree-mating)
    - [Two way cross](#two-way-cross)
    - [Three way cross](#three-way-cross)
    - [Four way cross](#four-way-cross)
    - [Back cross](#back-cross)
- [Output](#output)
    - [Population information](#population-information)
    - [Marker effects](#marker-effects)
    - [Trait information](#trait-information)
    - [Genotype](#genotype)
    - [Genotypic id](#genotypic-id)
    - [Genotypic map](#genotypic-map)
    - [Selection intensity](#selection-intensity)
- [Citation](#citation)
- [FAQ and Hints](#faq-and-hints)

<!-- /TOC -->

---

# Installation
**[back to top](#contents)**  

**WE STRONGLY RECOMMEND TO INSTALL SIMER ON Microsoft R Open(https://mran.microsoft.com/download/)**  

## Installation

``` r
# please wait
install.packages("simer")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("xiaolei-lab/SIMER")
```

After installed successfully, **SIMER** can be loaded by typing
```r
> library(simer)
```
Typing ```?simer``` could get the details of all parameters.

---

# Data Preparation

## Genotype
**[back to top](#contents)**  

Genotype data should be  Numeric (m rows and (2 * n) columns, m is the number of SNPs, n is the number of individuals) format. If you have genotype data in **PLINK Binary** format (details see http://zzz.bwh.harvard.edu/plink/data.shtml#bed), **VCF** or **Hapmap**, please convert them using "MVP.Data" function in the **rMVP**(https://github.com/xiaolei-lab/rMVP).

> `genotype.txt`

<table>
<tbody>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr></tbody></table>

## Genotypic map
**[back to top](#contents)**  
Map file is necessary in **SIMER**. It will generate genotype matrix according to number of markers in input map file. 

> `map.txt`

| SNP | Chrom | BP | REF | ALT |
| :---: | :---: |:---: |:---: |:---: |
| 1_10673082 | 1 | 10673082 | T | C |
| 1_10723065 | 1 | 10723065 | A | G |
| 1_11407894 | 1 | 11407894 | A | G |
| 1_11426075 | 1 | 11426075 | T | C |
| 1_13996200 | 1 | 13996200 | T | C |
| 1_14638936 | 1 | 14638936 | T | C |

## Pedigree
**[back to top](#contents)**  
**SIMER** supports user designed pedigree to control mating process. User designed pedigree is useful only in "userped" reproduction. Pedigree should at least start with generation 2. The first column is sample id, the sescond column is paternal id, and the third column is maternal id. Please make sure that paternal id and maternal id can be found in the last generation. 

> `userped.txt`

| Index | Sire | Dam |
| :---: | :---: |:---: |
| 41 | 1 | 11 |
| 42 | 1 | 11 |
| 43 | 1 | 11 |
| 44 | 1 | 11 |
| 45 | 2 | 12 |
| 46 | 2 | 12 |

---

# Data Input

## Basic
**[back to top](#contents)**  
At least you should prepare two datasets: genotypic map and genotype.  

**genotype**, genotype data in **Numeric** format (m * (2 * n), m rows and n columns, m is the number of SNPs, n is the number of individuals)
**genotypic map**, SNP map information, the first column is SNP name, the second column is Chromosome ID, the third column is phsical position, the fourth column is REF, and the fifth column is ALT  

```r
rawgeno <- read.table("genotype.txt")
input.map <- read.table("map.txt" , head = TRUE)
```

## Optional
**[back to top](#contents)**  
If you want to control mating process by user designed pedigree. 

**pedigree**, pedigree information, the first column is sample id, the second column is paternal id, and the third column is maternal id. Note that the individuals in the pedigree data file do not need to be sorted by the date of birth, and the missing value can be replaced by NA or 0.

```r
userped <- read.table("userped.txt", header = TRUE)
```

---

# Quick Start
**[back to top](#contents)** 

After obtaining genotypic map data and genotype data, we can start our simulation.

**num.gen**, number of generations in simulation  
**replication**, replication index of simulation  
**verbose**, whether to print detail  
**mrk.dense**, whether markers are dense  
**out path**, path of output files  
**out.format**, format of output, "numeric" or "plink"  
**seed.geno**, random seed of genotype matrix  
**seed.map**, random seed of map file  
**out.geno.gen**, indice of generation of output genotype   
**out.pheno.gen**, indice of generation of  output phenotype   
**rawgeno1**, extrinsic genotype matrix1  
**rawgeno2**, extrinsic genotype matrix2    
**rawgeno3**, extrinsic genotype matrix3  
**rawgeno4**, extrinsic genotype matrix4  
**num.ind**, population size of base population  
**prob**, weight of "0" and "1" in genotype matrix, the sum of element in vector equals 1  
**input.map**, map from outside  
**len.block**, length of every blocks  
**range.hot**, range of exchages in hot spot block  
**range.cold**, range of exchages in cold spot block  
**rate.mut**, mutation rate between 1e-8 and 1e-6  
**cal.model**, phenotype model with "A", "AD", "ADI"  
**h2.tr1**, heritability vector of trait1, corresponding to a, d, aXa, aXd, dXa, dXd  
**num.qtn.tr1**, integer or integer vector, the number of QTN in the trait1  
**var.tr1**, variances of different effects, the last 5 vector elements are corrresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a  
**dist.qtn.tr1**, distribution of QTN's effects with options: "normal", "geometry" and "gamma", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**eff.unit.tr1**, unit effect of geometric distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**shape.tr1**, shape of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**scale.tr1**, scale of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**multrait**, whether applying pair traits with overlapping, TRUE represents applying, FALSE represents not  
**num.qtn.trn**, QTN distribution matrix, diagnal elements are total QTN number of the trait, non-diagnal are QTN number of overlop qtn  
**eff.sd**, a matrix with the standard deviation of QTN effects  
**gnt.cov**, genetic covaiance matrix among all traits  
**env.cov**, environment covaiance matrix among all traits  
**qtn.spot**, QTN probability in every blocks  
**maf**, Minor Allele Frequency, markers selection range is from  maf to 0.5  
**sel.crit**, selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"  
**sel.on**, whether to add selection  
**mtd.reprod**, different reproduction methods with options: "clone", "dh", "selfpol", "singcro", "tricro", "doubcro", "backcro","randmate", "randexself" and "userped"  
**userped**, user-specific pedigree  
**num.prog**, litter size of dams  
**ratio**, ratio of males in all individuals  
**prog.tri**, litter size of the first single cross process in trible cross process  
**prog.doub**, litter size of the first two single cross process in double cross process  
**prog.back**, a vector with litter size in every generations  
**ps**, fraction selected in selection  
**decr**, whether to sorting with descreasing  
**sel.multi**, selection method of multi-trait with options: "tdm", "indcul" and "index"  
**index.wt**, economic weights of selection index method  
**index.tdm**, index represents which trait is being selected. NOT CONTROL BY USER  
**goal.perc**, percentage of goal more than mean of scores of individuals  
**pass.perc**, percentage of expected excellent individuals  
**sel.sing**, selection method of single trait with options: "ind", "fam", "infam" and "comb"  

```r
# simdata: rawgeno input.map userped
data(simdata)

simer.list <-
    simer(num.gen = 5,
          verbose = TRUE, 
          out = NULL,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          # num.ind = 100,
          mtd.reprod = "randmate",
          num.prog = 2,
          ratio = 0.5)
```

If you want to make replicated simulation

```r
rep <- 2
for (i in 1:rep) {
  simer.list <-
        simer(num.gen = 3,
              replication = i, # set index of replication
              verbose = TRUE, 
              out = NULL,
              input.map = input.map,
              rawgeno1 = rawgeno, # use your own genotype matrix
              # num.ind = 100,
              mtd.reprod = "randmate",
              num.prog = 2,
              ratio = 0.5)
}
```

---

# Genotype Simulation
**[back to top](#contents)** 

Genotype data in SIMER will be generated randomly or from outside genotype matrix. Chromosome crossovers and base mutations depend on block information and recombination informaion of map. 

## Gallery of genotype simulation input parameters
**[back to top](#contents)** 

`genotype()`, main function of genotype simulation:  
**rawgeno**, extrinsic genotype matrix  
**geno**, genotype matrix need dealing with  
**num.marker**, number of markers  
**num.ind**, population size of base population  
**prob**, weight of "0" and "1" in genotype matrix, the sum of element in vector equals 1  
**blk.rg**, it represents the started and ended position blocks  
**recom.spot**, whether to consider recombination in every blocks  
**range.hot**, range of exchages in hot spot block  
**range.cold**, range of exchages in cold spot block  
**rate.mut**, mutation rate between 1e-8 and 1e-6  
**verbose**, whether to print details  

 `check.map()`, add block id and combination information to genotypic map:   
**input.map**, map from outside  
**num.marker**, number of markers   
**len.block**, length of every blocks  

`cal.blk()`, get start position and end position of blocks:  
**pos.map**, map with block information and recombination information  

`input.geno()`, input sub-genotype matrix to total genotype matrix:  
**bigmtr**, total genotype matrix  
**mtr**, genotype matrix should be inputting  
**ed**, index of the last column in each process  
**mrk.dense**, whether markers are dense 

`simer()`, main function:  
**rawgeno1**, extrinsic genotype matrix1  
**rawgeno2**, extrinsic genotype matrix2    
**rawgeno3**, extrinsic genotype matrix3  
**rawgeno4**, extrinsic genotype matrix4 
 
## Generate genotype matrix of base population
**[back to top](#contents)** 

There are two different ways to generate genotype matrix of base population.

```r
# use num.marker and num.ind to generate a new genotype matrix
nmrk <- nrow(input.map)
basepop.geno <- genotype(num.marker = nmrk, num.ind = 40, verbose = verbose)

# use genotype matrix from outside
basepop.geno <- genotype(rawgeno = rawgeno, verbose = verbose)
```

## Set block information and recombination information
**[back to top](#contents)** 

Add block information and recombination information to map file. Calculate block ranges and recombination states of blocks.

```r
nmrk <- nrow(basepop.geno)
nind <- ncol(basepop.geno) / 2

# set block information and recombination information
pos.map <- check.map(input.map = input.map, num.marker = nmrk, len.block = 5e7)
blk.rg <- cal.blk(pos.map)
recom.spot <- as.numeric(pos.map[blk.rg[, 1], 7])
```

## Add chromosome crossovers and mutaions to genotype matrix
**[back to top](#contents)** 

After getting block information and recombination information, you can add chromosome crossovers and mutations to genotype matrix.  

```r
basepop.geno.em <-  # genotype matrix after cross and Mutation
    genotype(geno = basepop.geno,
             blk.rg = blk.rg,
             recom.spot = recom.spot,
             range.hot = 4:6,
             range.cold = 1:5,
             rate.mut = 1e-8, 
             verbose = verbose)
```

Note that recombination only exists in meiosis. Therefore, some reproduction methods such like "clone" do not have recombination process. You can set **recom.spot = NULL** to add only mutations to genotype matrix.  

```r
basepop.geno.em <-  # genotype matrix after crosses and mutations
    genotype(geno = basepop.geno,
             blk.rg = blk.rg,
             recom.spot = NULL, # only mutations on genotype matrix
             range.hot = 4:6,
             range.cold = 1:5,
             rate.mut = 1e-8, 
             verbose = verbose)
```             

---

# Phenotype Simulation
**[back to top](#contents)**  

Phenotype data in **SIMER** will be generated according to different phenotype model, QTN effect distribution and selection criteria. SIMER supports both single trait and multiple trait. In single traits, you could set different model by **cal.model** as your need. But in multiple traits, only "A" model can be used. In both single trait and multiple traits, you can set different amount, effect variance, effect distribution for QTNs and heritability for traits. In single trait, multiple groups QTN effects can also be under consideration. In multiple trait, you can set specific genetic correlation and environmental correlation. 

## Gallery of phenotype simulation input parameters
**[back to top](#contents)**  

`phenotype()`, main function of phenotype simulation:  
**effs**, a list with number of overlap markers, selected markers, effects of markers  
**pop**, population information of generation, family index, within-family index, index, sire, dam, sex  
**pop.geno**, genotype matrix of population, two columns represent a individual  
**pos.map**, marker information of population  
**h2.tr1**, heritability vector of trait1, corresponding to a, d, aXa, aXd, dXa, dXd  
**gnt.cov**, genetic covaiance matrix among all traits  
**env.cov**, environment covaiance matrix among all traits  
**sel.crit**, selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"  
**pop.total**, total population infarmation  
**sel.on**, whether to add selection
**inner.env**, environment of main function of simer  
**verbose**, whether to print details  

`cal.effs()`, calculate for marker effects: 
**pop.geno**, genotype matrix of population, two columns represent a individual  
**cal.model**, phenotype model with "A", "AD", "ADI"  
**num.qtn.tr1**, integer or integer vector, the number of QTN in the trait1  
**var.tr1**, variances of different effects, the last 5 vector elements are corrresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a  
**dist.qtn.tr1**, distribution of QTN's effects with options: "normal", "geometry" and "gamma", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**eff.unit.tr1**, unit effect of geometric distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**shape.tr1**, shape of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**scale.tr1**, scale of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**multrait**, whether applying pair traits with overlapping, TRUE represents applying, FALSE represents not  
**num.qtn.trn**, QTN distribution matrix, diagnal elements are total QTN number of the trait, non-diagnal are QTN number of overlop qtn  
**eff.sd**, a matrix with the standard deviation of QTN effects  
**qtn.spot**, QTN probability in every blocks  
**maf**, Minor Allele Frequency, markers selection range is from  maf to 0.5  
**verbose**, whether to print details  

`getpop()`, generate population information:  
**nind**, number of individuals in the population  
**from**, initial index of the population  
**ratio**, ratio of males in all individuals  

`set.pheno()`, add phenotype to population information:  
**pop**, population information of generation, family index, within-family index, index, sire, dam, sex  
**pop.pheno**, phenotype information
**sel.crit**, selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"  

## Generate base population information
**[back to top](#contents)** 

Generate base population information according to population size and sex ratio.  

```r
# base population for single trait
basepop1 <- getpop(nind = nind, from = 1, ratio = 0.1)

# base population for double traits
basepop2 <- getpop(nind = nind, from = nind + 1, ratio = 0.1)
```

## Generate phenotype of single trait by A model
**[back to top](#contents)** 

In "A" model, **SIMER** considers only Additive effect(a). Therefore, only the first elements of **var.tr1**, **dist.qtn.tr1**, **eff.unit.tr1**, **shape.tr1** , and **scale.tr1** are useful for "A" model. Add phenotypes of single trait to base population1 are displayed as follows: 

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A",
             num.qtn.tr1 = 18,
             var.tr1 = 2,
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop1.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop1.pheno, sel.crit = "pheno")
```

## Generate phenotype of single trait by AD model
**[back to top](#contents)**

In "AD" model, **SIMER** considers both Additive effect(a) and Dominance effect(d). Therefore, only the first two elements of **var.tr1**, **dist.qtn.tr1**, **eff.unit.tr1**, **shape.tr1**,  and **scale.tr1** are useful for "AD" model. Add phenotypes of single trait to base population1 are displayed as follows:

```r
####################
### single trait ###
# calculate for marker information
# Additive by Dominance model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "AD", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = c(0.4, 0.2), # (a, d)
             dist.qtn.tr1 = c("normal", "normal"),
             eff.unit.tr1 = c(0.5, 0.5),
             shape.tr1 = c(1, 1),
             scale.tr1 = c(1, 1),
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.5, 0.3),
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = TRUE,
               inner.env = NULL,
               verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

Note that real variance components ratio may not be consistent with expected herirability **h2.tr1**. You can set **sel.on = FALSE** and pass a **inner.env** to get a more accurate one. In addition, markers effects will be corrected in this case.

```r
####################
### single trait ###
# calculate for marker information
# Additive by Dominance model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "AD", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = c(0.4, 0.2), # (a, d)
             dist.qtn.tr1 = c("normal", "normal"),
             eff.unit.tr1 = c(0.5, 0.5),
             shape.tr1 = c(1, 1),
             scale.tr1 = c(1, 1),
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
inner.env <- environment()
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.5, 0.3),
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = FALSE,
               inner.env = inner.env,
               verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

## Generate phenotype of single trait by ADI model
**[back to top](#contents)**

In "ADI" model, **SIMER** considers not only Additive effect(a) and Dominance effect(d) but also interactive effects including Additive by Additive effect(aXa), Additive by Dominance effect(aXd), Dominance by Additive effect(dXa) and Dominance by Dominance effect(dXd). Therefore, all six elements of **var.tr1**, **dist.qtn.tr1**, **eff.unit.tr1**, **shape.tr1**, and **scale.tr1** are useful for "ADI" model. Meanwhile, QTN amount should be an even in "ADI" model. Add phenotypes of single trait to base population1 are displayed as follows:

```r
####################
### single trait ###
# calculate for marker information
# Additive by Dominance by Iteratiion model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "ADI", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(18),
             var.tr1 = c(2, 1, 0.5, 0.5, 0.5, 0.1),
             dist.qtn.tr1 = rep("normal", times = 6),
             eff.unit.tr1 = rep(0.5, 6), # (a, d, aXa, aXd, dXa, dXd)
             shape.tr1 = rep(1, times = 6),
             scale.tr1 =rep(1, times = 6),
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05),
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = TRUE,
               inner.env = NULL,
               verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

Note that real variance components ratio may not be consistent with expected herirability **h2.tr1**. You can set **sel.on = FALSE** and pass a **inner.env** to get a more accurate one. In addition, markers effects will be corrected in this case.

```r
####################
### single trait ###
# calculate for marker information
# Additive by Dominance model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "ADI", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(18),
             var.tr1 = c(2, 1, 0.5, 0.5, 0.5, 0.1),
             dist.qtn.tr1 = rep("normal", times = 6),
             eff.unit.tr1 = rep(0.5, 6), # (a, d, aXa, aXd, dXa, dXd)
             shape.tr1 = rep(1, times = 6),
             scale.tr1 =rep(1, times = 6),
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
inner.env <- environment()
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05),
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = FALSE,
               inner.env = inner.env,
               verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

## Generate phenotype of multiple traits
**[back to top](#contents)**  

In multiple traits, only "A" model is applied and **multrait** should be TRUE. If you want to generate multiple traits with specific genetic correlation and environmental correlation please assign genetic covariance matrix **gnt.cov** and environmental covariance matrix **env.cov**. But when **sel.on** is TRUE, these traits will have random genetic correlation. By designed gnt.cov and env.cov, heritability will be generated indirectly. For example, heritability of the first trait equals to gnt.cov[1, 1] / (gnt.cov[1, 1] + env.cov[1, 1]). Add phenotype of multiple  traits to base population2 are displayed as follows:

```r
#######################
### multiple traits ###
# calculate for marker information
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(2, 6, 10),
             var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
             dist.qtn.tr1 = rep("normal", 6),
             eff.unit.tr1 = rep(0.5, 6),
             shape.tr1 = rep(1, 6),
             scale.tr1 = rep(1, 6),
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = matrix(c(1, 0, 0, 2), 2, 2),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop2.pheno <-
     phenotype(effs = effs,
               pop = basepop2,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05),
               gnt.cov = matrix(c(14, 10, 10, 15), 2, 2),
               env.cov = matrix(c(6, 5, 5, 10), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop2,
               sel.on = TRUE, 
               inner.env = NULL,
               verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop2.pheno, sel.crit = "pheno")
```

Note that real genetic covariance matrix may not be consistent with expected genetic covariance matrix **gnt.cov**. You can set **sel.on = FALSE** and pass a **inner.env** to get a more accurate one. In addition, markers effects will be corrected in this case.

```r
#######################
### multiple traits ###
# calculate for marker information
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(2, 6, 10),
             var.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
             dist.qtn.tr1 = rep("normal", 6),
             eff.unit.tr1 = rep(0.5, 6),
             shape.tr1 = rep(1, 6),
             scale.tr1 = rep(1, 6),
             multrait = TRUE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = matrix(c(1, 0, 0, 2), 2, 2),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
inner.env <- environment()
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop2,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05),
               gnt.cov = matrix(c(14, 10, 10, 15), 2, 2),
               env.cov = matrix(c(6, 5, 5, 10), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop2,
               sel.on = FALSE, 
               inner.env = inner.env,
               verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "pheno")
```

## Different QTN effect distributions
**[back to top](#contents)**  

In different model, you can further set different QTN effect distributions of trait1 by \verb|dist.qtn.tr1|. The most common distribution is "normal" distribution. You can set different variances in "normal" distribution by **dist.qtn.tr1**.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

QTN effect distribution can be "geometry" distribution. You can set effect unit of "geometry" by **eff.unit.tr1**.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # a
             dist.qtn.tr1 = "geometry",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1, # effect unit of geomtry distribution
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

"Gamma" distribution is also a kind of QTN effect distribution. You can set shape and scale of "gamma" distribution by **shape.tr1** and **scale.tr1**. Note that default options of "gamma" distribution **shape.tr1 = 1** and **scale.tr1 = 1** exactly lead to exponential distribution.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # a
             dist.qtn.tr1 = "gamma",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1, # shape of gamma distribution
             scale.tr1 = 1, # scale of gamma distribution
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

## Different selection criteria
**[back to top](#contents)**  

In addition, "pheno" is not just phenotype. It can also be "TBV", "TGV", "pEBVs'", "gEBVs",  or "ssEBVs". "TBV" is True Breeding Value, represents only that part of genotypic value that can be transmitted from parent to offspring. 

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TBV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "TBV",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "TBV")
```

Phenotype of multiple traits can also be represented as "TBV".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TBV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "TBV",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "TBV")
```

"TGV" is True Genotypic Value, represents the sum of additive effect, dominance effect and epistatic effect.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TGV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "TGV",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "TGV")
```

Phenotype of multiple traits can also be represented as "TGV".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TGV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "TGV",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "TGV")
```

"pEBVs" is pedigree Estimated Breeding Values. It means that BLUP constructs kinship by pedigree. You can get "pEBVs" by "ABLUP" model in HIBLUP. 

```r
# call "hiblup" package
suppressMessages(library("hiblup"))

####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pEBVs"
# "pEBVs" needs hiblup package
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pEBVs")
```

Phenotype of multiple traits can also be represented as "pEBVs".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "pEBVs")
```

"gEBVs" is genomic Estimated Breeding Values. It means that BLUP constructs kinship by genotype matrix. You can get "gEBVs" by "GBLUP" model in HIBLUP.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "gEBVs"
# "gEBVs" needs hiblup package
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "gEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "gEBVs")
```

Phenotype of multiple traits can also be represented as "gEBVs".

```r
# Additive model
#######################
### multiple traits ###
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "gEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "gEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "gEBVs")
```

"ssEBVs" is single-step genomic Estimated Breeding Values. It means that BLUP constructs kinship by pedigree and genotype matrix. You can get "ssEBVs" by "SSBLUP" model in HIBLUP.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "ssEBVs"
# "ssEBVs" needs hiblup package
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "ssEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "ssEBVs")
```

Phenotype of multiple traits can also be represented as "ssEBVs".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "ssEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "ssEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "ssEBVs")
```

At last, "pheno" is phenotype including additive effect (and dominance effect) (and epistatic effect) and residual effect.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pheno"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

Phenotype of multiple traits can also be represented as "pheno".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             var.tr1 = 0.6, # variance of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pheno"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop,
              pop.geno = basepop.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              env.cov = matrix(c(10, 5, 5, 100), 2, 2),
              sel.crit = "pheno",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# add phenotype to basepop
basepop2 <- set.pheno(pop = basepop2, pop.pheno, sel.crit = "pheno")
```

## Multiple groups QTN effects
**[back to top](#contents)**  

Multiple groups QTN effects can be realized by setting different elements of **num.qtn.tr1**, every elements represent amount of QTNs affacting a effect. For example, **num.qtn.tr1 = c(2, 6, 10)** means that the additive effect of the trait is the sum of three QTN group effects. The first group has 2 QTNs, the second group has 6 QTNs and the third group has 10 QTNs. Because of the three groups, the first three elements of **var.tr1** mean the variances of the three QTN groups.

```r
#####################
### single traits ###
# calculate for marker information
effs <-
    cal.effs(pop.geno = basepop.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(2, 6, 10),
             var.tr1 = c(0.4, 0.2, 0.02),
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 1,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             eff.sd = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop.geno,
               pos.map = NULL,
               h2.tr1 = 0.3,
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               env.cov = matrix(c(10, 5, 5, 100), 2, 2),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = TRUE,
               inner.env = NULL,
               verbose = verbose)

# add phenotype to basepop
basepop1 <- set.pheno(pop = basepop1, pop.pheno, sel.crit = "pheno")
```

---

# Selection
**[back to top](#contents)**  

You can get ordered individuals indice according to phenotype in the populaton information. Fraction selected can be used to keep a certain amount of individuals. SIMER chooses automatically single trait selection or multiple traits selection according to number of columns of phenotype.

# Gallery of selection input parameters
**[back to top](#contents)**  

`selects()`, main function of selection:  
**pop**, population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye  
**decr**, whether to sorting with descreasing  
**sel.multi**, selection method of multi-trait with options: "tdm", "indcul" and "index"  
**index.wt**, economic weights of selection index method  
**index.tdm**, index represents which trait is being selected. NOT CONTROL BY USER  
**goal.perc**, percentage of goal more than mean of scores of individuals  
**pass.perc**, percentage of expected excellent individuals  
**sel.sing**, selection method of single trait with options: "ind", "fam", "infam" and "comb"  
**pop.total**, total population infarmation  
**pop.pheno**, list of all phenotype information  
**verbose**, whether to print detail  

`simer()`, main function:  
**ps**, fraction selected in selection  

# Individual selection on single trait
**[back to top](#contents)**  

Individual selection is a method of selecting according to the phenotype of individual traits, also known as mixed selection or collective selection. This selection method is simple and easy to used for traits with high heritability.

```r
# output index.tdm and ordered individuals indice
# single trait selection
# individual selection
ind.ordered <-
    selects(pop = basepop1, # population with single trait
            decr = TRUE, # sort individuals by descreasing
            sel.sing = "ind",
            pop.total = basepop1,
            pop.pheno = pop1.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Family selection on single trait
**[back to top](#contents)** 

Family selection is a method of selecting by family based on the average of the family. This selection method is used for traits with low heritability.

```r
# output index.tdm and ordered individuals indice
# single trait selection
# family selection
ind.ordered <-
    selects(pop = basepop1, # population with single trait
            decr = TRUE, # sort individuals by descreasing
            sel.sing = "fam",
            pop.total = basepop1,
            pop.pheno = pop1.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Within family selection on single trait
**[back to top](#contents)** 

Within-family selection is a method of selecting according to the deviation of individual phenotype and family mean value in each family. This selection method is used for traits with low heritability and small family.

```r
# output index.tdm and ordered individuals indice
# single trait selection
# within-family selection
ind.ordered <-
    selects(pop = basepop1, # population with single trait
            decr = TRUE, # sort individuals by descreasing
            sel.sing = "infam",
            pop.total = basepop1,
            pop.pheno = pop1.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Combined selection on single trait
**[back to top](#contents)**  

Combined selection is a method of selecting according to weighed combination of the deviation of individual phenotype and family mean value  and family mean value.

```r
# output index.tdm and ordered individuals indice
# single trait selection
# combined selection
ind.ordered <-
    selects(pop = basepop1, # population with single trait
            decr = TRUE, # sort individuals by descreasing
            sel.sing = "comb",
            pop.total = basepop1,
            pop.pheno = pop1.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Tandem selection on multiple traits
**[back to top](#contents)**  

Tandem selection is a method for sequentially selecting a plurality of target traits one by one. The index of selected trait is **index.tdm** and this parameter should not controlled by User.

```r
# output index.tdm and ordered individuals indice
# multiple traits selection
# tandem selection
ind.ordered <-
    selects(pop = basepop2, # population with multiple traits
            decr = TRUE,
            sel.multi = "tdm",
            index.wt = c(0.5, 0.5),
            index.tdm = 1,
            goal.perc = 0.1,
            pass.perc = 0.9,
            pop.total = basepop2,
            pop.pheno = pop2.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Independent culling selection on multiple traits
**[back to top](#contents)**  

After setting a minimum selection criterion for each target trait. Independent culling selection will eliminate this individual when the candidate's performance on any trait is lower than the corresponding criteria.

```r
# output index.tdm and ordered individuals indice
# multiple traits selection
# Independent culling selection
ind.ordered <-
    selects(pop = basepop2, # population with multiple traits
            decr = TRUE,
            sel.multi = "indcul",
            index.wt = c(0.5, 0.5),
            index.tdm = 1,
            goal.perc = 0.1,
            pass.perc = 0.9,
            pop.total = basepop2,
            pop.pheno = pop2.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

# Index selection on multiple traits
**[back to top](#contents)**  

Index selection is a comprehensive selection that will consider several traits based on their respective heritabilities, phenotypic variances, economic weights, and corresponding genetic correlations and phenotypes. Then calculate the index value of each body, and eliminate or select it according to its level. You can set the weight of Index selection by **index.wt**.

```r
# output index.tdm and ordered individuals indice
# multiple traits selection
# index selection
ind.ordered <-
    selects(pop = basepop2, # population with multiple traits
            decr = TRUE,
            sel.multi = "index",
            index.wt = c(0.5, 0.5), # trait1 and trait2
            index.tdm = 1,
            goal.perc = 0.1,
            pass.perc = 0.9,
            pop.total = basepop2,
            pop.pheno = pop2.pheno, 
            verbose = verbose)
index.tdm <- ind.ordered[1]
ind.ordered <- ind.ordered[-1]
ind.ordered
```

---

# Reproduction
**[back to top](#contents)**  

Different kind of reproduction methods need different preparation. Reproduction methods in **mtd.reprod** can be devided into single-breeds and multi-breeds methods. In every reproduction methods, you can set litter size **num.prog** and sex ratio **ratio** in freedom.

# Gallery of reproduction input parameters
**[back to top](#contents)** 

`simer()`, main function:  
**mtd.reprod**, different reproduction methods with options: "clone", "dh", "selfpol", "singcro", "tricro", "doubcro", "backcro","randmate", "randexself" and "userped"  
**userped**, user-specific pedigree  
**num.prog**, litter size of dams  
**ratio**, ratio of males in all individuals  
**prog.tri**, litter size of the first single cross process in trible cross process  
**prog.doub**, litter size of the first two single cross process in double cross process  
**prog.back**, a vector with litter size in every generations  

# Clone
**[back to top](#contents)** 

Asexual reproduction does not involve germ cells, and does not require a process of fertilization, directly forming a new individual's reproductive mode from a part of the mother. Sex of offspring will be "0" in clone. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# clone
simer.list <-
    simer(num.gen = 3,
          verbose = verbose, 
          out = out, 
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "clone",
          num.prog = 2,
          ratio = 0.5)
```

# Double haploid
**[back to top](#contents)**  

Breeding workers often use another culture in vitro to obtain haploid plants, and then artificially induced to double the number of chromosomes and restore the number of chromosomes in normal plants. This method is named double-haploid reproduction. Sex of offspring will be "0" in double-haploid. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# double-haploid
simer.list <-
    simer(num.gen = 4,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "dh",
          num.prog = 2,
          ratio = 0.5)
```

# Self pollination
**[back to top](#contents)** 

Self-pollination refers to the combination of male and female gametes from the same individual or between individuals from the same clonal breeding line. Sex of offspring will be "0" in self-pollination. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.  

```r
# self-pollination
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "selfpol",
          num.prog = 2,
          ratio = 0.5)
```

# Random mating
**[back to top](#contents)**  

In random mating, any female or male individual have the same probability to mate with any opposite sex in a sexually reproducing organism. Sex of offspring in random mating is up to sex of parents. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# random-mating
simer.list <-
    simer(num.gen = 4,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 100,
          mtd.reprod = "randmate",
          num.prog = 2,
          ratio = 0.5)
```

# Random mating without self pollination
**[back to top](#contents)**  

In random mating without self-pollination, a individual cannot mate to itself. Sex of offspring in random mating without self-pollination is up to sex of parents. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# random-mating without self-pollination
simer.list <-
    simer(num.gen = 3,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 100,
          mtd.reprod = "randexself",
          num.prog = 2,
          ratio = 0.5)
```

# User designed pedigree mating
**[back to top](#contents)**  

User-designed-pedigree mating needs a specific user designed pedigree to control mating process. Pedigree should at least start with generation 2. Please make sure that paternal id and maternal id can be found in the last generation. Note that the individuals in the pedigree data file do not need to be sorted by the date of birth, and the missing value can be replaced by NA or 0. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# user-designed-pedigree mating
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          # num.ind = 100,
          mtd.reprod = "userped",
          userped = userped) # input your own pedigree
```

# Two way cross
**[back to top](#contents)**  

Two-way cross method needs two genotype matrice of different two breeds. You can input your own genotype matrix by parameters **rawgeno1** and **rawgeno2**. If any of these two is NULL, **SIMER** will generates a random one.

```r
# two-way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          rawgeno2 = NULL,
          # num.ind = 100,
          mtd.reprod = "singcro",
          num.prog = 2,
          ratio = 0.5)
```

# Three way cross
**[back to top](#contents)** 

Three-way cross method needs three genotype matrice of different three breeds. You can input your own genotype matrix by parameters **rawgeno1**, **rawgeno2**, and **rawgeno3**. If any of these three is NULL, **SIMER** will generates a random one. In triple-cross method, you can set litter size of the single-cross of population2 and population3 by **prog.tri**.

```r
# three way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          rawgeno2 = NULL,
          rawgeno3 = NULL,
          # num.ind = 100,
          mtd.reprod = "tricro",
          num.prog = 2,
          prog.tri = 2,
          ratio = 0.5)
```

# Four way cross
**[back to top](#contents)**  

Four-way cross method needs four genotype matrice of different four breeds. You can input your own genotype matrix by parameters **rawgeno1**, **rawgeno2**, **rawgeno3**, and **rawgeno4**. If any of these four is NULL, **SIMER** will generates a random one. In four-way cross method, you can set litter size of the first two two-way cross by **prog.doub**.

```r
# four way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          rawgeno2 = NULL,
          rawgeno3 = NULL,
          rawgeno4 = NULL,
          # num.ind = 100,
          mtd.reprod = "doubcro",
          num.prog = 2,
          prog.doub = 2,
          ratio = 0.5)
```

# Back cross
**[back to top](#contents)**  

Back-cross method needs two different breeds. You can input your own genotype matrix by parameters **rawgeno1** and **rawgeno2**. If any of these two is NULL, **SIMER** will generates a random one. Back-cross is similar to two-way cross but with some differences: 1. Back-cross is multi-generation mating; 2. The first base population is fixed in every generations. In back-cross method, you can set litter size of two-way cross in every generation by a vector **prog.back**.

```r
# Back-cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          out = out,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          rawgeno2 = NULL,
          # num.ind = 100,
          mtd.reprod = "backcro",
          num.prog = 2,
          prog.back = rep(2, 5),
          ratio = 0.5)
```

---

# Output
**[back to top](#contents)**  

**SIMER** outputs data including population information, marker effects, trait information, genotype, genotypic id, genotypic map, and selection intensity. 

## Population information
**[back to top](#contents)**  

Population information contains generation,  individual indice, family indice, within-family indice, sire indice, dam indice, sex, and phenotpye. 

```r
> pop <- simer.list$pop
> head(pop)
  gen index fam infam sir dam sex     pheno
1   1     1   1     1   0   0   1 -1.043716
2   1     2   2     2   0   0   1  0.551076
3   1     3   3     3   0   0   1 -5.727581
4   1     4   4     4   0   0   1  5.825591
5   1     5   5     5   0   0   1  4.430751
6   1     6   6     6   0   0   1 15.059641
```

## Marker effects
**[back to top](#contents)**  

Marker effects is a list with selected markers and effects of markers. 

```r
> effs <- simer.list$effs
> str(effs)
List of 2
 $ mrk1: int [1:18] 4006 4950 10416 18070 21899 23817 26038 28886 29437 30012 ...
 $ eff1:List of 1
  ..$ eff.a: num [1:18] -0.0144 -5.7868 -1.2436 1.2436 -3.4893 ...
```

## Trait information
**[back to top](#contents)**  

Trait information is displayed by generation in single-breed reproduction or by breed in multiple-breeds reproduction. In every generation (or breed), it contains trait information (variance components and heritability), effect information (phenotype decomposition), phenotype information (TBV, TGV, pEBVs, gEBVs, ssEBVs, phenotype), and others.  

```r
> trait <- simer.list$trait
> str(trait)
List of 5
 $ gen1:List of 4
  ..$ info.tr   :List of 3
  .. ..$ Vg: num 23.3
  .. ..$ Ve: num 63
  .. ..$ h2: num 0.27
  ..$ info.eff  :'data.frame':	40 obs. of  2 variables:
  .. ..$ ind.a  : num [1:40] -6.46 6.46 6.46 6.46 0 ...
  .. ..$ ind.env: num [1:40] 5.419 -5.912 -12.191 -0.638 4.431 ...
  ..$ info.pheno:'data.frame':	40 obs. of  3 variables:
  .. ..$ TBV  : num [1:40] -6.46 6.46 6.46 6.46 0 ...
  .. ..$ TGV  : num [1:40] -6.46 6.46 6.46 6.46 0 ...
  .. ..$ pheno: num [1:40] -1.044 0.551 -5.728 5.826 4.431 ...
  ..$ qtn.a     : num [1:18, 1:40] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:40] "V1" "V3" "V5" "V7" ...
 $ gen2:List of 4
  ..$ info.tr   :List of 3
  .. ..$ Vg: num 41.7
  .. ..$ Ve: num 113
  .. ..$ h2: num 0.27
  ..$ info.eff  :'data.frame':	32 obs. of  2 variables:
  .. ..$ ind.a  : num [1:32] 5.17 1.29 6.46 6.46 6.46 ...
  .. ..$ ind.env: num [1:32] 10.8 -8.61 -2.41 1.8 -5.04 ...
  ..$ info.pheno:'data.frame':	32 obs. of  3 variables:
  .. ..$ TBV  : num [1:32] 5.17 1.29 6.46 6.46 6.46 ...
  .. ..$ TGV  : num [1:32] 5.17 1.29 6.46 6.46 6.46 ...
  .. ..$ pheno: num [1:32] 15.97 -7.31 4.06 8.26 1.43 ...
  ..$ qtn.a     : num [1:18, 1:32] 0 -1 0 0 0 0 0 0 0 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:32] "V17" "V18" "t1" "t2" ...
 $ gen3:List of 4
  ..$ info.tr   :List of 3
  .. ..$ Vg: num 17.5
  .. ..$ Ve: num 36.5
  .. ..$ h2: num 0.324
  ..$ info.eff  :'data.frame':	26 obs. of  2 variables:
  .. ..$ ind.a  : num [1:26] -1.525 -1.525 0 0.609 6.463 ...
  .. ..$ ind.env: num [1:26] 3.29 2.31 -6.94 5.97 15.28 ...
  ..$ info.pheno:'data.frame':	26 obs. of  3 variables:
  .. ..$ TBV  : num [1:26] -1.525 -1.525 0 0.609 6.463 ...
  .. ..$ TGV  : num [1:26] -1.525 -1.525 0 0.609 6.463 ...
  .. ..$ pheno: num [1:26] 1.77 0.781 -6.936 6.576 21.743 ...
  ..$ qtn.a     : num [1:18, 1:26] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:26] "t2" "t2" "t2" "t1" ...
 $ gen4:List of 4
  ..$ info.tr   :List of 3
  .. ..$ Vg: num 15.8
  .. ..$ Ve: num 32
  .. ..$ h2: num 0.33
  ..$ info.eff  :'data.frame':	20 obs. of  2 variables:
  .. ..$ ind.a  : num [1:20] -2.808 -0.242 9.085 -2.013 4.052 ...
  .. ..$ ind.env: num [1:20] 3.96 2.64 -2.9 -2.65 2.2 ...
  ..$ info.pheno:'data.frame':	20 obs. of  3 variables:
  .. ..$ TBV  : num [1:20] -2.808 -0.242 9.085 -2.013 4.052 ...
  .. ..$ TGV  : num [1:20] -2.808 -0.242 9.085 -2.013 4.052 ...
  .. ..$ pheno: num [1:20] 1.15 2.4 6.19 -4.66 6.26 ...
  ..$ qtn.a     : num [1:18, 1:20] 1 1 1 0 0 1 1 1 0 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:20] "V18" "t2" "t1" "t1" ...
 $ gen5:List of 4
  ..$ info.tr   :List of 3
  .. ..$ Vg: num 12.1
  .. ..$ Ve: num 28.4
  .. ..$ h2: num 0.299
  ..$ info.eff  :'data.frame':	16 obs. of  2 variables:
  .. ..$ ind.a  : num [1:16] 3.79 4.377 0.551 0.551 7.68 ...
  .. ..$ ind.env: num [1:16] -1.22 -5.5 -7.99 8.79 3.1 ...
  ..$ info.pheno:'data.frame':	16 obs. of  3 variables:
  .. ..$ TBV  : num [1:16] 3.79 4.377 0.551 0.551 7.68 ...
  .. ..$ TGV  : num [1:16] 3.79 4.377 0.551 0.551 7.68 ...
  .. ..$ pheno: num [1:16] 2.57 -1.12 -7.44 9.34 10.78 ...
  ..$ qtn.a     : num [1:18, 1:16] -1 -1 -1 0 0 0 -1 -1 -1 -1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:16] "t2" "t1" "t2" "t2" ...
```

## Genotype
**[back to top](#contents)**  

In genotype matrix, each row represents a marker and each column represents a individual. 

```r
> geno <- simer.list$geno
> geno[1:6, 1:6]
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    1    1    1    1    0    0
[2,]    1    1    1    1    0    0
[3,]    1    1    1    1    0    0
[4,]    1    1    1    1    0    0
[5,]    1    1    1    1    0    0
[6,]    1    1    1    1    0    0
```

## Genotypic id
**[back to top](#contents)**  

Genotypic id is the indice of genotyped individuals. 

```r
> genoid <- simer.list$genoid
> head(genoid)
[1] 73 74 75 76 77 78
```

## Genotypic map
**[back to top](#contents)**  

In SNP map information, the first column is SNP name, the second column is Chromosome ID, the third column is phsical position, the fourth column is REF, the fifth column is ALT, the sixth column is block ID, and the seventh is recombination information ("1" represents making recombination and "0" represents not).   


```r
> map <- simer.list$map
> head(map)
     SNP          Chrom BP          REF ALT block recom
[1,] "1_10673082" "1"   " 10673082" "T" "C" "1"   "1"  
[2,] "1_10723065" "1"   " 10723065" "A" "G" "1"   "1"  
[3,] "1_11407894" "1"   " 11407894" "A" "G" "1"   "1"  
[4,] "1_11426075" "1"   " 11426075" "T" "C" "1"   "1"  
[5,] "1_13996200" "1"   " 13996200" "T" "C" "1"   "1"  
[6,] "1_14638936" "1"   " 14638936" "T" "C" "1"   "1"  
```

## Selection intensity
**[back to top](#contents)**  

Selection intensity is calculated according to ratio of selected individuals. 

```r
> si <- simer.list$si
> si
[1] 0.3499524
```

---

# Citation
**[back to top](#contents)**  

```
For SIMER:
Hope it will be coming soon!

For ADI model:
Kao, Chenhung, et al. "Modeling Epistasis of Quantitative Trait Loci Using Cockerham's Model." Genetics 160.3 (2002): 1243-1261.

For build.cov:
B. D. Ripley "Stochastic Simulation." Wiley-Interscience (1987): Page 98.
```
---

# FAQ and Hints
**[back to top](#contents)**  

:sos: **Question1:** Failing to install "devtools":

***ERROR: configuration failed for package ‘git2r’***

***removing ‘/Users/acer/R/3.4/library/git2r’***

***ERROR: dependency ‘git2r’ is not available for package ‘devtools’***

***removing ‘/Users/acer/R/3.4/library/devtools’***

:yum: **Answer:** Please try following codes in terminal:
```ssh
apt-get install libssl-dev/unstable
```

---

:sos: **Question2:** When installing packages from Github with "devtools", an error occurred:
 
 ***Error in curl::curl_fetch_disk(url, x$path, handle = handle): Problem with the SSL CA cert (path? access rights?)***
 
:yum: **Answer:** Please try following codes and then try agian.
```r
library(httr)
set_config(config(ssl_verifypeer = 0L))
```

---

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/xiaolei-lab/SIMER/issues)
 
**[back to top](#contents)**  
