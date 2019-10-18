# SIMER [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/XiaoleiLiuBio/rMVP/issues/new) 

## Data Simulation for Life Science and Breeding

<p align="center">
<a href="https://raw.githubusercontent.com/xiaolei-lab/SIMER/master/results/simer_logo.png">
<img src="results/mvp_logo.png" height="180px" width="370px">
</a>
</p>

### Authors:
Designed and Maintained by [Dong Yin](https://github.com/Foredawnbio), [Lilin Yin](https://github.com/YinLiLin), [Haohao Zhang](https://github.com/hyacz), [Xuanning Zhang](https://github.com/zxnyaonuli), and [**Xiaolei Liu**](https://github.com/XiaoleiLiuBio).  
Contributors: Zhenshuang Tang, Jingya Xu, Xinyun Li, Mengjin Zhu, Xiaohui Yuan,  Shuhong Zhao

Questions, suggestions, and bug reports are welcome and appreciated: [xiaoleiliu@mail.hzau.edu.cn](mailto:xiaoleiliu@mail.hzau.edu.cn)

### Contents
<!-- TOC updateOnSave:false -->

- [Installation](#installation)
- [Data Preparation](#data-preparation)
    - [Genotypic map](#genotypic-map)
    - [Genotype](#genotype)
    - [Designed pedigree](#designed-pedigree)
- [Data Input](#data-input)
    - [Basic](#basic)
    - [Selective](#selective)
- [Start Simulation](#start-simulation)
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

**WE STRONGLY RECOMMEND TO INSTALL MVP ON Microsoft R Open(https://mran.microsoft.com/download/)**  

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

**For more help on Windows installation, see the [wiki page](https://github.com/XiaoleiLiuBio/rMVP/wiki/R-%E8%AF%AD%E8%A8%80%E7%8E%AF%E5%A2%83%E9%85%8D%E7%BD%AE%E6%8C%87%E5%8D%97%EF%BC%88Windows%EF%BC%89) (Chinese)**

---

# Data Preparation
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

## Genotype
**[back to top](#contents)**  

Genotype data should be  Numeric (m * n, m rows and n columns, m is the number of SNPs, n is the number of individuals) format. If you have genotype data in **PLINK Binary** format (details see http://zzz.bwh.harvard.edu/plink/data.shtml#bed), **VCF** or **Hapmap**, please convert them using "MVP.Data" function in the **rMVP**(https://github.com/xiaolei-lab/rMVP).

> `genotype.txt`

<table>
<tbody>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
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
<td align="center">2</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
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

## Designed pedigree
**[back to top](#contents)**  
**SIMER** supports designed pedigree to control mating process. Designed pedigree is useful only in "userped" reproduction. 

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

**genotypic map**, SNP map information, the first column is SNP name, the second column is Chromosome ID, the third column is phsical position, the fourth column is REF, the fifth column is ALT  
**genotype**, genotype data in **Numeric** format (m * n, m rows and n columns, m is the number of SNPs, n is the number of individuals)

```r
input.map <- read.table("map.txt" , head = TRUE)
rawgeno <- read.table("genotype.txt")
```

## Selective
**[back to top](#contents)**  
If you want to control mating process by designed pedigree. 

```r
userped <- read.table("userped.txt", header = TRUE)
```

---

# Start Simulation
**[back to top](#contents)**  

After obtaining genotypic map data and genotype data, we can start our simulation.

**num.gen**, number of generations in simulation  
**replication**, replication index of simulation  
**verbose** whether to print detail  
**mrk.dense** whether markers are dense  
**out path** of output files  
**out.format** format of output, "numeric" or "plink"  
**seed.geno** random seed of genotype matrix  
**seed.map** random seed of map file  
**out.geno.gen** indice of generation of output genotype   
**out.pheno.gen** indice of generation of  output phenotype   
**rawgeno1** extrinsic genotype matrix1  
**rawgeno2** extrinsic genotype matrix2    
**rawgeno3** extrinsic genotype matrix3  
**rawgeno4** extrinsic genotype matrix4  
**num.ind** population size of base population  
**prob** weight of "0" and "1" in genotype matrix, the sum of element in vector equals 1  
**input.map** map from outside  
**len.block** length of every blocks  
**range.hot** range of exchages in hot spot block  
**range.cold** range of exchages in cold spot block  
**rate.mut** mutation rate between 1e-8 and 1e-6  
**cal.model** phenotype model with "A", "AD", "ADI"  
**h2.tr1** heritability vector of trait1, corresponding to a, d, aXa, aXd, dXa, dXd  
**num.qtn.tr1** integer or integer vector, the number of QTN in the trait1  
**var.tr1** variances of different effects, the last 5 vector elements are corrresponding to d, aXa, aXd, dXa, dXd respectively and the rest elements are corresponding to a  
**dist.qtn.tr1** distribution of QTN's effects with options: "normal", "geometry" and "gamma", vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**eff.unit.tr1** unit effect of geometric distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**shape.tr1** shape of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**scale.tr1** scale of gamma distribution of trait1, vector elements are corresponding to a, d, aXa, aXd, dXa, dXd respectively  
**multrait** whether applying pair traits with overlapping, TRUE represents applying, FALSE represents not  
**num.qtn.trn** QTN distribution matrix, diagnal elements are total QTN number of the trait, non-diagnal are QTN number of overlop qtn  
**eff.sd** a matrix with the standard deviation of QTN effects  
**gnt.cov** genetic covaiance matrix among all traits  
**env.cov** environment covaiance matrix among all traits  
**qtn.spot** QTN probability in every blocks  
**maf** Minor Allele Frequency, markers selection range is from  maf to 0.5  
**sel.crit** selection criteria with options: "TGV", "TBV", "pEBVs", "gEBVs", "ssEBVs", "pheno"  
**sel.on** whether to add selection  
**mtd.reprod** different reproduction methods with options: "clone", "dh", "selfpol", "singcro", "tricro", "doubcro", "backcro","randmate", "randexself" and "userped"  
**userped** user-specific pedigree  
**num.prog** litter size of dams  
**ratio** ratio of males in all individuals  
**prog.tri** litter size of the first single cross process in trible cross process  
**prog.doub** litter size of the first two single cross process in double cross process  
**prog.back** a vector with litter size in every generations  
**ps** fraction selected in selection  
**decr** whether to sorting with descreasing  
**sel.multi** selection method of multi-trait with options: "tdm", "indcul" and "index"  
**index.wt** economic weights of selection index method  
**index.tdm** index represents which trait is being selected. NOT CONTROL BY USER  
**goal.perc** percentage of goal more than mean of scores of individuals  
**pass.perc** percentage of expected excellent individuals  
**sel.sing** selection method of single trait with options: "ind", "fam", "infam" and "comb"  

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
              verbose = verbose, 
              out = out,
              input.map = input.map,
              rawgeno1 = rawgeno, # use your own genotype matrix
              # num.ind = 100,
              mtd.reprod = "randmate",
              num.prog = 2,
              ratio = 0.5)
}
```

---

# Output
**[back to top](#contents)**  
**SIMER** outputs data including population information, marker effects, trait information, genotype, genotypic id, genotypic map, and selection intensity. 

## Population information
**[back to top](#contents)**  

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

```r
> genoid <- simer.list$genoid
> head(genoid)
[1] 73 74 75 76 77 78
```

## Genotypic map
**[back to top](#contents)**  

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
B. D. Ripley (1987) <ISBN:9780470009604>.
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

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/xiaolei-lab/SIMERP/issues)
 
**[back to top](#contents)**  
