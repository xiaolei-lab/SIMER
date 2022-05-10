# SIMER 
[![GitHub issues](https://img.shields.io/github/issues/xiaolei-lab/SIMER?color=green)](https://github.com/xiaolei-lab/SIMER/issues/new) [![](https://img.shields.io/badge/GitHub-0.9.0.0-blueviolet.svg)](https://github.com/xiaolei-lab/SIMER) <a href="https://hits.seeyoufarm.com"/><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fxiaolei-lab%2FSIMER"/></a>

## Data Simulation for Life Science and Breeding

### Authors:
***Design and Maintenance:*** Dong Yin, Xuanning Zhang, Lilin Yin ,Haohao Zhang, and ***Xiaolei Liu***.  
***Contributors:*** Zhenshuang Tang, Jingya Xu, Xiaohui Yuan, Xinyun Li, and Shuhong Zhao.

***If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/xiaolei-lab/SIMER/issues/new):point_left:.***

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
    - [Quick Start for Complete Simulation](#quick-start-for-complete-simulation)
    - [Quick Start for Annotation Simulation](#quick-start-for-annotation-simulation)
    - [Quick Start for Genotype Simulation](#quick-start-for-genotype-simulation)
    - [Quick Start for Phenotype Simulation](#quick-start-for-phenotype-simulation)
- [Genotype Simulation](#genotype-simulation)
    - [Gallery of genotype simulation parameters](#gallery-of-genotype-simulation-parameters)
    - [Generate genotype matrix of base population](#generate-genotype-matrix-of-base-population)
    - [Add chromosome crossovers and mutations to genotype matrix](#add-chromosome-crossovers-and-mutations-to-genotype-matrix)
- [Phenotype Simulation](#phenotype-simulation)
    - [Gallery of phenotype simulation parameters](#gallery-of-phenotype-simulation-parameters)
    - [Generate phenotype by A model](#generate-phenotype-by-A-model)
    - [Generate phenotype by AD model](#generate-phenotype-by-AD-model)
    - [Generate phenotype by GxG model](#generate-phenotype-by-GxG-model)
    - [Generate phenotype by Repeated Record model](#generate-phenotype-by-repeated-record-model)
    - [Generate phenotype controlled by QTNs subject to Normal distribution](#generate-phenotype-controlled-by-QTNs-subject-to-normal-distribution)
    - [Generate phenotype controlled by QTNs subject to Geometric distribution](#generate-phenotype-controlled-by-QTNs-subject-to-geometric-distribution)
    - [Generate phenotype controlled by QTNs subject to Gamma distribution](#generate-phenotype-controlled-by-QTNs-subject-to-gamma-distribution)
    - [Generate phenotype controlled by QTNs subject to Beta distribution](#generate-phenotype-controlled-by-QTNs-subject-to-beta-distribution)
    - [Generate phenotype with fixed effect and environmental random effect](#generate-phenotype-with-fixed-effect-and-environmental-random-effect)
    - [Generate phenotype by GxE model](#generate-phenotype-by-GxE-model)
- [Population Simulation of Multiple-Generation with Genotype and Phenotype](#population-simulation-of-multiple-generation-with-genotype-and-phenotype)
    - [Gallery of population simulation parameters](#gallery-of-population-simulation-parameters)
    - [Individual selection on single trait](#individual-selection-on-single-trait)
    - [Family selection on single trait](#family-selection-on-single-trait)
    - [Within family selection on single trait](#within-family-selection-on-single-trait)
    - [Combined selection on single trait](#combined-selection-on-single-trait)
    - [Tandem selection on multiple traits](#tandem-selection-on-multiple-traits)
    - [Independent culling selection on multiple traits](#independent-culling-selection-on-multiple-traits)
    - [Index selection on multiple traits](#index-selection-on-multiple-traits)
    - [Clone for plant](#clone-for-plant)
    - [Double haploid for plant](#double-haploid-for-plant)
    - [Self pollination for plant and micro-organism](#self-pollination-for-plant-and-micro-organism)
    - [Random mating for plant and animal](#random-mating-for-plant-and-animal)
    - [Random mating excluding self pollination for animal](#random-mating-excluding-self-pollination-for-animal)
    - [Two way cross for animal](#two-way-cross-for-animal)
    - [Three way cross for animal](#three-way-cross-for-animal)
    - [Four way cross for animal](#four-way-cross-for-animal)
    - [Back cross for animal](#back-cross-for-animal)
    - [User designed pedigree mating for plant and animal](#user-designed-pedigree-mating-for-plant-and-animal)
    - [AN EASY WAY TO GENERATE A POPULATION](#an-easy-way-to-generate-a-population)
- [Breeding Program Design](#breeding-program-design)
    - [Gallery of breeding program design parameters](#gallery-of-breeding-program-design-parameters)
    - [Breeding program design preparation](#breeding-program-design-preparation)
    - [Breeding plan comparison](#breeding-plan-comparison)
- [Global Options](#global-options)
    - [Gallery of global input parameters](#gallery-of-global-input-parameters)
    - [Counts of total population size](#counts-of-total-population-size)
    - [Simulation of multiple populations](#simulation-of-multiple-populations)
    - [File output](#file-output)
    - [Generation selective output](#generation-selective-output)
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

**WE STRONGLY RECOMMEND TO INSTALL SIMER ON Microsoft R Open (https://mran.microsoft.com/download/)**.  

## Installation

* The latest version: 
```R
devtools::install_github("xiaolei-lab/SIMER")
```

After installed successfully, **```SIMER```** can be loaded by typing
```r
> library(simer)
```
Typing ```?simer``` could get the details of all parameters.

---

# Data Preparation

## Genotype
**[back to top](#contents)**  

***Genotype data*** should be ***Numeric*** format (***m*** rows and ***n*** columns, ***m*** is the number of SNPs, ***n*** is the number of individuals). Other ***genotype data*** such as ***PLINK Binary*** format (details see http://zzz.bwh.harvard.edu/plink/data.shtml#bed), ***VCF*** or ***Hapmap*** can be converted to ***Numeric*** format using ```MVP.Data``` function in the **```rMVP```** (https://github.com/xiaolei-lab/rMVP).

> `genotype.txt`

<table>
<tbody>
<tr>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">2</td>
<td align="center">…</td>
<td align="center">0</td>
</tr></tbody></table>

## Genotypic map
**[back to top](#contents)**  
***Genotypic Map*** is necessary in **```SIMER```**. The first column is ***SNP name***, the second column is ***Chromosome ID***, the third column is ***physical position***, the fourth column is ***REF***, and the fifth column is ***ALT***.  It will be used to generate ***annotation data***, ***genotype data***, and ***phenotype data***. 

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
**```SIMER```** supports ***user designed pedigree*** to control mating process. ***User designed pedigree*** is useful only in ```userped``` reproduction. The first column is ***sample id***, the second column is ***paternal id***, and the third column is ***maternal id***. Please make sure that ***paternal id*** and ***maternal id*** can match to ***genotype data***. 

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
At least users should prepare two datasets: ***genotypic map*** and ***genotype data***.  

***genotype data***, ***Numeric*** format (***m*** rows and ***n*** columns, ***m*** is the number of SNPs, ***n*** is the number of individuals)  
***genotypic map***, SNP map information, the first column is ***SNP name***, the second column is ***Chromosome ID***, the third column is ***physical position***, the fourth column is ***REF***, and the fifth column is ***ALT***.  

```r
pop.geno <- read.table("genotype.txt")
pop.map <- read.table("map.txt" , head = TRUE)
```

## Optional
**[back to top](#contents)**  
Mating process can be designed by ***user designed pedigree***. 

***pedigree***, pedigree information, the first column is ***sample id***, the second column is ***paternal id***, and the third column is ***maternal id***. Note that the individuals in the ***pedigree*** do not need to be sorted by the date of birth, and the missing value can be replaced by NA or 0.

```r
userped <- read.table("userped.txt", header = TRUE)
```

---

# Quick Start
**[back to top](#contents)** 

All simulation processes can be devided into 2 steps: ***1) generate simulation parameters***; ***2) run simulation process***.

## Quick Start for Complete Simulation
**[back to top](#contents)**

A quick start for ***Complete Simulation*** is shown below:

```r
# Generate all simulation parameters
SP <- param.simer(out = "simer")

# Run Simer
SP <- simer(SP)
```

## Quick Start for Annotation Simulation
**[back to top](#contents)** 

A quick start for ***Annotation Simulation*** is shown below:

```r
# Real genotypic map
# pop.map <- read.table("Real_Genotypic_map.txt", header = TRUE)
# Simulated genotypic map
pop.map <- generate.map(pop.marker = 1e4)

# Generate annotation simulation parameters
SP <- param.annot(pop.map = pop.map, qtn.num = 10)

# Run annotation simulation
SP <- annotation(SP)
```

## Quick Start for Genotype Simulation
**[back to top](#contents)** 

A quick start for ***Genotype Simulation*** is shown below:

```r
# Generate genotype simulation parameters
SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)

# Run genotype simulation
SP <- genotype(SP)
```

## Quick Start for Phenotype Simulation
**[back to top](#contents)** 

A quick start for ***Phenotype Simulation*** is shown below:

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

---

# Genotype Simulation
**[back to top](#contents)** 

***Genotype data*** in **```SIMER```** will be generated randomly or from outside genotype matrix. Chromosome crossovers and base mutations depend on block information and recombination information of ***Annotation data***. 

## Gallery of genotype simulation parameters
**[back to top](#contents)** 

```genotype```, main function of ***Genotype Simulation***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>pop.geno</b></td>
<td>NULL</td>
<td>big.matrix or matrix</td>
<td>the genotype data.</td>
</tr>
<tr>
<td><b>incols</b></td>
<td>1</td>
<td>1 or 2</td>
<td>'1': one-column genotype represents an individual; '2': two-column genotype represents an individual.</td>
</tr>
<tr>
<td><b>pop.marker</b></td>
<td>1e4</td>
<td>num</td>
<td>the number of markers.</td>
</tr>
<tr>
<td><b>pop.ind</b></td>
<td>1e2</td>
<td>num</td>
<td>the number of individuals in the base population.</td>
</tr>
<tr>
<td><b>prob</b></td>
<td>NULL</td>
<td>num vector</td>
<td>the genotype code probability.</td>
</tr>
<tr>
<td><b>rate.mut</b></td>
<td>1e-8</td>
<td>num</td>
<td>the mutation rate of the genotype data.</td>
</tr>
</tbody>
</table>

```annotation```, main function of ***Annotation Simulation***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>recom.spot</b></td>
<td>FALSE</td>
<td>TRUE or FALSE</td>
<td>whether to generate recombination events.</td>
</tr>
<tr>
<td><b>range.hot</b></td>
<td>4:6</td>
<td>num vector</td>
<td>the recombination times range in the hot spot.</td>
</tr>
<tr>
<td><b>range.cold</b></td>
<td>1:5</td>
<td>num vector</td>
<td>the recombination times range in the cold spot.</td>
</tr>
</tbody>
</table>

## Generate genotype matrix of base population
**[back to top](#contents)** 

There are two different ways to generate genotype matrix of base population.

```r
### 01 Use Genotype Data from Outside ###
# Create a genotype matrix
pop.geno <- matrix(0, nrow = 1e4, ncol = 1e2)

# Generate genotype simulation parameters
SP <- param.geno(pop.geno = pop.geno)

# Run genotype simulation
SP <- genotype(SP)

### 02 Create Genotype Data Randomly ###
# Generate genotype simulation parameters
SP <- param.geno(pop.marker = 1e4, pop.ind = 1e2)

# Run genotype simulation
SP <- genotype(SP)
```

## Add chromosome crossovers and mutations to genotype matrix
**[back to top](#contents)** 

With ***annotation data***, chromosome crossovers and mutations can be added to genotype matrix.  

```r
# Generate annotation simulation parameters
# If recom.spot = TRUE, chromsome crossovers will be added to genotype matrix
SP <- param.annot(recom.spot = TRUE)
# Generate genotype simulation parameters
# Base mutation rate is 1e8
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2, rate.mut = 1e-8)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
```

Note that recombination only exists in meiosis. Therefore, some reproduction methods such as ```clone``` do not have recombination process. Users can set ```recom.spot = FALSE``` to add only mutations to genotype matrix.  

```r
# Generate annotation simulation parameters
# If recom.spot = FALSE, chromsome crossovers will not be added to genotype matrix
SP <- param.annot(recom.spot = FALSE)
# Generate genotype simulation parameters
# Base mutation rate is 1e8
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2, rate.mut = 1e-8)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
```             

---

# Phenotype Simulation
**[back to top](#contents)**  

***Phenotype data*** in **```SIMER```** will be generated according to different models, including:  
**(1)** Single-Trait Model  
**(2)** Multiple-Trait Model  
**(3)** Repeated Record Model  
**(4)** Genetic Effect Model (***A***dditive effect, ***D***ominant effect, and ***G***enetic-***G***enetic interaction effect)  
**(5)** Effect Distribution Model (QTN effect distribution: ****Norm***al distribution, ***Geom***etric distribution, ***Gamma*** distribution, and ***Beta*** distribution)  
**(6)** Linear Mixed Model (***F***ixed effect, ***E***nvironmental ***R***andom effect, ***G***enetic ***R***andom effect, and ***G***enetic-***E***nvironmental effect)  

## Gallery of phenotype simulation parameters
**[back to top](#contents)**  

```phenotype```, main function of ***Phenotype Simulation***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>pop</b></td>
<td>NULL</td>
<td>data.frame</td>
<td>the population information containing environmental factors and other effects.</td>
</tr>
<tr>
<td><b>pop.ind</b></td>
<td>100</td>
<td>num</td>
<td>the number of individuals in the base population.</td>
</tr>
<tr>
<td><b>pop.rep</b></td>
<td>1</td>
<td>num</td>
<td>the repeated times of repeated records.</td>
</tr>
<tr>
<td><b>pop.rep.bal</b></td>
<td>TRUE</td>
<td>TRUE or FALSE</td>
<td>whether repeated records are balanced.</td>
</tr>
<tr>
<td><b>pop.env</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of environmental factors setting.</td>
</tr>
<tr>
<td><b>phe.model</b></td>
<td>list(tr1 = "T1 = A + E")</td>
<td>list</td>
<td>a list of genetic model of phenotype such as "T1 = A + E".</td>
</tr>
<tr>
<td><b>phe.h2A</b></td>
<td>list(tr1 = 0.3)</td>
<td>list</td>
<td>a list of additive heritability.</td>
</tr>
<tr>
<td><b>phe.h2D</b></td>
<td>list(tr1 = 0.1)</td>
<td>list</td>
<td>a list of dominant heritability.</td>
</tr>
<tr>
<td><b>phe.h2GxG</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of GxG interaction heritability.</td>
</tr>
<tr>
<td><b>phe.h2GxE</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of GxE interaction heritability.</td>
</tr>
<tr>
<td><b>phe.h2PE</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of permanent environmental heritability.</td>
</tr>
<tr>
<td><b>phe.var</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of phenotype variance.</td>
</tr>
<tr>
<td><b>phe.corA</b></td>
<td>NULL</td>
<td>matrix</td>
<td>the additive genetic correlation matrix.</td>
</tr>
<tr>
<td><b>phe.corD</b></td>
<td>NULL</td>
<td>matrix</td>
<td>the dominant genetic correlation matrix.</td>
</tr>
<tr>
<td><b>phe.corGxG</b></td>
<td>NULL</td>
<td>list</td>
<td>a list of the GxG genetic correlation matrix.</td>
</tr>
<tr>
<td><b>phe.corPE</b></td>
<td>NULL</td>
<td>matrix</td>
<td>the permanent environmental correlation matrix.</td>
</tr>
<tr>
<td><b>phe.corE</b></td>
<td>NULL</td>
<td>matrix</td>
<td>the residual correlation matrix.</td>
</tr>
</tbody>
</table>

```annotation```, main function of ***Annotation Simulation***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>pop.map</b></td>
<td>NULL</td>
<td>data.frame</td>
<td>the map data with annotation information.</td>
</tr>
<tr>
<td><b>qtn.num</b></td>
<td>10</td>
<td>num, vector or matrix</td>
<td>integer: the QTN number of single trait; vector: the multiple group QTN number of single trait; matrix: the QTN number of multiple traits.</td>
</tr>
<tr>
<td><b>qtn.model</b></td>
<td>'A'</td>
<td>character</td>
<td>the genetic model of QTN such as 'A + D'.</td>
</tr>
<tr>
<td><b>qtn.dist</b></td>
<td>list(tr1 = 'norm')</td>
<td>list</td>
<td>the QTN distribution containing 'norm', 'geom', 'gamma' or 'beta'.</td>
</tr>
<tr>
<td><b>qtn.sd</b></td>
<td>list(tr1 = 1)</td>
<td>list</td>
<td>the standard deviations for normal distribution.</td>
</tr>
<tr>
<td><b>qtn.prob</b></td>
<td>NULL</td>
<td>list</td>
<td>the probability of success for geometric distribution.</td>
</tr>
<tr>
<td><b>qtn.shape</b></td>
<td>NULL</td>
<td>list</td>
<td>the shape parameter for gamma distribution.</td>
</tr>
<tr>
<td><b>qtn.scale</b></td>
<td>NULL</td>
<td>list</td>
<td>the scale parameter for gamma distribution.</td>
</tr>
<tr>
<td><b>qtn.shape1</b></td>
<td>NULL</td>
<td>list</td>
<td>the shape1 parameter for beta distribution.</td>
</tr>
<tr>
<td><b>qtn.shape2</b></td>
<td>NULL</td>
<td>list</td>
<td>the shape2 parameter for beta distribution.</td>
</tr>
<tr>
<td><b>qtn.ncp</b></td>
<td>NULL</td>
<td>list</td>
<td>the ncp parameter for beta distribution.</td>
</tr>
<tr>
<td><b>qtn.spot</b></td>
<td>NULL</td>
<td>list</td>
<td>the QTN distribution probability in each block.</td>
</tr>
<tr>
<td><b>len.block</b></td>
<td>5e7</td>
<td>num</td>
<td>the block length.</td>
</tr>
<tr>
<td><b>maf</b></td>
<td>NULL</td>
<td>num</td>
<td>the maf threshold, markers less than this threshold will be exclude.</td>
</tr>
</tbody>
</table>

## Generate phenotype by A model
**[back to top](#contents)** 

In ***A*** model, **```SIMER```** only considers ***A***dditive effect as genetic effect. Users should prepare ***A***dditive ***QTN*** effect in the ***Annotation data*** for generating ***A***dditive ***I***ndividual effect. ***A***dditive single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

In multiple-trait simulation, **```SIMER```** can build ***accurate Additive genetic correlation*** between multiple traits. ***A***dditive multiple-trait simulation is displayed as follows: 


```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

In single-trait simulation, the trait can be controlled by ***multiple-group QTNs***. An example of single-trait controlled by two-group QTNs is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = c(2, 8), qtn.model = "A") # Group1: 2 QTNs; Group 2: 8 QTNs
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype by AD model
**[back to top](#contents)** 

In ***AD*** model, **```SIMER```** considers ***A***dditive effect and ***D***ominant effect as genetic effect. Users should prepare ***A***dditive ***QTN*** effect and ***D***ominant ***QTN*** effect in the ***Annotation data*** for generating ***A***dditive ***I***ndividual effect and ***D***ominant ***I***ndividual effect. ***A***dditive and ***D***ominant single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A + D") # Additive effect and Dominant effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + D + E"), # "T1" (Trait 1) consists of Additive effect, Dominant effect, and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3),
  phe.h2D = list(tr1 = 0.1)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

In multiple-trait simulation, **```SIMER```** can build ***accurate Additive genetic correlation*** and ***accurate Dominant genetic correlation*** between multiple traits. ***A***dditive and ***D***ominant multiple-trait simulation is displayed as follows: 


```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A + D") # Additive effect and Dominant effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + D + E", # "T1" (Trait 1) consists of Additive effect, Dominant effect, and Residual effect
    tr2 = "T2 = A + D + E"  # "T2" (Trait 2) consists of Additive effect, Dominant effect, and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.h2D = list(tr1 = 0.1, tr2 = 0.1),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2), # Additive genetic correlation
  phe.corD = matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Dominant genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype by GxG model
**[back to top](#contents)** 

In ***GxG*** model, **```SIMER```** considers ***G***enetic-***G***enetic effect as genetic effect. Users should prepare ***G***enetic-***G***enetic ***QTN*** effect in the ***Annotation data*** for generating ***G***enetic-***G***enetic ***I***ndividual effect. An example of ***A***dditive-***D***ominant interaction in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A + D + A:D") # Additive effect, Dominant effect, and Additive-Dominant interaction effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + D + A:D + E"), # "T1" (Trait 1) consists of Additive effect, Dominant effect, Additive-Dominant interaction effect, and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3),
  phe.h2D = list(tr1 = 0.1),
  phe.h2GxG = list(tr1 = list("A:D" = 0.1))
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

In multiple-trait simulation, **```SIMER```** can build ***accurate Genetic-Genetic interaction correlation*** between multiple traits. An example of ***A***dditive-***D***ominant interaction in multiple-trait simulation is displayed as follows: 


```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A + D + A:D") # Additive effect, Dominant effect, and Additive-Dominant interaction effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + D + A:D + E", # "T1" (Trait 1) consists of Additive effect, Dominant effect, Additive-Dominant interaction effect, and Residual effect
    tr2 = "T2 = A + D + A:D + E"  # "T2" (Trait 2) consists of Additive effect, Dominant effect, Additive-Dominant interaction effect, and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.h2D = list(tr1 = 0.1, tr2 = 0.1),
  phe.h2GxG = list(tr1 = list("A:D" = 0.1), tr2 = list("A:D" = 0.1)),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2),                 # Additive genetic correlation
  phe.corD = matrix(c(1, 0.5, 0.5, 1), 2, 2),                 # Dominant genetic correlation
  phe.corGxG = list("A:D" = matrix(c(1, 0.5, 0.5, 1), 2, 2))  # Additive-Dominant interaction genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype by Repeated Record model
**[back to top](#contents)** 

In ***Repeated Record*** model, **```SIMER```** adds ***PE*** (***P***ermanent ***E***nvironmental) effect to the phenotype. The number of repeated records can be set by ```pop.rep```. In the meantime, ```pop.rep.bal``` can be used to determine whether repeated records are balanced. ***Repeated Record*** in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  pop.rep = 2,                          # The number of repeated records is 2
  pop.rep.bal = TRUE,                   # Repeated records are balanced
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

In multiple-trait simulation, **```SIMER```** can build ***accurate Permanent Environmental correlation*** between multiple traits. ***Repeated Record*** in multiple-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  pop.rep = 2,          # The number of repeated records is 2
  pop.rep.bal = TRUE,   # Repeated records are balanced
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2), # Additive genetic correlation
  phe.corPE = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Permanent Environmental correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype controlled by QTNs subject to Normal distribution
**[back to top](#contents)** 

***Norm***al distribution is the most common QTN effect distribution. Phenotype controlled by QTNs subject to ***Norm***al distribution in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = 10,
  qtn.model = "A",
  qtn.dist = list(tr1 = "norm"),
  qtn.sd = list(tr1 = 1)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

Phenotype controlled by QTNs subject to ***Norm***al distribution in multiple-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = matrix(c(6, 4, 4, 6), 2, 2),
  qtn.model = "A",
  qtn.dist = list(tr1 = "norm", tr2 = "norm"),
  qtn.sd = list(tr1 = 1, tr2 = 1)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype controlled by QTNs subject to Geometric distribution
**[back to top](#contents)** 

***Geom***etric distribution is the probability of success for the first time obtained only after K trials among the N Bernoulli trials. ***Geom***etric distribution can be used as a QTN effect distribution. Phenotype controlled by QTNs subject to ***Geom***etric distribution in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = 10,
  qtn.model = "A",
  qtn.dist = list(tr1 = "geom"),
  qtn.prob = list(tr1 = 0.5)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

Phenotype controlled by QTNs subject to ***Geom***etric distribution in multiple-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = matrix(c(6, 4, 4, 6), 2, 2),
  qtn.model = "A",
  qtn.dist = list(tr1 = "geom", tr2 = "geom"),
  qtn.prob = list(tr1 = 0.5, tr2 = 0.5)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype controlled by QTNs subject to Gamma distribution
**[back to top](#contents)** 

***Gamma*** distribution is the sum of N independent exponential random variables. Note that ***Exp***onential distribution is a special form of ***Gamma*** distribution when ```qtn.shape = 1``` and ```qtn.scale = 1```. Phenotype controlled by QTNs subject to ***Gamma*** distribution in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = 10,
  qtn.model = "A",
  qtn.dist = list(tr1 = "gamma"),
  qtn.shape = list(tr1 = 1),
  qtn.scale = list(tr1 = 1)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

Phenotype controlled by QTNs subject to ***Gamma*** distribution in multiple-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = matrix(c(6, 4, 4, 6), 2, 2),
  qtn.model = "A",
  qtn.dist = list(tr1 = "gamma", tr2 = "gamma"),
  qtn.shape = list(tr1 = 1, tr2 = 1),
  qtn.scale = list(tr1 = 1, tr2 = 1)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype controlled by QTNs subject to Beta distribution
**[back to top](#contents)** 

***Beta*** distribution is a density function of conjugate prior distribution as Bernoulli distribution and Binomial distribution. Phenotype controlled by QTNs subject to ***Beta*** distribution in single-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = 10,
  qtn.model = "A",
  qtn.dist = list(tr1 = "beta"),
  qtn.shape1 = list(tr1 = 1),
  qtn.shape2 = list(tr1 = 1),
  qtn.ncp = list(tr1 = 0)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(tr1 = "T1 = A + E"), # "T1" (Trait 1) consists of Additive effect and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

Phenotype controlled by QTNs subject to ***Beta*** distribution in multiple-trait simulation is displayed as follows: 

```r
# Generate annotation simulation parameters
SP <- param.annot(
  qtn.num = matrix(c(6, 4, 4, 6), 2, 2),
  qtn.model = "A",
  qtn.dist = list(tr1 = "beta", tr2 = "beta"),
  qtn.shape1 = list(tr1 = 1, tr2 = 1),
  qtn.shape2 = list(tr1 = 1, tr2 = 1),
  qtn.ncp = list(tr1 = 0, tr2 = 0)
)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  phe.model = list(
    tr1 = "T1 = A + E", # "T1" (Trait 1) consists of Additive effect and Residual effect
    tr2 = "T2 = A + E"  # "T2" (Trait 2) consists of Additive effect and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype with fixed effect and environmental random effect
**[back to top](#contents)** 

**```SIMER```** supports add ***F***ixed effects and ***E***nvironmental ***R***andom effects to phenotype. Users should prepare a list of environmental factors setting. The effect value of ***F***ixed effect is determined by ```eff``` and the variance ratio to phenotype variance (similar to heritability) of ***E***nvironmental ***R***andom is set by ```ratio```. Phenotype with ***F***ixed effect and ***E***nvironmental ***R***andom effect in single-trait simulation is displayed as follows: 

```r
# Prepare environmental factor list
pop.env <- list(
  F1 = list( # fixed effect 1
    level = c("1", "2"),
    eff = list(tr1 = c(50, 30))
  ), 
  F2 = list( # fixed effect 2
    level = c("d1", "d2", "d3"),
    eff = list(tr1 = c(10, 20, 30))
  ),
  R1 = list( # random effect 1
    level = c("l1", "l2", "l3"),
    ratio = list(tr1 = 0.1)
  )
)

# Generate genotype simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A")
# Generate annotation simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP, 
  pop.ind = 100,
  pop.env = pop.env,
  phe.model = list(tr1 = "T1 = A + F1 + F2 + R1 + E"), # "T1" (Trait 1) consists of Additive effect, F1, F2, R1, and Residual effect
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3)
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

Phenotype with ***F***ixed effect and ***E***nvironmental ***R***andom effect in multiple-trait simulation is displayed as follows: 

```r
# Prepare environmental factor list
pop.env <- list(
  F1 = list( # fixed effect 1
    level = c("1", "2"),
    eff = list(tr1 = c(50, 30), tr2 = c(50, 30))
  ), 
  F2 = list( # fixed effect 2
    level = c("d1", "d2", "d3"),
    eff = list(tr1 = c(10, 20, 30), tr2 = c(10, 20, 30))
  ),
  R1 = list( # random effect 1
    level = c("l1", "l2", "l3"),
    ratio = list(tr1 = 0.1, tr2 = 0.1)
  )
)

# Generate genotype simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A")
# Generate annotation simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP, 
  pop.ind = 100,
  pop.env = pop.env,
  phe.model = list(
    tr1 = "T1 = A + F1 + F2 + R1 + E", # "T1" (Trait 1) consists of Additive effect, F1, F2, R1, and Residual effect
    tr2 = "T2 = A + F1 + F2 + R1 + E"  # "T2" (Trait 1) consists of Additive effect, F1, F2, R1, and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

## Generate phenotype by GxE model
**[back to top](#contents)** 

In ***GxE*** model, **```SIMER```** adds ***G***enetic-***E***nvironmental interaction effect to phenotype. Users should prepare ***G***enetic ***QTN*** effect in the ***Annotation data*** and environmental factor by ```pop.env``` for generating ***G***enetic-***E***nvironmental ***I***ndividual effect. An example of ***G***enetic-***E***nvironmental interaction in single-trait simulation is displayed as follows: 

```r
# Prepare environmental factor list
pop.env <- list(
  F1 = list( # fixed effect 1
    level = c("1", "2"),
    eff = list(tr1 = c(50, 30))
  ), 
  F2 = list( # fixed effect 2
    level = c("d1", "d2", "d3"),
    eff = list(tr1 = c(10, 20, 30))
  ),
  R1 = list( # random effect 1
    level = c("l1", "l2", "l3"),
    ratio = list(tr1 = 0.1)
  )
)

# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10, qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  pop.env = pop.env,
  phe.model = list(
    tr1 = "T1 = A + F1 + F2 + R1 + A:F1 + E" # "T1" (Trait 1) consists of Additive effect, F1, F2, R1, Additive-F1 interaction effect, and Residual effect
  ),
  # phe.var = list(tr1 = 100),
  phe.h2A = list(tr1 = 0.3),
  phe.h2GxE = list(tr1 = list("A:F1" = 0.1))
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

An example of ***G***enetic-***E***nvironmental interaction in multiple-trait simulation is displayed as follows: 

```r
# Prepare environmental factor list
pop.env <- list(
  F1 = list( # fixed effect 1
    level = c("1", "2"),
    eff = list(tr1 = c(50, 30), tr2 = c(50, 30))
  ), 
  F2 = list( # fixed effect 2
    level = c("d1", "d2", "d3"),
    eff = list(tr1 = c(10, 20, 30), tr2 = c(10, 20, 30))
  ),
  R1 = list( # random effect 1
    level = c("l1", "l2", "l3"),
    ratio = list(tr1 = 0.1, tr2 = 0.1)
  )
)

# Generate annotation simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A") # Additive effect
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP,
  pop.ind = 100,
  pop.env = pop.env,
  phe.model = list(
    tr1 = "T1 = A + F1 + F2 + R1 + A:F1 + E", # "T1" (Trait 1) consists of Additive effect, F1, F2, R1, Additive-F1 interaction effect, and Residual effect
    tr2 = "T2 = A + F1 + F2 + R1 + A:F1 + E"  # "T2" (Trait 2) consists of Additive effect, F1, F2, R1, Additive-F1 interaction effect, and Residual effect
  ),
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.h2A = list(tr1 = 0.3, tr2 = 0.3),
  phe.corA = matrix(c(1, 0.5, 0.5, 1), 2, 2) # Additive genetic correlation
)

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
```

---

# Population Simulation of Multiple-Generation with Genotype and Phenotype
**[back to top](#contents)**  

**```SIMER```** imitates the reproductive process of organisms to generate ***Multiple-Generation*** population. The ***genotype data*** and ***phenotype data*** of the population are screened by ***single-trait selection*** or ***multiple-trait selection***, and then amplified by ***species-specific reproduction***.

## Gallery of population simulation parameters
**[back to top](#contents)**  

```selects```, main function of ***Selection***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>pop.sel</b></td>
<td>NULL</td>
<td>list</td>
<td>the selected males and females.</td>
</tr>
<tr>
<td><b>ps</b></td>
<td>c(0.8, 0.8)</td>
<td>num vector</td>
<td>if ps <= 1, fraction selected in selection of males and females; if ps > 1, ps is number of selected males and females.</td>
</tr>
<tr>
<td><b>decr</b></td>
<td>TRUE</td>
<td>TRUE or FALSE</td>
<td>whether the sort order is decreasing.</td>
</tr>
<tr>
<td><b>sel.crit</b></td>
<td>'pheno'</td>
<td>character</td>
<td>the selection criteria, it can be 'TBV', 'TGV', and 'pheno'.</td>
</tr>
<tr>
<td><b>sel.single</b></td>
<td>'comb'</td>
<td>character</td>
<td>the single-trait selection method, it can be 'ind', 'fam', 'infam', and 'comb'.</td>
</tr>
<tr>
<td><b>sel.multi</b></td>
<td>'index'</td>
<td>character</td>
<td>the multiple-trait selection method, it can be 'index', 'indcul', and 'tmd'.</td>
</tr>
<tr>
<td><b>index.wt</b></td>
<td>c(0.5, 0.5)</td>
<td>num vector</td>
<td>the weight of each trait for multiple-trait selection.</td>
</tr>
<tr>
<td><b>index.tdm</b></td>
<td>1</td>
<td>num</td>
<td>the index of tandem selection for multiple-trait selection.</td>
</tr>
<tr>
<td><b>goal.perc</b></td>
<td>0.1</td>
<td>num</td>
<td>the percentage of goal more than the mean of scores of individuals.</td>
</tr>
<tr>
<td><b>pass.perc</b></td>
<td>0.9</td>
<td>num</td>
<td>the percentage of expected excellent individuals.</td>
</tr>
</tbody>
</table>

```reproduces```, main function of ***Reproduction***:  

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>pop.gen</b></td>
<td>2</td>
<td>num</td>
<td>the generations of simulated population.</td>
</tr>
<tr>
<td><b>reprod.way</b></td>
<td>'randmate'</td>
<td>character</td>
<td>reproduction method, it consists of 'clone', 'dh', 'selfpol', 'randmate', 'randexself', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.</td>
</tr>
<tr>
<td><b>sex.rate</b></td>
<td>0.5</td>
<td>num</td>
<td>the male rate in the population.</td>
</tr>
<tr>
<td><b>prog</b></td>
<td>2</td>
<td>num</td>
<td>the progeny number of an individual.</td>
</tr>
</tbody>
</table>

## Individual selection on single trait
**[back to top](#contents)**  

***Individual selection*** is a selecting method according to the ***phenotype*** of individual traits, also known as mixed selection or collective selection. This selection method is simple and easy to used for traits with ***high heritability***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "ind")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Family selection on single trait
**[back to top](#contents)** 

***Family selection*** is a selecting method by family based on the ***average of the family***. This selection method is used for traits with ***low heritability***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "fam")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Within family selection on single trait
**[back to top](#contents)** 

***Within-family*** selection is a selecting method according to the ***deviation of individual phenotype and family mean value in each family***. This selection method is used for traits with ***low heritability and small family***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "infam")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Combined selection on single trait
**[back to top](#contents)**  

***Combined selection*** is a selecting method according to ***weighed combination of the deviation of individual phenotype and family mean value***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Tandem selection on multiple traits
**[back to top](#contents)**  

***Tandem selection*** is a method for ***sequentially selecting a plurality of target traits one by one***. The index of selected trait is ```index.tdm``` and this parameter should ***not controlled by Users***.

```r
# Generate genotype simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A")
# Generate annotation simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP, 
  pop.ind = 100,
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.model = list(
    tr1 = "T1 = A + E",
    tr2 = "T2 = A + E"
  )
)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.multi = "tdm")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Independent culling selection on multiple traits
**[back to top](#contents)**  

After setting a ***minimum selection criterion*** for each target trait. ***Independent culling selection*** will ***eliminate*** this individual when the candidate's performance on any trait is ***lower than the corresponding criteria***.

```r
# Generate genotype simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A")
# Generate annotation simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP, 
  pop.ind = 100,
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.model = list(
    tr1 = "T1 = A + E",
    tr2 = "T2 = A + E"
  )
)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.multi = "indcul")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Index selection on multiple traits
**[back to top](#contents)**  

***Index selection*** is a comprehensive selection that will consider several traits based on their respective ***heritabilities***, ***phenotypic variances***, ***economic weights***, corresponding ***genetic correlations***, and ***phenotypes***. Then calculate the ***index value of each trait***, and eliminate or select it according to its level. Users can set the weight of each trait by ```index.wt```.

```r
# Generate genotype simulation parameters
SP <- param.annot(qtn.num = matrix(c(6, 4, 4, 6), 2, 2), qtn.model = "A")
# Generate annotation simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(
  SP = SP, 
  pop.ind = 100,
  # phe.var = list(tr1 = 100, tr2 = 100),
  phe.model = list(
    tr1 = "T1 = A + E",
    tr2 = "T2 = A + E"
  )
)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.multi = "index")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
```

## Clone for plant
**[back to top](#contents)** 

***Clone*** is an sexual reproduction method which does not involve germ cells, and does not require a process of fertilization, directly forming a new individual's reproductive mode from a part of the mother. ***Sex*** of offspring will be ***0*** in ```clone```. 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "clone")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Double haploid for plant
**[back to top](#contents)**  

***Double haploid*** is a reproduction method for breeding workers to obtain haploid plants. It induced to double the number of chromosomes and restore the number of chromosomes in normal plants. ***Sex*** of offspring will be ***0*** in ```dh```. 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "dh")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Self pollination for plant and micro-organism
**[back to top](#contents)** 

***Self pollination*** refers to the combination of male and female gametes from the same individual or between individuals from the same clonal breeding line. ***Sex*** of offspring will be ***0*** in ```selfpol```. 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "selfpol")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Random mating for plant and animal
**[back to top](#contents)**  

In ***random mating***, any female or male individual have the same probability to mate with any opposite sex in a sexually reproducing organism. ***Sex*** of offspring in random mating is controlled by ```sex.ratio``` in ```randmate```. 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "randmate")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Random mating excluding self pollination for animal
**[back to top](#contents)**  

In ***random mating excluding self pollination***, an individual cannot mate to itself. ***Sex*** of offspring in random mating is controlled by ```sex.ratio``` in ```randexself```. 

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "randexself")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Two-way cross for animal
**[back to top](#contents)**  

***Two-way cross*** method needs to use ***sex*** to distinguish ***two*** different breeds, in which the ***first breed*** is ***sire*** and the ***second breed*** is ***dam***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "2waycro")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Two different breeds are cut by sex
SP$pheno$pop$gen1$sex <- rep(c(1, 2), c(50, 50))
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Three-way cross for animal
**[back to top](#contents)** 

***Three-way cross*** method needs to use ***sex*** to distinguish ***three*** different breeds, in which the ***first breed*** is ***sire*** and the ***second breed*** is ***dam*** in the ***first two-way cross***, the ***third breed*** is termimal ***sire***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "3waycro")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Three different breeds are cut by sex
SP$pheno$pop$gen1$sex <- rep(c(1, 2, 1), c(30, 30, 40))
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Four-way cross for animal
**[back to top](#contents)**  

***Four-way cross*** method needs to use ***sex*** to distinguish ***four*** different breeds, in which the ***first breed*** is ***sire*** and the ***second breed*** is ***dam*** in the ***first two-way cross***, the ***third breed*** is ***sire*** and the ***fourth breed*** is ***dam*** in the ***second two-way cross***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "4waycro")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Three different breeds are cut by sex
SP$pheno$pop$gen1$sex <- rep(c(1, 2, 1, 2), c(25, 25, 25, 25))
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## Back cross for animal
**[back to top](#contents)**  

***Back cross*** method needs to use ***sex*** to distinguish ***two*** different breeds, in which the ***first breed*** is always ***sire*** in each generation and the ***second breed*** is ***dam*** in the ***first two-way cross***.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate selection parameters
SP <- param.sel(SP = SP, sel.single = "comb")
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "backcro")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Two different breeds are cut by sex
SP$pheno$pop$gen1$sex <- rep(c(1, 2), c(50, 50))
# Run selection
SP <- selects(SP)
# Run reproduction
SP <- reproduces(SP)
```

## User designed pedigree mating for plant and animal
**[back to top](#contents)**  

***User designed pedigree mating*** needs a specific ***user designed pedigree*** to control mating process. The first column is ***sample id***, the second column is ***paternal id***, and the third column is ***maternal id***. Please make sure that ***paternal id*** and ***maternal id*** can match to genotype data.

```r
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)
# Generate genotype simulation parameters
SP <- param.geno(SP = SP, pop.marker = 1e4, pop.ind = 1e2)
# Generate phenotype simulation parameters
SP <- param.pheno(SP = SP, pop.ind = 100)
# Generate reproduction parameters
SP <- param.reprod(SP = SP, reprod.way = "userped")

# Run annotation simulation
SP <- annotation(SP)
# Run genotype simulation
SP <- genotype(SP)
# Run phenotype simulation
SP <- phenotype(SP)
# Run reproduction
SP <- reproduces(SP)
```

## AN EASY WAY TO GENERATE A POPULATION
**[back to top](#contents)**  

The above methods are to generate population ***step by step***, which are ***easy to understand***. Actually, **```SIMER```** can directly generate a population by a ***MORE CONVENIENT*** way. 

```r
# Generate all simulation parameters
SP <- param.simer(qtn.num = 10, pop.marker = 1e4, pop.ind = 1e2, sel.single = "comb", reprod.way = "randmate")

# Run Simer
SP <- simer(SP)
```

---

# Breeding Program Design
**[back to top](#contents)**  

After generating a population, further work can be done. Breeders wish to evaluate their ***Breeding Program Design***. To save a lot of money and time, **```SIMER```** can assist breeders to evaluate their ***Breeding Program Design*** by simulation. 

## Gallery of breeding program design parameters
**[back to top](#contents)**

```simer.Data.Json```, main function of ***Breeding Program Design***:

<table text-algn="left">
<thead>
<tr>
<td><i><b>Paramater</b></i></td>
<td><i><b>Default</b></i></td>
<td><i><b>Options</b></i></td>
<td><i><b>Description</b></i></td>
</tr>
</thead>
<tbody>
<tr>
<td><b>jsonFile</b></td>
<td>NULL</td>
<td>character</td>
<td>the path of JSON file.</td>
</tr>
<tr>
<td><b>out</b></td>
<td>'simer.qc'</td>
<td>character</td>
<td>the prefix of output files.</td>
</tr>
<tr>
<td><b>dataQC</b></td>
<td>TRUE</td>
<td>TRUE or FALSE</td>
<td>whether to make data quality control.</td>
</tr>
<tr>
<td><b>buildModel</b></td>
<td>TRUE</td>
<td>TRUR or FALSE</td>
<td>whether to build EBV model.</td>
</tr>
<tr>
<td><b>buildIndex</b></td>
<td>TRUE</td>
<td>TRUR or FALSE</td>
<td>whether to build Selection Index.</td>
</tr>
<tr>
<td><b>ncpus</b></td>
<td>10</td>
<td>num</td>
<td>the number of threads used, if NULL, (logical core number - 1) is automatically used.</td>
</tr>
<tr>
<td><b>verbose</b></td>
<td>TRUE</td>
<td>TRUE or FALSE</td>
<td>whether to print detail.</td>
</tr>
</tbody>
</table>

## Breeding program design preparation
**[back to top](#contents)**

***Breeding program design*** should be stored on a JSON file. 
> ***plan1.json***  
>>  **genotype**: the path of genotype data  
>>  **pedigree**: the filename of pedigree data  
>>  **selection_index**: the economic weight of phenotype for each trait  
>>  **breeding_value_index**: the economic weight of breeding value for each trait
>>  **quality_control_plan**: the quality control plan for genotype, pedigree, and phenotype  

```json
{
    "genotype": ["/home/yindong/R/x86_64-pc-linux-gnu-library/4.0/simer/extdata/02plinkb"],
    "pedigree": ["/home/yindong/R/x86_64-pc-linux-gnu-library/4.0/simer/extdata/05others/pedigree.txt"],
    "selection_index": [],
    "breeding_value_index": "0.2 * T1 + 0.8 * T2",
    "quality_control_plan": {
        "genotype_quality_control":{
            "filter": ["F1 == 'Male'"],
            "filter_geno": 0.1,
            "filter_mind": 0.1,
            "filter_maf": 0.05,
            "filter_hwe": 0.001
        },
        "pedigree_quality_control":{
            "standard_ID": false,
            "candidate_sire_file": [],
            "candidate_dam_file": [],
            "exclude_threshold": 0.01, 
            "assign_threshold": 0.005
        },
        "phenotype_quality_control":[
            {
                "job_name": "Data Quality Control Demo",
                "sample_info": "/home/yindong/R/x86_64-pc-linux-gnu-library/4.0/simer/extdata/05others/phenotype.txt",
                "repeated_records": false,
                "multi_trait": true,
                "filter": ["F1 == 'Male'"],
                "job_traits": [
                    {
                        "traits": "T1",
                        "definition": "T1",
                        "range": []
                    },
                    {
                        "traits": "T2",
                        "definition": "T2",
                        "range": []
                    }
                ]
            }
        ]
    },
    "analysis_plan":[
        {
            "job_name": "EBV Model Demo",
            "sample_info": "/home/yindong/R/x86_64-pc-linux-gnu-library/4.0/simer/extdata/05others/phenotype.txt",
            "repeated_records": false,
            "multi_trait": true,
            "random_ratio": 0.05,
            "job_traits": [
                {
                    "traits": "T1",
                    "covariates": [],
                    "fixed_effects": ["F1", "F2"],
                    "random_effects": ["R1"]
                },
                {
                    "traits": "T2",
                    "covariates": [],
                    "fixed_effects": ["F1", "F2"],
                    "random_effects": ["R1"]
                }
            ]
        }
    ]
}
```

## Breeding plan comparison
**[back to top](#contents)**  

After a total reproduction process, breeding plans comparison can be done. 

```r
# random-mating
simer.list <-
    simer(num.gen = 10,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          # num.ind = 100,
          mtd.reprod = "randmate",
          num.prog = 2,
          ratio = 0.5)

pop <- simer.list$pop
out.pop <- read.selgeno(pop = pop, selPath = selPath, outpath = outpath)

# make comparison on breeding plans
plan1 <- out.pop$plan1
plan2 <- out.pop$plan2
plan3 <- out.pop$plan3
summary(plan1)
summary(plan2)
summary(plan3)
```

---

# Global Options
**[back to top](#contents)**  

In this part, calculation of population size and different ourput methods will be introduced. 

## Gallery of global input parameters
**[back to top](#contents)** 

`simer()`, main function:  
**num.gen**, number of generations in simulation  
**replication**, replication index of simulation  
**verbose**, whether to print detail  
**out**, prefix of output file name  
**outpath**, path of output files  
**out.format**, format of output, "numeric" or "plink"  
**seed.sim**, random seed of a simulation process  
**seed.map**, random seed of map file  
**out.geno.gen**, indice of generations of output genotype  
**out.pheno.gen**, indice of generations of output phenotype  

## Counts of total population size
**[back to top](#contents)**  

The following is the method of obtaining population size of every generation. Every elements in **cound.ind** are population size in this generation respectively. 

```r
# parameters that controls population size of every generations 
nind <- 40
num.gen <- 10
ratio <- 0.5
sel.on <- TRUE # whether selection exsits
ps <- 0.8
ps <- ifelse(sel.on, ps, 1)
num.prog <- 2
mtd.reprod <- "randmate"

# populations of the first generation
basepop <- getpop(nind = nind, from = 1, ratio = ratio)
pop2 <- getpop(nind = nind, from = 41, ratio = ratio)
pop3 <- getpop(nind = nind, from = 81, ratio = ratio)
pop4 <- getpop(nind = nind, from = 121, ratio = ratio)

# calculate number of individuals in every generation
count.ind <- rep(nind, num.gen)
if (mtd.reprod == "clone" || mtd.reprod == "dh" || mtd.reprod == "selfpol") {
  if (num.gen > 1) {
    for(i in 2:num.gen) {
      count.ind[i] <- round(count.ind[i-1] * (1-ratio) * ps) * num.prog
    }
  }
  
} else if (mtd.reprod == "randmate" || mtd.reprod == "randexself") {
  if (num.gen > 1) {
    for(i in 2:num.gen) {
      count.ind[i] <- round(count.ind[i-1] * (1-ratio) * ps) * num.prog
    }
  }

} else if (mtd.reprod == "singcro") {
  sing.ind <- round(nrow(pop2) * ps) * num.prog
  count.ind <- c(nrow(basepop), nrow(pop2), sing.ind)

} else if (mtd.reprod == "tricro") {
  dam21.ind <- round(nrow(pop2) * ps) * prog.tri
  tri.ind <- round(dam21.ind * (1-ratio) * ps) * num.prog
  count.ind <- c(nrow(basepop), nrow(pop2), nrow(pop3), dam21.ind, tri.ind)

} else if (mtd.reprod == "doubcro") {
  sir11.ind <- round(nrow(pop2) * ps) * prog.doub
  dam22.ind <- round(nrow(pop4) * ps) * prog.doub
  doub.ind <- round(dam22.ind * (1-ratio) * ps) * num.prog
  count.ind <- c(nrow(basepop), nrow(pop2), nrow(pop3), nrow(pop4), sir11.ind, dam22.ind, doub.ind)

} else if (mtd.reprod == "backcro") {
  count.ind[1] <- nrow(basepop) + nrow(pop2)
  if (num.gen > 1) {
    count.ind[2] <- round(nrow(pop2) * ps) * num.prog
    for(i in 3:num.gen) {
      count.ind[i] <- round(count.ind[i-1] * (1-ratio) * ps) * num.prog
    }
  }
}
```

## Simulation of multiple populations
**[back to top](#contents)**

Simulation of multiple populations can be realized by "for" in **R** and **replication**. Random seed of simulation is random in every replication and random seed of map is fixed. In every replication, you can set your own random seed of simulation and random seed of map by **seed.geno**| and **seed.map** respectively. 

```r
# random-mating 
rep <- 2
for (i in 1:rep) {
  simer.list <-
    simer(num.gen = 3,
          replication = i, # set index of replication
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 100,
          mtd.reprod = "randmate",
          num.prog = 2,
          ratio = 0.5)
}
```

## File output
**[back to top](#contents)** 

**SIMER** won't output files by default. A series of files with prefix `out = simer` output when specify a exsited out path by **outpath**. File output format can be "numeric" or "plink" by **out.format**. 

```r
# set a output path
outpath = getwd()

# "numeric" format
simer.list <-
  simer(num.gen = 3,
        verbose = verbose, 
        outpath = outpath,
        out.format = "numeric", 
        input.map = input.map,
        # rawgeno1 = rawgeno, # use your own genotype matrix
        num.ind = 100,
        mtd.reprod = "randmate",
        num.prog = 2,
        ratio = 0.5)

# "plink" format
simer.list <-
  simer(num.gen = 3,
        verbose = verbose, 
        outpath = outpath,
        out.format = "plink", 
        input.map = input.map,
        # rawgeno1 = rawgeno, # use your own genotype matrix
        num.ind = 100,
        mtd.reprod = "randmate",
        num.prog = 2,
        ratio = 0.5)
```

## Generation selective output
**[back to top](#contents)**  

Output of genotype and phenotype can be generation-selective by **out.geno.gen** and **out.pheno.gen**. For example, **out.geno.gen = 3:5** and **out.pheno.gen = 1:5** represent outputting genotype from generation3 to generation5 and outputting phenotype from generation1 to generation5.

```r
# set output generations of genotype and phenotype
out.geno.gen <- 3:5
out.pheno.gen <- 1:5
simer.list <-
  simer(num.gen = 5,
        verbose = verbose, 
        outpath = outpath,
        out.geno.gen = out.geno.gen, 
        out.pheno.gen = out.pheno.gen, 
        input.map = input.map,
        # rawgeno1 = rawgeno, # use your own genotype matrix
        num.ind = 100,
        mtd.reprod = "randmate",
        num.prog = 2,
        ratio = 0.5)
```

---

# Output
**[back to top](#contents)**  

**SIMER** outputs data including population information, marker effects, trait information, genotype, genotypic id, genotypic map, and selection intensity. Note that several result files with prefix `out = "simer"` will be generated. Firstly, a result path will be generated. The number at the beginning represents the total individuals number and the ending character is the format of output ("num_Simer_Data_format"). Secondly, you will see a path named "replication1". It is the first replication of simulation and you can get numerous different replications by setting parameter **replication**. In "replication1", results are following:  
`simer.geno.id` contains indice of genotyped individuals  
`simer.geno.bin` and `simer.geno.desc` contain genotype matrix of all individuals  
`simer.map` contains input map with block information and recombination information  
`simer.ped` contains pedigree of individuals  
`simer.phe` contains phenotype of individuals  

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
