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
    - [Generate base population information](#generate-base-population-information)  
    - [Generate phenotype of single trait by A model](#generate-phenotype-of-single-trait-by-A-model)  
    - [Generate phenotype of single trait by AD model](#generate-phenotype-of-single-trait-by-AD-model)  
    - [Generate phenotype of single trait by ADI model](#generate-phenotype-of-single-trait-by-ADI-model)  
    - [Generate phenotype of multiple traits](#generate-phenotype-of-multiple-traits)  
    - [Add fixed effects and random effects to phenotype](#add-fixed-effects-and-random-effects-to-phenotype)
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
- [Comparison on Breeding Plans](#comparison-on-breeding-plans) 
    - [Gallery of comparison input parameters](#gallery-of-comparison-input-parameters)
    - [Breeding plan preparation](#breeding-plan-preparation)
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
At least user should prepare two datasets: ***genotypic map*** and ***genotype data***.  

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
# Generate annotation simulation parameters
SP <- param.annot(qtn.num = 10)

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

Note that recombination only exists in meiosis. Therefore, some reproduction methods such as ```clone``` do not have recombination process. User can set ```recom.spot = FALSE``` to add only mutations to genotype matrix.  

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
<td>'norm'</td>
<td>character</td>
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

In ***A*** model, **```SIMER```** only considers ***A***dditive effect as genetic effect. User should prepare ***A***dditive ***QTN*** effect in the ***Annotation data*** for generating ***A***dditive ***I***ndividual effect. ***A***dditive single trait simulation is displayed as follows: 

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

In multiple-trait simulation, **```SIMER```** can build ***accurate Additive genetic correlation*** between multiple traits. ***A***dditive multiple trait simulation is displayed as follows: 


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

In single trait simulation, the trait can be controlled by ***multiple-group QTNs***. An example of single trait controlled by two-group QTNs is displayed as follows: 

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

In ***AD*** model, **```SIMER```** considers ***A***dditive effect and ***D***ominant effect as genetic effect. User should prepare ***A***dditive ***QTN*** effect and ***D***ominant ***QTN*** effect in the ***Annotation data*** for generating ***A***dditive ***I***ndividual effect and ***D***ominant ***I***ndividual effect. ***A***dditive and ***D***ominant single trait simulation is displayed as follows: 

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

In multiple-trait simulation, **```SIMER```** can build ***accurate Additive genetic correlation*** and ***accurate Dominant genetic correlation*** between multiple traits. ***A***dditive and ***D***ominant multiple trait simulation is displayed as follows: 


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

In ***GxG*** model, **```SIMER```** considers ***G***enetic-***G***enetic effect as genetic effect. User should prepare ***G***enetic-***G***enetic ***QTN*** effect in the ***Annotation data*** for generating ***G***enetic-***G***enetic ***I***ndividual effect. An example of ***A***dditive-***D***ominant interaction is displayed as follows: 

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

In multiple-trait simulation, **```SIMER```** can build ***accurate Genetic-Genetic interaction correlation*** between multiple traits. An example of ***A***dditive-***D***ominant interaction multiple trait simulation is displayed as follows: 


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

In ***Repeated Record*** model, **```SIMER```** adds ***PE*** (***P***ermanent ***E***nvironmental) effect to the phenotype. The number of repeated records can be set by ```pop.rep`` and ```pop.rep.bal``` can be used to determine whether repeated records are balanced. ***Repeated Record*** single trait simulation is displayed as follows: 

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

In multiple-trait simulation, **```SIMER```** can build ***accurate Permanent Environmental correlation*** between multiple traits. ***Repeated Record*** multiple trait simulation is displayed as follows: 

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

## Add fixed effects and random effects to phenotype
**[back to top](#contents)** 

**SIMER** supports add fixed effects and random effects to phenotype. Users should prepare a list with fixed effects, random effects, and their combination first. In the list, the attribution name are fixed. "cmb.fix" contains different combinations of fixed factors in different traits. "fix" contains different fixed factors which have their levels and effects. "cmb.rand" contains different combinations of random factors. Unlike the fixed effect, the user needs to specify the proportion of the variance of the random effect variance by "ratio". Besides, users can also set related random effects for different traits by "cr". Users can use attribute names in population information, meaning that users can assign specific levels to different individuals. At the same time, if the attribute name is not in the population information, **SIMER** will add its level to the population information. 

```r
# specific levels to different individuals
a <- sample(c("a1", "a2", "a3"), nind, replace = TRUE)
b <- sample(c("b1", "b2", "b3"), nind, replace = TRUE)
basepop1$a <- a # load your fixed  effects
basepop1$b <- b # load your random effects

# combination of fixed effects
cmb.fix <- list(tr1 = c("mu", "gen", "sex", "a"), # trait 1
                tr2 = c("mu", "diet", "season"))  # trait 2		

# available fixed effects
fix <- list(
         mu = list(level = "mu", eff = 2),
        gen = list(level = "1", eff = 1), # inherent fixed effect in simer
        sex = list(level = c("1", "2"), eff = c(5, 3)), # inherent fixed effect in simer
       diet = list(level = c("d1", "d2", "d3"), eff = c(1, 2, 3)),
     season = list(level = c("s1", "s2", "s3", "s4"), eff = c(1, 2, 3, 2)), 
          a = list(level = c("a1", "a2", "a3"), eff = c(1, 2, 3)))

# combination and ralation of random effects
# rn, random effect name
# ratio, phenotype variance proportion of the random effects
# cr, corelation of the random effects
tr1 <- list(rn = c("sir", "dam", "b"), ratio = c(0.03, 0.05, 0.03))
tr2 <- list(rn = c("PE", "litter"), ratio = c(0.01, 0.03),
            cr = matrix(c(1, 0.5, 0.5, 1), 2, 2))
cmb.rand <- list(tr1 = tr1, tr2 = tr2)

# available random effects
rand <- list(
        sir = list(mean = 0, sd = 1), # sir and dam are inherent random effect in simer
        dam = list(mean = 0, sd = 1), # control mean and sd only
         PE = list(level = c("p1", "p2", "p3"), eff = c(1, 2, 3)),
     litter = list(level = c("l1", "l2"), eff = c(1, 2)), 
          b = list(level = c("b1", "b2", "b3"), eff = c(1, 2, 3)))

FR <- list(cmb.fix = cmb.fix, fix = fix, cmb.rand = cmb.rand, rand = rand)
```

After preparation above, you can get phenotype with fixed effects and random effects but without additive effects by setting **pop.geno = NULL**. 

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = NULL, # set pop.geno = NULL
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              FR = FR, # input fixed effects and random effects
              pop = basepop1,
              pop.geno = NULL, # set pop.geno = NULL
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

If inputting **pop.geno** a genotype matrix, you will get phenotype with fixed effects, random effects, and genetic effects. 

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno, # input genotype matrix
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              FR = FR, # input fixed effects and random effects
              pop = basepop1,
              pop.geno = basepop1.geno, # input genotype matrix
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Users can also add fixed effects and random effects to different traits respectively by multiple traits mode. 

```r
#######################
### multiple traits ###
# calculate for marker information
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(2, 6, 10),
             sd.tr1 = c(0.4, 0.2, 0.02, 0.02, 0.02, 0.02, 0.02, 0.001),
             dist.qtn.tr1 = rep("normal", 6),
             eff.unit.tr1 = rep(0.5, 6),
             shape.tr1 = rep(1, 6),
             scale.tr1 = rep(1, 6),
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = matrix(c(1, 0, 0, 2), 2, 2),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop2.pheno <-
     phenotype(effs = effs,
               FR = FR, 
               pop = basepop2,
               pop.geno = basepop2.geno,
               pos.map = NULL,
               h2.tr1 = c(0.4, 0.2, 0.1, 0.1, 0.1, 0.05),
               gnt.cov = matrix(c(14, 10, 10, 15), 2, 2),
               h2.trn = c(0.3, 0.5),
               sel.crit = "pheno",
               pop.total = basepop2,
               sel.on = TRUE, 
               inner.env = NULL,
               verbose = verbose)

# get population with phenotype
basepop2 <- pop2.pheno$pop
pop2.pheno$pop <- NULL
```

## Different QTN effect distributions
**[back to top](#contents)**  

In different model, you can further set different QTN effect distributions of trait1 by **dist.qtn.tr1*|**. The most common distribution is "normal" distribution. You can set different variances in "normal" distribution by **dist.qtn.tr1**.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

QTN effect distribution can be "geometry" distribution. You can set effect unit of "geometry" by **eff.unit.tr1**.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # a
             dist.qtn.tr1 = "geometry",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1, # effect unit of geomtry distribution
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

"Gamma" distribution is also a kind of QTN effect distribution. You can set shape and scale of "gamma" distribution by **shape.tr1** and **scale.tr1**. Note that default options of "gamma" distribution **shape.tr1 = 1** and **scale.tr1 = 1** exactly lead to exponential distribution.

```r
####################
### single trait ###
# calculate for marker information
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # a
             dist.qtn.tr1 = "gamma",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1, # shape of gamma distribution
             scale.tr1 = 1, # scale of gamma distribution
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

## Different selection criteria
**[back to top](#contents)**  

In addition, "pheno" is not just phenotype. It can also be "TBV", "TGV", "pEBVs'", "gEBVs",  or "ssEBVs". "TBV" is True Breeding Value, represents only that part of genotypic value that can be transmitted from parent to offspring. 

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TBV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "TBV",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "TBV".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TBV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "TBV",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

"TGV" is True Genotypic Value, represents the sum of additive effect, dominance effect and epistatic effect.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TGV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "TGV",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "TGV".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "TGV"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "TGV",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

"pEBVs" is pedigree Estimated Breeding Values. It means that BLUP constructs kinship by pedigree. You can get "pEBVs" by "ABLUP" model in HIBLUP. 

```r
# call "hiblup" package
suppressMessages(library("hiblup"))

####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
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
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "pEBVs".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

"gEBVs" is genomic Estimated Breeding Values. It means that BLUP constructs kinship by genotype matrix. You can get "gEBVs" by "GBLUP" model in HIBLUP.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
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
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "gEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "gEBVs".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "gEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "gEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

"ssEBVs" is single-step genomic Estimated Breeding Values. It means that BLUP constructs kinship by pedigree and genotype matrix. You can get "ssEBVs" by "SSBLUP" model in HIBLUP.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
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
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "ssEBVs",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "ssEBVs".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "ssEBVs"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "ssEBVs",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

At last, "pheno" is phenotype including additive effect (and dominance effect) (and epistatic effect) and residual effect.

```r
####################
### single trait ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pheno"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop1,
              pop.geno = basepop1.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop1,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

Phenotype of multiple traits can also be represented as "pheno".

```r
#######################
### multiple traits ###
# Additive model
effs <-
    cal.effs(pop.geno = basepop2.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = 18,
             sd.tr1 = 0.6, # standard deviation of normal distribution
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 0.5,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = TRUE, # multiple traits
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
# selection criterion is "pheno"
pop.pheno <-
    phenotype(effs = effs,
              pop = basepop2,
              pop.geno = basepop2.geno,
              pos.map = NULL,
              h2.tr1 = 0.8,
              gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
              h2.trn = c(0.3, 0.5),
              sel.crit = "pheno",
              pop.total = basepop2,
              sel.on = TRUE,
              inner.env = NULL,
              verbose = verbose)

# get population with phenotype
basepop2 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

## Multiple groups QTN effects
**[back to top](#contents)**  

Multiple groups QTN effects can be realized by setting different elements of **num.qtn.tr1**, every elements represent amount of QTNs affacting a effect. For example, **num.qtn.tr1 = c(2, 6, 10)** means that the additive effect of the trait is the sum of three QTN group effects. The first group has 2 QTNs, the second group has 6 QTNs and the third group has 10 QTNs. Because of the three groups, the first three elements of **sd.tr1** mean the variances of the three QTN groups.

```r
#####################
### single traits ###
# calculate for marker information
effs <-
    cal.effs(pop.geno = basepop1.geno,
             cal.model = "A", # it can be"A", "AD" or "ADI"
             num.qtn.tr1 = c(2, 6, 10),
             sd.tr1 = c(0.4, 0.2, 0.02),
             dist.qtn.tr1 = "normal",
             eff.unit.tr1 = 1,
             shape.tr1 = 1,
             scale.tr1 = 1,
             multrait = FALSE, # single trait
             num.qtn.trn = matrix(c(18, 10, 10, 20), 2, 2),
             sd.trn = diag(c(1, 0.5)),
             qtn.spot = rep(0.1, 10),
             maf = 0,
             verbose = verbose)

# generate phenotype
# generate single trait or multiple traits according to effs
pop.pheno <-
     phenotype(effs = effs,
               pop = basepop1,
               pop.geno = basepop1.geno,
               pos.map = NULL,
               h2.tr1 = 0.3,
               gnt.cov = matrix(c(1, 2, 2, 15), 2, 2),
               h2.trn = c(0.3, 0.5),
               sel.crit = "pheno",
               pop.total = basepop1,
               sel.on = TRUE,
               inner.env = NULL,
               verbose = verbose)

# get population with phenotype
basepop1 <- pop.pheno$pop
pop.pheno$pop <- NULL
```

---

# Selection
**[back to top](#contents)**  

You can get ordered individuals indice according to phenotype in the populaton information. Fraction selected can be used to keep a certain amount of individuals. SIMER chooses automatically single trait selection or multiple traits selection according to number of columns of phenotype.

## Gallery of selection input parameters
**[back to top](#contents)**  

`selects()`, main function of selection:  
**pop**, population information of generation, family index, within-family index, index, sire, dam, sex, phenotpye  
**decr**, whether to sort by descreasing  
**sel.multi**, selection method of multiple traits with options: "tdm", "indcul" and "index"  
**index.wt**, economic weights of selection index method, its length should equals to the number of traits  
**index.tdm**, index represents which trait is being selected  
**goal.perc**, percentage of goal more than mean of scores of individuals  
**pass.perc**, percentage of expected excellent individuals  
**sel.sing**, selection method of single trait with options: "ind", "fam", "infam" and "comb"  
**pop.total**, total population infarmation  
**pop.pheno**, list of all phenotype information  
**verbose**, whether to print detail  

`simer()`, main function:  
**ps**, fraction selected in selection  

## Individual selection on single trait
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

## Family selection on single trait
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

## Within family selection on single trait
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

## Combined selection on single trait
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

## Tandem selection on multiple traits
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

## Independent culling selection on multiple traits
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

## Index selection on multiple traits
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

## Gallery of reproduction input parameters
**[back to top](#contents)** 

`simer()`, main function:  
**mtd.reprod**, different reproduction methods with options: "clone", "dh", "selfpol", "singcro", "tricro", "doubcro", "backcro","randmate", "randexself" and "userped"  
**userped**, user-designed pedigree to control mating process  
**num.prog**, litter size of dams  
**ratio**, ratio of males in all individuals  
**prog.tri**, litter size of the first single cross process in trible cross process  
**prog.doub**, litter size of the first two single cross process in double cross process  
**prog.back**, a vector with litter size in every generation of back-cross  

## Clone
**[back to top](#contents)** 

Asexual reproduction does not involve germ cells, and does not require a process of fertilization, directly forming a new individual's reproductive mode from a part of the mother. Sex of offspring will be "0" in clone. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# clone
simer.list <-
    simer(num.gen = 3,
          verbose = verbose, 
          outpath = outpath, 
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "clone",
          num.prog = 2,
          ratio = 0.5)
```

## Double haploid
**[back to top](#contents)**  

Breeding workers often use another culture in vitro to obtain haploid plants, and then artificially induced to double the number of chromosomes and restore the number of chromosomes in normal plants. This method is named double-haploid reproduction. Sex of offspring will be "0" in double-haploid. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# double-haploid
simer.list <-
    simer(num.gen = 4,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "dh",
          num.prog = 2,
          ratio = 0.5)
```

## Self pollination
**[back to top](#contents)** 

Self-pollination refers to the combination of male and female gametes from the same individual or between individuals from the same clonal breeding line. Sex of offspring will be "0" in self-pollination. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.  

```r
# self-pollination
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 50,
          mtd.reprod = "selfpol",
          num.prog = 2,
          ratio = 0.5)
```

## Random mating
**[back to top](#contents)**  

In random mating, any female or male individual have the same probability to mate with any opposite sex in a sexually reproducing organism. Sex of offspring in random mating is up to sex of parents. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# random-mating
simer.list <-
    simer(num.gen = 4,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 100,
          mtd.reprod = "randmate",
          num.prog = 2,
          ratio = 0.5)
```

## Random mating without self pollination
**[back to top](#contents)**  

In random mating without self-pollination, a individual cannot mate to itself. Sex of offspring in random mating without self-pollination is up to sex of parents. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# random-mating without self-pollination
simer.list <-
    simer(num.gen = 3,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          # rawgeno1 = rawgeno, # use your own genotype matrix
          num.ind = 100,
          mtd.reprod = "randexself",
          num.prog = 2,
          ratio = 0.5)
```

## User designed pedigree mating
**[back to top](#contents)**  

User-designed-pedigree mating needs a specific user designed pedigree to control mating process. Pedigree should at least start with generation 2. Please make sure that paternal id and maternal id can be found in the last generation. Note that the individuals in the pedigree data file do not need to be sorted by the date of birth, and the missing value can be replaced by NA or 0. Single-breed methods needs only one genotype matrix, you can generate a random genotype matrix by setting **num.prog** or use your own genotype matrix by setting **rawgeno1 = your_own_matrix**.

```r
# user-designed-pedigree mating
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          # num.ind = 100,
          mtd.reprod = "userped",
          userped = userped) # input your own pedigree
```

## Two way cross
**[back to top](#contents)**  

Two-way cross method needs two genotype matrice of different two breeds. You can input your own genotype matrix by parameters **rawgeno1** and **rawgeno2**. If any of these two is NULL, **SIMER** will generates a random one.

```r
# two-way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
          input.map = input.map,
          rawgeno1 = rawgeno, # use your own genotype matrix
          rawgeno2 = NULL,
          # num.ind = 100,
          mtd.reprod = "singcro",
          num.prog = 2,
          ratio = 0.5)
```

## Three way cross
**[back to top](#contents)** 

Three-way cross method needs three genotype matrice of different three breeds. You can input your own genotype matrix by parameters **rawgeno1**, **rawgeno2**, and **rawgeno3**. If any of these three is NULL, **SIMER** will generates a random one. In triple-cross method, you can set litter size of the single-cross of population2 and population3 by **prog.tri**.

```r
# three way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
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

## Four way cross
**[back to top](#contents)**  

Four-way cross method needs four genotype matrice of different four breeds. You can input your own genotype matrix by parameters **rawgeno1**, **rawgeno2**, **rawgeno3**, and **rawgeno4**. If any of these four is NULL, **SIMER** will generates a random one. In four-way cross method, you can set litter size of the first two two-way cross by **prog.doub**.

```r
# four way cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
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

## Back cross
**[back to top](#contents)**  

Back-cross method needs two different breeds. You can input your own genotype matrix by parameters **rawgeno1** and **rawgeno2**. If any of these two is NULL, **SIMER** will generates a random one. Back-cross is similar to two-way cross but with some differences: 1. Back-cross is multi-generation mating; 2. The first base population is fixed in every generations. In back-cross method, you can set litter size of two-way cross in every generation by a vector **prog.back**.

```r
# Back-cross
simer.list <-
    simer(num.gen = 5,
          verbose = verbose, 
          outpath = outpath,
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

# Comparison on Breeding Plans
**[back to top](#contents)**  

After a total reproduction process, further work can be done. For breeders, they want to predict their breeding plans. To save a lot of money and time, simulation will assist them to make comparison on different breeding plans. 

## Gallery of comparison input parameters
**[back to top](#contents)**

`read.selgeno()`, function to make comparison on breeding plans:  
**pop**, total population information  
**selPath**, the path of breeding plans  
**outpath**, path of output files  

## Breeding plan preparation
**[back to top](#contents)**

Breeding plans should be stored on different files respectively. Filenames must begin with breeding_plan, such like breeding_plan01.txt. For now, **SIMER** supports different breeding plans on generation, family_index, within_family_index, and sex. 

```r
> # get path of breeding plan
> selPath <- system.file("extdata", "01breeding_plan", package = "simer")

> # show filename format
> dir(selPath)
[1] "breeding_plan01.txt" "breeding_plan02.txt" "breeding_plan03.txt" "ReadMe.txt" 

> # show file contents format
> fn1 <- file.path(selPath, "breeding_plan01.txt")
> bp1 <- read.delim(fn1, header = TRUE)

> # ###NOTE###
> # 1. any generation can be choosed within num.gen.
> # 2. any sex(1 represents sir 2 represents dam) can be choosed.
> # 3. any family index can be choosed in all family. 1:5 means 
> #    that the first 5 family will be chosen.
> # 4. any infamily index can be choosed in every family.
> # 5. numbers should be separated by comma(",").
> # 6. ":" can choose serial number, "1:3,7" represents "1 2 3 7".
> # 7. "all" will choose all of options in a category.
> bp1
  generation family_index within_family_index sex
1          1        1:3,7                 all   1
2          2        1:3,7                 all   1
3          3        1:3,7                 all   1
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
