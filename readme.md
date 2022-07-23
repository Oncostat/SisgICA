


# SisgICA: Survival - infinite sparse graphical Independent Component Analysis (Variational Inference)


These R codes are the online material for the article:
  
A nonparametric Bayesian joint model for latent individual molecular profiles and survival in oncology. \
By Sarah-Laure Rincourt, Stefan Michiels and Damien Drubay

Correspondance author : `damien.drubay@gustaveroussy.fr`

## Purpose 

Identifying the molecular biomarkers that are associated with a patient's prognosis is essential in stratified medicine. Even though several individual biomarkers (e.g., mutations, copy number alterations or gene expression values) with strong effects on the survival of cancer patients have been identified in the past, the synergy of the weaker statistical signals of the different biomarkers involved in intricate molecular pathways is more complex to identify.


Non-observed structures in omics data are usually identified using latent variable models for dimension reduction. In oncology, different tumors of the same organ can result from various molecular mechanisms in different patients, which complicates the development of stratified or targeted treatment.

To capture the individual patient's heterogeneity, we proposed a joint model of survival analysis and the proposed infinite sparse graphical ICA (isgICA), which is referred to SisgICA. We assumed that tumor expression results from a mixture of a subset of independent signatures. We deconvoluted the omics data using a non-parametric independent component analysis with a double sparseness structure for the source and the weight matrices, corresponding to the gene-component and individual-component associations, respectively.

The proposed algorithm provides a new insight into the individual molecular heterogeneity that is associated with patient prognosis to better understand the complex tumor mechanisms.


We present two folders: one simulation study and one application in early breast cancer.


## How to use ?

### Prerequisites


The R setup is available on the CRAN website (https://cran.r-project.org/).

The following extra-packages are required:
  

```
library(zeallot)  
library(ggplot2)
library(gridExtra)
library(reshape2)
library(GPfit) 
library(ParBayesianOptimization)
library(doParallel)
library(compiler)
library(class)
library(RcppHungarian)
library(psych)
library(survival)
library(base)
library(mvtnorm)
library(invgamma)
library(stats) 
library("profvis")
library(tensorflow)
library(tfprobability)
```




For installing tensorflow


```
devtools::install_github("rstudio/tensorflow")
library(tensorflow)
install_tensorflow()
```


For installing tfprobability


```
devtools::install_github("rstudio/tfprobability")
library(tfprobability)
install_tfprobability()
```




### Run  

#### Pre-processing

The isgICA model takes in input normalized and center covariables (in row the covariables and in column the individuals).

#### Example

The code of the SisgICA model are in the `00._SisgICA_model.r` file. The `data` used is in the Breast folder and the simulation study is in the simulation folder. 
To launch the codes, simply set the path in the launch files. The object `X` is the expression data with in columns and in rows, the individuals and the genes respectively.


```
setwd("...")
tryCatch(dir.create("result"))
tryCatch(dir.create("./result/inprogress"))

initPar = list(y=data$time,
               d=sum(data$status),delta=data$status,
               X=X, N=dim(X)[2], P=dim(X)[1],
               Kmax=100L,
               maxit = 1000,
               seed=1, 
               a_prior=1,
               cc=1,
               dd=1,
               e=1e-6,
               f=1e-6,
               vec_seed = rep(1,2),
               value_var_error=NULL,
               type=NULL, 
               typeZ=NULL,
               sim = NULL,
               namefile = "result",
               ic_freq2 = NULL,
               paraCI = NULL)

source("00._SisgICA_model.r")
source("00._GP_additionnal_function.R")
fit = Do_Opt(initPar)

```

To observe to results, we can run :

```
source("02._plot_model.r")

# Figure: isgICA sub-model
grid.arrange(p2,p3, p4, ncol=3)  

# Figure: survival sub-model
p5

# Table: survival sub-model coefficient
list_recap2
```









