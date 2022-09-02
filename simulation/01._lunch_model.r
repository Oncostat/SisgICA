# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

cat("\014") # clear console
rm(list=ls())
graphics.off()

#####################################################################
########################### library #################################

library(zeallot) # to use :  %<-%
library(glue)
library(GPfit)
library(magrittr) # for %>%
library(ggplot2)
library(gridExtra)
library(reshape2) #melt 
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
library(CholWishart) # dWishart
library("profvis")

#####################################################################
################# Path to the code directory ########################


# setwd("...")

tryCatch(dir.create("result"))
tryCatch(dir.create("./result/inprogress"))

#####################################################################
###################### Initialization ###############################


source("00._SisgICA_model.r")
source("00._additionnal_function.r")


###################### Init data ###################### 

source("00._simulation_SisgICA.r")

simu_SisgICA_save()
load("simu_N_600_P_1200_K_10_censured_0.3_var_error_1_seed_100_typeZ_random_Z_0.314_Zphi_0.264666666666667.rdata")

seed = 1

X = sim$X
P = dim(X)[1]
N = dim(X)[2]
Kmax = 100L
Kmax = min(Kmax,dim(X))
seed_val = 1
maxit = 1000

initPar = list(y=sim$time,d=sum(sim$delta),delta=sim$delta,
               X=X, N=N, P=P,Kmax=Kmax,
               maxit = maxit ,
               seed=seed_val, 
               a_prior=1,
               cc=1,
               dd=1,
               e=1e-6,
               f=1e-6,
               vec_seed = rep(seed_val,2),
               value_var_error=NULL,
               type=NULL, 
               typeZ=NULL,
               sim = NULL,
               namefile = "Breast_filtration",
               ic_freq2 = NULL,
               paraCI = NULL)

#####################################################################
###################### Analyse function #############################

print("GO1")


analyse2 <- function(initPar){
  source("00._SisgICA_model.r")
  source("00._GP_additionnal_function.R")
  print(c(dim(initPar$X), initPar$vec_seed[1]))
  fit = Do_Opt(initPar)
  return(fit)
}

# Execution of the program
print("GO2")
start <- Sys.time()

z = analyse2(initPar)

end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
print(t.parallel)



source("02._plot_model.r")


# Figure isgICA sub-model
grid.arrange(p2,p3, p4, ncol=3)  


# coefficient survival sub-model
list_recap2




