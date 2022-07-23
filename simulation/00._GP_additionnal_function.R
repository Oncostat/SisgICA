getTestFun <- function(initPar) {
  fBO=function(U){
    mu = qbeta(U,0.25,2.5)
    b_prior=1/mu-1
    initPar$b_prior = b_prior
    y<<-initPar$y;delta <<- initPar$delta;d<<- initPar$d
    fit =  Weib_isgICA_model(y = initPar$y,delta = initPar$delta,d = initPar$d,
                             initPar$X,Kmax=as.integer(initPar$Kmax),
                             maxit = initPar$maxit,
                             seed=initPar$seed,
                             iniVarW = 1, iniVarPhi = 1,a=1,b=initPar$b_prior,
                             initPar$cc,initPar$dd,initPar$e,initPar$f)
    
    parmsBase = list(rateC = initPar$sim$censth,
                     var_value = initPar$value_var_error,
                     b_prior = initPar$b_prior,
                     seed = initPar$vec_seed[1],
                     Kth = length(initPar$sim$coef_reg))
    
    save(fit,initPar,
         parmsBase,
         file = paste0("result/_N_",initPar$N,
                       "_P_",initPar$P,"_a_prior_",initPar$a_prior, 
                       "_b_prior_",initPar$b_prior,
                        "_ite",initPar$maxit,  ".rdata"    ))
         
    return(list(Score=fit$lastElbo))
  }
  fBO
}



Do_Opt=function(initPar) { 
  Test_Fun <- getTestFun(initPar)
  beta = rev(c(1,10,100,1000,10000,100000) )
  mu = 1 / (beta +1)
  U = pbeta(mu, 0.25,2.5) 
  bornes = pbeta(1 / (c(1e6,1e-6)+1), 0.25,2.5)
  
  tryCatch( 
    Results <- bayesOpt(FUN = Test_Fun, bounds = list(U=bornes),  
                        initGrid = data.frame(U), 
                        iters.n = 9,iters.k=1, 
                        gsPoints = 1000,parallel=F,errorHandling="continue"),
    error = function(e){
      write(as.character(paste0("vec_seed=", initPar$vec_seed ,
                                "_K_",length(initPar$sim$coef_reg),
                                "_delta_",initPar$sim$censth,
                                " value_var_error=",
                                initPar$value_var_error,
                                # " typeZ=", initPar$typeZ,
                                " Erreur: stop bayesOpt")),
            "test2.txt", append=TRUE)
      write(as.character(e), "test2.txt", append=TRUE)
    })
  if(exists("Results")){
    save(Results,initPar, file=paste0(
      "result/result_all","_N_",initPar$N,
      "_P_",initPar$P, 
      "_ite_",initPar$maxit, ".rdata"))
    return(Results)
  }else{
    return(NA)
  }
}


