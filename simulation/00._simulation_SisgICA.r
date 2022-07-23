# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

#####################################################################
########################### library #################################

library(mvtnorm)

#####################################################################
########################### generate data############################

simu_SisgICA <- function(P=50, N=100, Kmax=50,
                         w_var_prior=1,
                         phi_var_prior=1, error_var_prior=0.1,
                         rho=NULL, seed=1,Zth=NULL, K=10,
                         mean_phi=NULL,mean_W=NULL, Z_Phi=NULL,
                         val.shape.th=NULL, val.scale.thC=NULL,
                         coef_reg = NULL, censR = NULL){
  set.seed(seed)
  nbSimFactors = K
  if(is.null(Zth)){ 
    Zth = 1*(matrix(rnorm(K*N), nrow=K)>0)
  }
  if( K!= dim(Zth)[1]){
    stop("Dimension error between K and the row of Zth")
  }
  
  if(is.null(mean_phi)){ mean_phi <-  rep(0, nbSimFactors) }
  if(is.null(mean_W)){ mean_W <-  rep(0, nbSimFactors) }
  w_var <- rep(w_var_prior,nbSimFactors) # weight variance
  eta_var <- error_var_prior #error variance
  phi_var <- rep(phi_var_prior,nbSimFactors)
  
  Wth <- W <- t(rmvnorm(n=N, mean=mean_W, sigma = diag(w_var)))
  Eth <- E <- t(matrix(rnorm(N*P,0,eta_var),nrow=N,ncol=P))
  Phith <- phi<-matrix(rnorm(K*P,mean_phi,phi_var),nrow=P,ncol=K,byrow = T)
  
  if(is.null(Z_Phi)){
    B<-matrix(0,nrow = P,ncol = K)
    offset=0.05*P
    length=round(P/(K-1)+offset,1)
    end<-1
    for(i in (1:K)){
      start<-end-(i!=1)*offset
      end<-start+length-1
      B[start:min(end,nrow(B)),i]<-phi[start:min(end,nrow(B)),i]
    }
    end<-1
    for(i in (K:1)){
      start<-end-(i!=K)*offset
      end<-start+length-1
      B[start:min(end,nrow(B)),i]<-phi[start:min(end,nrow(B)),i]
    }
  }else{
    B = Z_Phi
    B = B*phi
  }
  Phith <- phi<- B
  
  X <- Xth <- Phith%*%(Wth*Zth)+Eth
  
  
  if(is.null(val.shape.th)){ val.shape.th =2 }
  
  if(is.null(val.scale.thC) & (! (is.null(censR))) ){ 
    
    val.scale.thC =  2
    stockScaleC = val.scale.thC
    d = 1
    
    while( abs(censR  - d) > 1e-2 ){
      # print(stockScaleC)
      
      if(d < censR ){
        stockScaleC = stockScaleC*2
      }else{
        stockScaleC = stockScaleC/3
      }
      
      set.seed(seed)
      
      if(is.null(coef_reg)){ coef_reg = rnorm(dim(Zth)[1],0,0.5) }else{
        if(length(coef_reg) != dim(Zth)[1]){
          stop("Dimension error between coefficients of survival model and the components of Zth") }
      }
      
      truetime<-rweibull(N,shape = val.shape.th, scale = exp(t(coef_reg) %*%Zth))
      cens<-rweibull(N,shape =val.shape.th, scale =stockScaleC)# distribution to adapt the censured rate
      time = pmin(truetime,cens) 
      delta<-1*I(truetime<cens) 
      status = delta
      d = sum(delta) / length(delta)
    }
    val.scale.thC = stockScaleC
    print(d)
    # print(stockScaleC)
    
  }else{
    
    set.seed(seed)
    
    if(is.null(coef_reg)){ coef_reg = rnorm(dim(Zth)[1],0,0.5) }else{
      if(length(coef_reg) != dim(Zth)[1]){
        stop("Dimension error between coefficients of survival model and the components of Zth") }
    }
    
    
    if(is.null(val.scale.thC) & is.null(censR)){ val.scale.thC =2 }
    
    truetime<-rweibull(N,shape = val.shape.th, scale = exp(t(coef_reg) %*%Zth))
    cens<-rweibull(N,shape =val.shape.th, scale =val.scale.thC) 
    time = pmin(truetime,cens)  
    delta<-1*I(truetime<cens) 
    status = delta
    d = sum(delta)
  }
  
  return(list("X"=X, 
              "Z"=Zth,
              "W"=Wth,
              "ZW"=Zth*W,
              "phi"=Phith,
              "w_var"= w_var_prior,
              "eta_var"= error_var_prior,
              "phi_var"= phi_var_prior,
              "K_sim"=K,
              "val.shape.th"=val.shape.th,
              "val.scale.thC" = val.scale.thC,
              "coef_reg"= coef_reg,
              "delta"=delta,
              "time" = time,
              "censth"=censR))
}



simu_SisgICA_save <- function( N=600,P=1200,K = 10,vec_var_error = 1,
                               seed.coef = 0,vec_i=100,
                               cens = 0.3){
  
  ###################### Init data ###################### 
  
  listZ = list()
  vec_seed = rep(vec_i,2)
  vec_sd_error = sqrt(vec_var_error)
  value_sd_error <<- sqrt(vec_var_error)
  value_var_error <<- vec_var_error
  
  ###################### ALGORITHM #######################
  
  set.seed(1)
  BZ<-matrix(0,nrow = K,ncol = N)
  end<-1
  
  
  ############ for K=10
  set.seed(vec_seed[1]); (group = rbinom(n = K,size = 15,prob = 0.5)+1)
  pinit = sapply(group, function(i) sample(1:N,size = i,replace = T)) 
  offset = sapply(group, function(i) sample(1:40,size = i,replace = T)) #p150
  
  
  for(i in (1:K)){
    group[1]
    pinit[[1]]
    offset[[1]]
    for(j in 1:group[i]){
      BZ[i,pinit[[i]][j]:min(ncol(BZ),pinit[[i]][j]+offset[[i]][j] )] = 1
    }
  }
  BZ[1,] = 1
  Zth2 = BZ
  table(c(BZ))[2] / sum(table(c(BZ)))
  
  Zth = Zth2
  type = 2 
  typeZ = "random"
  typeZ <<- c("diag","random")[type]
  print(dim(Zth))
  BZ = Zth 
  table(c(BZ))
  
  
  meltedZth <- melt(Zth > 0.5)
  ggplot(meltedZth, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +
    scale_fill_manual(values = c("white", "black")) +
    theme_bw() +
    theme(legend.position = "none")+ ggtitle("Real Z")+
    ylab("Facteurs latents") + xlab("Individus")
  
  listZ[[1]] = Zth
  K = nrow(listZ[[1]])
  val = 3
  mean_phi = seq(-K*val/K,K*val/K, length=K)
  mean_W = rep(0,K) 
  mean_W = mean_W[order(abs(mean_W))]
  
  B<-matrix(0,nrow = P,ncol = K)
  end<-1
  set.seed(vec_seed[2]); (group = rbinom(n = K,size = 10,prob = 0.5)+1)
  pinit = sapply(group, function(i) sample(1:P,size = i,replace = T)) 
  offset = sapply(group, function(i) sample(1:150,size = i,replace = T))  
  
  for(i in (1:K)){
    group[1]
    pinit[[1]]
    offset[[1]]
    for(j in 1:group[i]){
      B[pinit[[i]][j]:min(nrow(B),pinit[[i]][j]+offset[[i]][j] ) ,i] = 1
    }
  }
  
  
  set.seed(seed.coef)
  coef_reg = c(sample(c(-1,-0.5,-0.1,0.1,0.5,1, rep(0, length=K-6-1)),replace = F),0)
  
  sim = simu_SisgICA(N = N,P = P,
                     Kmax = Kmax,seed = vec_seed[2],
                     Zth=Zth,K=K, mean_phi = mean_phi,
                     mean_W = mean_W,Z_Phi = B,
                     error_var_prior =value_sd_error,
                     coef_reg = coef_reg, cens = cens)
  
  X=sim$X
  Zth=(sim$Z>0.5)*1
  Wth=sim$W
  ZWth=sim$ZW
  phith=sim$phi
  w_var_Sim= sim$w_var
  eta_var_Sim= sim$eta_var
  phi_var_Sim= sim$phi_var
  K_sim=sim$K_sim
  Kmax = min(100,P,N)
  maxit = 1000
  
  
  val.shape.th=sim$val.shape.th
  val.shape.thC = sim$val.shape.thC
  coef_reg= sim$coef_reg
  delta <<- sim$delta
  (d <<- sum(delta))
  y <<- sim$time
  
  summary(fit2<-survreg(Surv(y,delta) ~ t(Zth),dist="weibull")) # intercept = Z[,1]
  betas_weib<-fit2$coefficients #extraction de la valeur des betas
  lambda_weib2<-exp(betas_weib) # HRs
  scale_weib2<-fit2$scale
  sigma_weib2 = 1/scale_weib2 #
  
  ic_freq2 = data.frame(IC_inf = c(betas_weib,scale_weib2) #coef(fit2) 
                        -1.96*sqrt(diag(vcov(fit2))), #[-11],
                        mean =c(betas_weib,scale_weib2),# coef(fit2),
                        IC_sup = c(betas_weib,scale_weib2) #coef(fit2) 
                        +1.96*sqrt(diag(vcov(fit2)))) #[-11])
  ic_freq2$sign = ifelse((ic_freq2$IC_inf > 0 & ic_freq2$IC_sup > 0) | (ic_freq2$IC_inf < 0 & ic_freq2$IC_sup < 0), "*", ".")
  
  ic_freq2 = ic_freq2[-2,]
  ic_freq2$th_inside = ifelse((ic_freq2$IC_inf <  c(coef_reg,log(val.shape.th)) & ic_freq2$IC_sup >  c(coef_reg,log(val.shape.th))), "1", ".")
  ic_freq2$th = c(coef_reg,log(val.shape.th))
  ic_freq2
  
  table(c(B))[2] / sum(table(c(B)))
  
  print(c("Z" = table(c(BZ))[2] / sum(table(c(BZ))), 
          "Phi_Z" = table(c(B))[2] / sum(table(c(B)))))
  
  save(sim,B,BZ,vec_sd_error,vec_i,type, typeZ,value_var_error,value_sd_error, 
       fit2,ic_freq2,y,d,delta,coef_reg,val.shape.th,val.shape.thC,
       file = paste0( "simu_N_",N,"_P_",P,"_K_",K,
                      "_censured_",cens,
                      "_var_error_",value_var_error,
                      "_seed_",vec_seed[1],
                      "_typeZ_",typeZ,"_Z_",table(c(BZ))[2] / sum(table(c(BZ))), "_Zphi_",table(c(B))[2] / sum(table(c(B))),".rdata"))
  
}