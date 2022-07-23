# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

#####################################################################
########################### library #################################

library(fBasics) # tr()
library(MASS) 
library(invgamma) # digamma()


library(tensorflow)
library(tfprobability)
tryCatch({tfd = tfp$distributions}, error = function(e) e)
tryCatch({tfd = tfp$distributions}, error = function(e) e)
tryCatch({tfd = tfp$distributions}, error = function(e) e)
# tfd = tfp$distributions
tryCatch({test1 = tf$constant(1)}, error = function(e) e) 

#####################################################################
########################### Code ####################################

# The following arguments are required to customize your model:
#   
# - the survival response: `y`
# - the status: `delta`
# - the number of events: `d` 
# - the data in the dimension (P,N): `X` with P the number of genes and N the number of observations
# - the allowed number of latents factors: `Kmax`
# - the number of iterations of the variational inference algorithm: `maxit`
# - the seed of the algorithm: `seed`
# - the variance hyperparameter of W matrix: `inivarW`
# - the variance hyperparameter of $\Phi$ matrix: `iniVarPhi`
# - the variance hyperparameter of $\E$ matrix: `tau_error`
# - the variance hyperparameter of the gamma distribution: for $\Phi$ (`cc` and `dd`) ; for W (`e0` and `f0`)
# - the variance hyperparameter of the beta distribution: `a` and `b`
# - the number of sampling in the tensorflow algorithm: `S` 
# - the learning rate in the tensorflow algorithm: `learnings_rate` 

Weib_isgICA_model <- function(y,delta,d, X,Kmax,maxit, seed=1, 
                              iniVarW, iniVarPhi,a=1,b=Kmax,
                              cc=10^-1,dd=10^-1,e0=10^-1,f0=10^-1,tau_error=1,
                              S = 20, learning_rate= 0.01){
  ### Note : result$stock : too large, keep only last update
  EWTW <- c()
  J <- c()
  W_sigma <- list()
  Phi_sigma <- list()
  set.seed(seed)
  eps=2e-8
  P = dim(X)[1]
  N = dim(X)[2]
  initKmax = Kmax=min(Kmax,N,P)
  
  ######## hyperparameter for A matrix (covariance)  ########
  c0 =cc*matrix(1,nrow = P, ncol = Kmax) 
  d0 = dd*matrix(1,nrow = P, ncol = Kmax) 
  c = c0; d = d0;
  lgc0 = lgamma(c0) 
  
  e=e0;f=f0;
  
  # standardize the data
  X = t(scale(t(X)))
  tX = t(X)
  
  ######## Initialize A and S   ######## 
  #  X = U D V'
  
  SVD <- svd(X)
  Phi = SVD$u %*% diag(SVD$d)
  Phi = Phi[,1:Kmax]
  W = t(SVD$v)
  W = W[1:Kmax,]
  
  ##### Init random
  
  #Phi = matrix(rnorm(Kmax*P, sd=2), nrow = P, ncol = Kmax)
  #W = matrix(rnorm(N*Kmax, sd=2), nrow = Kmax, ncol = N)
  
  Z=matrix(T,nrow = Kmax, ncol = N) 
  Z <<- Z
  Zb <<- Z > 0.5
  ZW=Z*W;  
  sigmav = matrix(1, nrow = P, ncol = Kmax)
  
  tZW = t(ZW)
  
  V_error = matrix(tau_error, nrow = P, ncol = 1)
  vec_tr=c()
  
  ######## initialize pi  ######## 
  
  KmaxZ = Kmax-1 # mise en place du profil moyen
  
  pia0=c(1,a/initKmax* matrix(1, nrow = 1,ncol = KmaxZ) )
  pib0=c(1, b*(initKmax-1)/initKmax*matrix(1, nrow = 1,ncol = KmaxZ) )
  pia=pia0;pib=pib0;
  pia0_profil = pia0[-1]
  pib0_profil = pib0[-1]
  lgpia0=lgamma(pia0_profil)
  lgpib0=lgamma(pib0_profil)
  lgpiab0=lgamma(pia0_profil+pib0_profil)
  pai= 0.1*matrix(1, nrow = 1,ncol = Kmax) 
  
  Phi_sigma=matrix(1e-6, nrow = P,ncol = Kmax) #convariance of Phi_jk
  EPhiPhi = Phi*Phi+Phi_sigma;
  PhiV_error = Phi*matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax)
  PhiTPhi=t(PhiV_error)%*%Phi
  EPhiTPhi=t(Phi)%*%Phi+diag(colSums(Phi_sigma))
  sigmaw=1
  DWvar=matrix(0,Kmax,N)
  
  
  ############################# TENSORFLOW INITIALIZATION #############################
  # K = dim(Z)[2]
  # N = dim(Z)[1]
  y2 = log(y)
  dimK = Kmax
  Xtensor = tf$convert_to_tensor(t(Z),dtype="float64")
  yTensor = tf$reshape(tf$math$log(tf$convert_to_tensor(y, dtype="float64")),  shape(N,1))
  evtTensor = tf$convert_to_tensor(delta,dtype="float64")
  set.seed(seed)
  tf$random$set_seed(seed)
  #optimizer choice
  opt=tf$keras$optimizers$Adam(learning_rate=learning_rate,epsilon=1e-8)
  ############################# Init pars
  rhoLocs=tf$Variable(tf$zeros(1L,dtype=tf$dtypes$float64))
  rhoScs=tf$Variable(tf$constant(-5L,dtype=tf$dtypes$float64))
  nuLocs = tf$Variable(tf$constant(dimK,dtype=tf$dtypes$float64)) # based on conjugate prior / posterior value
  kappaLocs = tf$Variable(tf$constant(dimK,dtype=tf$dtypes$float64)) # based on conjugate prior / posterior value 
  TauScs = tf$Variable(tf$constant(exp(diag(dimK)#*2.2
  ),shape=shape(dimK,dimK),dtype=tf$dtypes$float64))
  
  betaLocs=tf$Variable(tf$zeros(dimK,dtype=tf$dtypes$float64))# betaScs = tf$linalg$pinv(kappaLocs*nuLocs*TauScs) # kappa * nu * V
  
  VLocs = tf$Variable(tf$matmul(TauScs, tf$transpose(TauScs))) # V = L %*% t(L) cholesky decomposition
  # diag(as.matrix(VLocs))
  betaScs = (tf$linalg$pinv(kappaLocs*nuLocs*VLocs)) # kappa * nu * V
  
  ############################# Priors
  tfp.multinorm = tfp$distributions$MultivariateNormalFullCovariance
  V0 = diag(dimK) * 0.1#
  invV0 =  tf$linalg$pinv(V0)
  nu0 = dimK
  log_sum_gamma0 = log(sum(gamma( (nu0+1-(1:dimK))/2)))
  chol_scale2 = tf$linalg$cholesky(tf$constant(V0,dtype=tf$dtypes$float64))
  TauPrior = tfp$distributions$WishartTriL(df=tf$constant(nu0,dtype=tf$dtypes$float64), scale_tril=chol_scale2)
  kappa0 = 1 #dimK
  S0 = tf$linalg$pinv(diag(dimK)) * 0.1 #kappa0*nu0*V0) #
  betaPrior = tfp.multinorm(loc = tf$reshape(tf$zeros(dimK,dtype=tf$dtypes$float64),shape(dimK,1)),
                            covariance_matrix = S0)
  rhoPrior=tfd$LogNormal(tf$zeros(1L,dtype=tf$dtypes$float64),tf$Variable(tf$ones(1L,dtype=tf$dtypes$float64)))
  
  ############################# VI dist
  betas = tfp.multinorm(betaLocs, betaScs) # betas=tfd$Normal(betaLocs,exp(betaScs))
  rhos=tfd$LogNormal(rhoLocs,exp(rhoScs))
  chol_scaleTau = tf$linalg$cholesky(VLocs)
  taus = tfp$distributions$WishartTriL(df=nuLocs, scale_tril=chol_scaleTau)
  
  
  r = rep(0, dim(betaLocs)[1] )#stock beta mean
  rsc=rep(0,dimK) #np.vstack((asc,betaScs))
  rsc_all = list(matrix(0, dimK, dimK)) #stock beta scale
  w = 0#stock rho mean
  wsc = 0#stock rho scale
  
  ############################# Loop ############################# 
  
  result <- list(
    W=list(),
    Z=list(),
    ZW=list(),
    Phi=list(),
    pai=list(),
    V_error=list(),
    mse=numeric(maxit), #c(),
    K=numeric(maxit), # c(),
    elbo=numeric(maxit), # c(),
    sigmav = list(),
    sigmaw = list(),
    hyperparameter = list(),
    stockSurv = list(),
    lastElbo = c()
  )
  iter = 0; change = 1; 
  start_time <- Sys.time()
  
  ######## While ########
  while (iter<maxit){
    iter = iter + 1;
    
    ######## update Z  ########
    ############ ICA information
    Zcnt = rowSums(Z) 
    exp_pi = cbind(digamma(pia0 + Zcnt) - digamma(pia0 + pib0 + N),
                   digamma(pib0 + N - Zcnt) - digamma(pia0 + pib0 + N) )
    
    PhiV_error = Phi*matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax)
    PhiTPhi=t(PhiV_error)%*%Phi
    EPhiPhiV_error = colSums(EPhiPhi * matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax))
    nonZeroPos=which(Zcnt> eps)
    Z[which(Zcnt<=eps),]=0
    nonZeroPos=nonZeroPos[nonZeroPos!=1]
    # print(c(dim(Z),length(nonZeroPos), dim(PhiV_error)))
    
    if(length(nonZeroPos)!=0){
      if(length(nonZeroPos)>1){
        obj=tX%*%PhiV_error[,nonZeroPos] - tZW%*%PhiTPhi[,nonZeroPos]  + t(diag(PhiTPhi)[nonZeroPos]*ZW[nonZeroPos,])
      }else{
        #insure the matrix object format
        obj=tX%*%matrix(PhiV_error[,nonZeroPos], ncol=1) - tZW%*%matrix(PhiTPhi[,nonZeroPos], ncol=1)  + matrix(t(diag(PhiTPhi)[nonZeroPos]*ZW[nonZeroPos,]), ncol=1)
      }
    }
    Z <<- Z
    ############ survival information
    alpha = c(betaLocs$numpy())
    rho = c(exp(rhoLocs$numpy()))
    alphaZ = as.matrix(alpha,nrow=1)
    lp =  t(Z) %*% alphaZ
    lpk = sapply(1:nrow(Z), function(k) y2 - t(Zb[-k,]) %*% alphaZ[-k] )
    lpk = exp(t(lpk) / rho)
    
    lpi =  sapply(1:ncol(Z), function(i) t(Zb[,i]) * t(alphaZ)*rho ) # b * Zi / rho : N vector of K value
    dm = matrix(rep(delta,Kmax),ncol=N, byrow=T)
    obj2 = - dm * lpi - exp(lpk)*exp(-lpi)
    
    # (nonZeroPos,N)
    t1 = exp_pi[nonZeroPos,2]-(exp_pi[nonZeroPos,1] - 0.5*(( (W[nonZeroPos,]^2) + DWvar[nonZeroPos,])*EPhiPhiV_error[nonZeroPos]) + t(obj)*W[nonZeroPos,])
    t1 = t1 - obj2[nonZeroPos, ] #- matrix(rep(log_L,nrow(t1)),ncol=N, byrow=T)
    Z[nonZeroPos,]=1/(1+exp(t1))  
    
    Zb <<- Z > 0.5
    
    ######## Update for W    ######## 
    logdetWsig=0;
    t2=PhiTPhi
    diag(t2)=EPhiPhiV_error+(sigmaw+1);
    
    midval=tryCatch(solve(t2), 
                    error = function(e){
                      return(ginv(t2))
                    }) 
    
    t1 = t(PhiV_error)%*%X 
    diagMidval=diag(midval)
    DWvar=diagMidval*matrix(1,nrow=Kmax,ncol=N)
    W = midval%*%t1
    WW = W*W
    EWTW=colSums(WW) + sum(diagMidval) 
    logdetWsig=0.5*N*log(2*pi*(det(midval) +eps)) 
    
    
    ZW=Z*W
    tZW = t(ZW)
    ZWZW = ZW%*%tZW
    
    ######## Update  Phi_jk    ######## 
    
    magWZ = rowSums(Z*(WW + DWvar))
    Phi_sigma=1/(V_error%*%magWZ+sigmav)
    for(k in 1:Kmax){
      Phi[,k]=0
      Xmid = X%*%ZW[k,] - Phi%*%ZWZW[,k]
      Phi[,k]=Phi_sigma[,k]*V_error*Xmid 
    }
    EPhiTPhi = t(Phi)%*%Phi + diag(colSums(Phi_sigma))
    EPhiPhi = Phi*Phi + Phi_sigma
    
    
    ######## Update survival coefficients    ######## 
    
    
    Xtensor = tf$convert_to_tensor(t(Z),dtype="float64")
    with(tf$GradientTape() %as% tape, {
      VLocs = tf$matmul(TauScs, tf$transpose(TauScs)) # V = L %*% t(L) cholesky decomposition
      chol_scaleTau = tf$linalg$cholesky(VLocs)
      taus = tfp$distributions$WishartTriL(df=nuLocs, scale_tril=chol_scaleTau)
      betaScs = tf$linalg$pinv(kappaLocs*nuLocs*VLocs) # kappa * nu * V
      betas=tfp.multinorm(betaLocs,betaScs)     # betas=tfd$Normal(betaLocs,tf$exp(betaScs))  
      rhos=tfd$LogNormal(rhoLocs,tf$exp(rhoScs))
      beta_samples=betas$sample(S)
      rho_samples=rhos$sample(S)
      lp=tf$linalg$matvec(beta_samples,Xtensor)
      rho=tf$transpose(rho_samples)
      diff= yTensor-lp
      loglik=tf$math$reduce_sum(tf$transpose(tf$math$multiply(evtTensor,tf$transpose(tf$math$log(rho)+rho*diff)))-tf$math$exp(diff)**rho,1L)
      betaKL=tf$reduce_sum(tfd$kl_divergence(betas,betaPrior))
      rhoKL=tf$reduce_sum(tfd$kl_divergence(rhos,rhoPrior))
      
      tmp_wishart = tf$matmul(invV0,VLocs) 
      tauKL = - nu0/2*log(abs(tf$linalg$det(tmp_wishart))) + nuLocs/2*(tf$linalg$trace(tmp_wishart) - dimK) + log_sum_gamma0 - log(sum(gamma( (nuLocs$numpy()+1-(1:dimK))/2))) + (nuLocs - nu0)/2*sum(digamma(nuLocs$numpy()/2 + (1-1:dimK)/2 ))
      loss=(rhoKL+betaKL+tauKL)/N-tf$math$reduce_mean(loglik)
    })
    pars = list(betaLocs,#
                nuLocs, # need nu > dimK
                kappaLocs,
                TauScs,rhoLocs,rhoScs) 
    r= rbind(r, t(as.matrix(betaLocs))) 
    rsc=rbind(rsc, diag(t(betaScs$numpy() ))) #as.matrix(betaScs)))) 
    rsc_all=append(rsc_all, list(betaScs$numpy())) #as.matrix(betaScs))) #
    w=rbind(w, t(rhoLocs$numpy())) #as.matrix(rhoLocs))) 
    wsc=rbind(wsc, t(rhoScs$numpy() )) #as.matrix(rhoScs))) 
    grads=tape$gradient(loss,pars)
    opt$apply_gradients(purrr::transpose(list(grads, pars  )))
    
    ######## Update  pi_k    ######## 
    
    pia=rowSums(Z)+pia0
    pib=N-rowSums(Z)+pib0
    pai=pia/(pia+pib)
    
    ######## Update update sigma Phi    ######## 
    
    c=c0+0.5;
    d=d0+0.5*EPhiPhi;
    sigmav=c/d;  
    
    ######## Update sigma W    ######## 
    
    e=e0 + Kmax*N/2; #sum(apply(Z>=0.5,1, sum)>=1)*N/2; #
    f=sum(EWTW)/2 + f0;
    sigmaw = e/f;
    Zb = Z > 0.5
    
    ######## Save samples.   ########
    
    res=X-Phi%*%ZW
    result$mse[iter] = mean(sqrt(colSums(res*res)))
    result$K[iter] = sum(rowSums(Zb) >0) #sum(Zcnt > 0)
    
    # if(iter>49 & ((iter%%100) ==0)){ 
    if(((iter%%100) ==0)){
      result$W = W
      result$Phi = Phi
      result$Z = Z
      result$ZW = ZW
      result$pai = pai
      result$V_error = V_error
      result$sigmav = sigmav
      result$sigmaw = sigmaw
      result$hyperparameter = list(pia,pib,c,d,e,f)
      pars.numpy = list(betaLocs = betaLocs$numpy(),#
                        nuLocs = nuLocs$numpy(), # need nu > dimK
                        kappaLocs = kappaLocs$numpy(),
                        TauScs= TauScs$numpy(),
                        rhoLocs=rhoLocs$numpy(),rhoScs= rhoScs$numpy()) 
      result$stockSurv = list( data = list(Z=Z,y=y,delta=delta),
                               # stock = list(betaLocs = r[iter,],
                               #              betaScs_var = rsc[iter,],
                               #              betaScs = rsc_all[iter],
                               #              rhoLocs = w[iter,],
                               #              rhoScs = wsc[iter,]),
                               pars = pars.numpy)
      save(result,X,a,b,cc,dd,e0,f0,y,d,delta, file = paste0("./result/inprogress/", "N_",N,"_P_",P, "_Kmax_",Kmax,"_a_prior_",a, "_b_prior_",b, "_ongoing", ".rdata"))
      
      ############################ Resize the objects #######################################
      
      
      Zcnt = rowSums(Z)
      # print(Zcnt)
      valZ = 1e-3
      nline = which(Zcnt<valZ )
      line = which(Zcnt>=valZ )
      # initKmax
    }
    
    if( (iter%%20) == 0){
      print(paste("Elapsed time: ", round(Sys.time()-start_time,4),"; Iteration", iter))
      print(paste("K", result$K[iter],"; mse", round(result$mse[iter],4)))
      start_time <- Sys.time()
      # print(i)
      outTf = cbind(w,r) #means from tf
      outTfscs = cbind(wsc, rsc) #scales from tf
      burnin = 10
      vecI = burnin:iter
      par(mfrow=c(4,7),mar=c(2,2,0.5,0.5))
      tryCatch({
        for(j in c(1:28)){ #ncol(outTf))[outTf[iter,] != 0] ){
          # print(outTfscs[vecI,j])
          if(j ==1){
            vecCI = c(outTf[vecI,j]-1.96*sqrt(exp(outTfscs[vecI,j])),rev(outTf[vecI,j]+1.96*sqrt(exp(outTfscs[vecI,j]))))
          }else{
            vecCI = c(outTf[vecI,j]-1.96*sqrt(outTfscs[vecI,j]),rev(outTf[vecI,j]+1.96*sqrt(outTfscs[vecI,j])))
          }
          plot(vecI,outTf[vecI,j],type="l",ylim=c(min(unlist(c(vecCI))),
                                                  max(unlist(c(vecCI)) ) ))
          polygon(c(vecI,rev(vecI)), vecCI, col = rgb(0,0,0,0.2), border=NA)
          abline(h=0, col=rgb(0,0,0,0.2), lty=1)
        }
      })
      
    }
  }
  print(dim(rsc))
  
  
  result$W = W
  result$Phi = Phi
  result$Z = Z
  result$ZW = ZW
  result$pai = pai
  result$V_error = V_error
  result$sigmav = sigmav
  result$sigmaw = sigmaw
  result$hyperparameter = list(pia,pib,c,d,e,f)
  pars.numpy = list(betaLocs = betaLocs$numpy(),#
                    nuLocs = nuLocs$numpy(), # need nu > dimK
                    kappaLocs = kappaLocs$numpy(),
                    TauScs= TauScs$numpy(),
                    rhoLocs=rhoLocs$numpy(),rhoScs= rhoScs$numpy()) 
  result$stockSurv = list( data = list(Z=Z,y=y,delta=delta),
                           stock = list(betaLocs = r[maxit,],
                                        betaScs_var = rsc[maxit,],
                                        betaScs = rsc_all[maxit],
                                        rhoLocs = w[maxit,],
                                        rhoScs = wsc[maxit,]),
                           pars = pars.numpy) 
  
  ######## Lower Bound    ########
  J9 = 0
  tmp2 = digamma(c) - log(d)
  tmp3 = digamma(e) - log(f)      
  
  pia_profil = pia[-1]
  pib_profil = pib[-1]
  pai_profil = pai[-1]
  
  # log likelihood 
  predX = Phi%*%ZW
  midval = apply(X*(X-2*predX),1,sum)
  tmp5 = ZWZW-diag(diag(ZWZW))+diag(magWZ)
  vec_tr = c()
  for(j in 1:P){
    tmp4 = Phi[j,]%*%t(Phi[j,]) # attention au probleme de R sur les dimensions
    tmp4 = tmp4-diag(diag(tmp4))+diag(EPhiPhi[j,])
    trtmp = tr(tmp5%*%tmp4)
    vec_tr[j] = 0.5*(trtmp +midval[j])
  }
  j1bis = -P*N/2*log(2*pi)+sum(N*P/2*log(V_error)-V_error*vec_tr) 
  
  # KLD of gamma distribution
  j3 = - sum( c*log(d) - c0*log(d0) - lgamma(c) +lgc0 + (c-c0)*(tmp2) - c*(1-d0/d) )
  j11 = -(e*log(f)-e0*log(f0)-lgamma(e)+lgamma(e0)+(e-e0)*tmp3-e*(1-f0/f))
  
  
  # KL de Pi_k
  piab=pia_profil+pib_profil
  dgpia=digamma(pia_profil)
  dgpib=digamma(pib_profil)
  dgpiab=digamma(piab)
  j5 = sum(lgpiab0-lgamma(piab)+lgamma(pia_profil)-lgpia0+lgamma(pib_profil)-lgpib0-(pia_profil-pia0_profil)*(dgpia-dgpiab)-(pib_profil-pib0_profil)*(dgpib-dgpiab))
  j4 = (dgpia-dgpiab)%*%rowSums(Z[-1,])+(dgpib-dgpiab)%*%(N-rowSums(Z[-1,]));
  tmpZ = Z
  tmpZ[which(Z<eps)]=eps
  j9 = sum(tmpZ*log(tmpZ))
  
  j2 =-0.5*P*Kmax*log(2*pi) + 0.5*sum(tmp2) - 0.5*sum((c/d)*(d-d0));
  j6 = (N*Kmax/2)*tmp3-0.5*sigmaw*sum(EWTW)
  j8 = 0.5*sum(log(2*pi*abs(Phi_sigma)))
  j10 = logdetWsig
  
  # surv KL and surv loss 
  betaKL=tf$reduce_sum(tfd$kl_divergence(betas,betaPrior))
  rhoKL=tf$reduce_sum(tfd$kl_divergence(rhos,rhoPrior))
  tmp_wishart = tf$matmul(invV0,VLocs) 
  tauKL = - nu0/2*log(abs(tf$linalg$det(tmp_wishart))) + nuLocs/2*(tf$linalg$trace(tmp_wishart) - dimK) + log_sum_gamma0 - log(sum(gamma( (nuLocs$numpy()+1-(1:dimK))/2))) + (nuLocs - nu0)/2*sum(digamma(nuLocs$numpy()/2 + (1-1:dimK)/2 ))
  lossSurv=(rhoKL+betaKL+tauKL)/N-tf$math$reduce_mean(loglik)
  
  result$lastElbo = j2+j3+j4+j5+j6+j8+j9+j10+j11 + j1bis + lossSurv$numpy() 
  result$ICAElbo = j2+j3+j4+j5+j6+j8+j9+j10+j11 + j1bis
  result$SurvElbo = lossSurv$numpy()+j4+j5+j9
  
  return(result)
}
