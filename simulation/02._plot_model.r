
# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

cat("\014") # clear console
rm(list=ls())
graphics.off()

#####################################################################
################# Path to the code directory ########################

# setwd("...")

#####################################################################
########################### library #################################


source("00._additionnal_function.r")
library(plotly)
library(seqinr)
library(ggplot2)
library(gridExtra)
library(reshape2) # melt
library("dcov") # dcor()
library(stringr)
library(survival)
library(ggpubr) # ggarrange

#####################################################################
############################# Results ###############################

##### evall object : a model summary
##### list_recap2 : summary of the survival submodel, comparison between IC with variational inference and frequentist inference


load( list.files(path="result/", recursive = F, full.names = T, pattern = "result_all")) # with Results and InitPar
Results$scoreSummary$beta = 1/ qbeta(Results$scoreSummary$U,0.25,2.5) - 1
selectBeta = Results$scoreSummary$beta[  which.max(Results$scoreSummary$Score)]

ll = list.files(path="result/", recursive = F, full.names = T, pattern = ".rdata")
file_simu = ll[str_detect(ll, pattern = paste(c("b_prior",selectBeta, "ite"),collapse = "_" )) ]
load(file_simu) 


evall = data.frame(matrix(NA, ncol=10, nrow = 1)) 
colnames(evall) = c(
  "a_prior","b_prior",
  "c","d","e","f", 
  "elbo", # !!!!!!!!!!! survELBO
  "K", 
  "ite", "mse")

evall[1,1:ncol(evall)] = c(
  initPar$a_prior,    # initPar$b_prior,
  selectBeta,
  initPar$cc,
  initPar$dd,
  initPar$e,
  initPar$f,
  fit$lastElbo,
  fit$K[length(fit$K)],
  length(fit$K), 
  fit$mse[length(fit$K)]
)
evall


zz = fit$Z
ligne = round(apply(zz>0.5,1,sum),1) 
zz1 = (zz[ligne>=1,] > 0.5)*1
w = fit$W
w = w[ligne>=1,]
zw = fit$ZW
zw = zw[ligne>=1,] * zz1
pp = fit$Phi
pp = pp[,which(ligne>=1)] 


list_recap2 = list()
fitSto = fit$stockSurv
dimZ = dim(fitSto$data$Z)

# IC for beta parameters (variational inference)
m = fitSto$stock$betaLocs
sd = sqrt(fitSto$stock$betaScs_var)
ic = cbind(m-1.96*sd,
           m,
           m+1.96*sd)
ic = ic[ligne > 0,] 

# IC for rho parameter (variational inference)
m_se = fitSto$stock$rhoLocs
v_se = exp(fitSto$stock$rhoScs)
ic_se = cbind(m_se-1.96*sqrt(v_se),
              m_se,
              m_se+1.96*sqrt(v_se))

# IC rbind
CI_VI = rbind(ic_se,ic) 
selection_VI = (CI_VI[,1] < 0) != (CI_VI[,3] > 0)


# Extract  IC (frequentist inference)
summary(fit2<-survreg(Surv(initPar$y,initPar$delta) ~ t(zz1[,]),dist="weibull")) # intercept = Z[,1]
betas_weib<-fit2$coefficients #extraction of beta parameters
lambda_weib2<-exp(betas_weib) # HRs
scale_weib2<-fit2$scale
sigma_weib2 = 1/scale_weib2 # The "log(scale)" equals to 1/σ (1/sigma)
ic_freq2 = data.frame(IC_inf = c(betas_weib,scale_weib2) #coef(fit2) 
                      -1.96*sqrt(diag(vcov(fit2))), #[-11],
                      mean =c(betas_weib,scale_weib2),# coef(fit2),
                      IC_sup = c(betas_weib,scale_weib2) #coef(fit2) 
                      +1.96*sqrt(diag(vcov(fit2)))) #[-11])
ic_freq2$sign = (ic_freq2$IC_inf > 0 & ic_freq2$IC_sup > 0) | (ic_freq2$IC_inf < 0 & ic_freq2$IC_sup < 0)
ic_freq2 = ic_freq2[-2,] 
paraCI = ic_freq2[c(nrow(ic_freq2), 1:nrow(ic_freq2)-1),]
nbparms = nrow(paraCI)-1
paraCI$parms = c("log_se",paste0("coef",1:nbparms))

list_recap2[[1]] = cbind(CI_VI, #coverage_VI.th,
                         selection_VI,paraCI,
                         N = dimZ[2], #K = initPar$sim$K_sim, 
                         seed =initPar$vec_seed[1] )  







#####################################################################
########################### Table HR  #################################

forme_IC_95 = function(x){
  # cbind( IC inf, mean, IC supp)
  val = round(x,3)
  return(paste0(val[,2], " [",val[,1],", ", val[,3],"]"))
}

see1 = list_recap2[[1]]$selection_VI==T

see = cbind(see1, c("log_se","Baseline profile", paste0("Alteration component ", 1:(length(see1)-2))  ))

see = see[see[,1]==T,]

datt = forme_IC_95(list_recap2[[1]][see1,1:3])
datt2 = forme_IC_95(exp(list_recap2[[1]][see1,1:3]))
datt = cbind(see[,2],datt,datt2)

colnames(datt) = c("Components", "Parameters" , "Hazard ratio")

datt


#####################################################################
########################### Figures  #################################

v1<-c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02")
v2<-c("#7FCDBB","#40B6C4","#2C7FB8" ,"#253494")  
pal1<-colorRampPalette(c(rev(v2),"white",v1))
nb_color_need = 100
ColorRamp <- pal1(nb_color_need)




borne = round(as.numeric(abs(quantile(c(zw), 0.1))),2)
borne2 = round(as.numeric(abs(quantile(c(zw), 0.001))),2)
tzw = t(zw)
cat = (tzw > borne) + (tzw > - borne)
tzw_tri = tzw[order(cat[,1],cat[,2],cat[,3],cat[,4],
                    cat[,5],cat[,6],cat[,7],cat[,8],
                    cat[,9]),]
rownames(tzw_tri)=NULL
Mpp = melt( t(tzw_tri) )
rrr = summary(c(tzw_tri),-2)

p3 = ggplot(Mpp, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  theme_bw() +
  scale_fill_gradientn(colours = ColorRamp, trans = ggallin::pseudolog10_trans, 
                       limits=c(-max(abs(c(tzw))),max(abs(c(tzw)))), breaks=c(-borne2,-borne,borne,borne2))+ 
  ggtitle("Reconstructed ZW")+
  ylab("Components") + xlab("Individuals")+
  theme(legend.position="bottom")


p2 = ggplot( melt( t(tzw_tri!=0) ), 
            aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  scale_fill_manual(values = c("white", "black")) +
  theme_bw() +
  # theme(legend.position = "none")+ 
  ggtitle("Reconstructed Z")+
  ylab("Components") + xlab("Individuals")+
  theme(legend.position="bottom",
        legend.justification="center",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(10,10,10,10))+
  theme(legend.key = element_rect(color="black") )


borne = round(as.numeric(abs(quantile(c(pp), 0.1))),0)
borne2 = round(as.numeric(abs(quantile(c(pp), 0.0001))),0)
cat = (pp > borne) + (pp > - borne)
pppp2 = pp[order(cat[,1],cat[,2],cat[,3],cat[,4],
                 cat[,5],cat[,6],cat[,7],cat[,8],
                 cat[,9]),]
Mpp = melt((pppp2))
rrr = summary(c(pppp2),-2)

p4 = ggplot(Mpp, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  theme_bw() +
  scale_fill_gradientn(colours = ColorRamp, trans = ggallin::pseudolog10_trans, 
                       limits=c(-max(abs(c(pp))),max(abs(c(pp)))), breaks=c(-borne2,-borne,borne,borne2))+
  # theme(legend.position = "none")+ 
  ggtitle(bquote("Reconstructed " *Phi))+ #"Reconstructed Phi"
  xlab("Components") + ylab("Gene expressions") + 
  theme(legend.position="bottom")+ coord_flip()

 

# grid.arrange(p2,p3, p4, ncol=3)   

################################################################################

