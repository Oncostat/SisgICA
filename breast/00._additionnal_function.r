library(class)
library(RcppHungarian)

sorting <- function(TH,new){
  # TH = phith # theorical Phi
  # new = pp # estimated Phi
  colnames(TH) = paste0("K",1:ncol(TH))
  colnames(new) = paste0("K",1:ncol(new))
  cor <- psych::corr.test( TH, new,method="pearson", adjust="none")
  cost = abs(cor$r)
  soln <- HungarianSolver(1-cost)
  return(list(new[soln$pairs[,2],],
              soln$pairs[,2]))
  }


accuracy_Z <- function(TH, new){
  # TH = Zth # theorical Z
  # new = zz2 # estimated Z
  # format to ignore the baseline
  TH=c(TH[-1,])
  new=c(new[-1,])
  
  recap = table(new,TH)
  
  acc = (recap[1,1]+recap[2,2]) / sum(recap)
  return(acc)
}


myImagePlot<-function(x,lines){
  x<-x+min(x)
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(1,0,length=100),  # Red
                    seq(1,0,length=100),  # Green
                    seq(1,1,length=100))  # Blue
  
  v1<-c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02")
  v2<-c("#7FCDBB","#40B6C4","#2C7FB8" ,"#253494") 
  pal1<-colorRampPalette(c(rev(v2),v1))
  nb_color_need = 100
  ColorRamp <- pal1(nb_color_need)
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  ColorLevels <- ColorRamp
  
  # Data Map
  image(1:ncol(x),1:nrow(x),t(x)[,nrow(x):1], col=ColorLevels, xlab="",
        ylab="", axes=F, zlim=c(min,max))
  par(xpd=T)
  text(1:ncol(x),-0.2,labels=xLabels,srt = 45, adj = 1)
  text(0.1,1:nrow(x),labels=yLabels)
  par(xpd=F)
  
  if(lines){ 
    for (i in (1:ncol(x))){
      abline(v=0.5+i,lty=2,col="grey")
    }	}
  box(lwd=2) 
  
}