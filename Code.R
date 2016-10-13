#Código Xtraining
#Author: Said Muñoz Montero
#13/Oct/2016

#Installing packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("oligo","pd.hg.u95av2"))

#Funciones adicionales
volcanoplot2<-function(TT,M,B){
  lfc.status = TT$logFC
  B.status = TT$B
  fitCB=data.frame(coef=lfc.status,lods=B.status)
  x0 = min(lfc.status) -.5
  x1 = max(lfc.status) +.5
  y0 = min(B.status) -.5
  y1 = max(B.status) +.5
  plot(lfc.status,B.status,col="black", ylim=c(y0,y1), 
       xlim=c(x0,x1), main=contraste, pch=16,cex=.35,
       cex.lab=1.2,xlab="Log Fold Change",
       ylab="B statistic")
  par(new=T)
  abline(v=-M, col="brown", ylab="", xlab="")
  par(new=T)
  abline(v=M, col="brown", ylab="", xlab="")
  par(new=T)
  abline(h=B, col="black", ylab="", xlab="")
  par(new=T)
  ind1 = abs(fitCB$coef)>M
  ind2 = fitCB$lods >B
  ind3 = (fitCB$coef>M & fitCB$lods>B)
  ind4 = (fitCB$coef< -M & fitCB$lods>B)
  x = as.matrix(fitCB$coef[ind1])
  y = as.matrix(fitCB$lods[ind1])
  plot(x, y, col="magenta",ylim=c(y0,y1), 
       xlim=c(x0,x1),main="", pch = 20, xlab="", 
       ylab="",cex.lab=1.2)
  x = as.matrix(fitCB$coef[ind2])
  y = as.matrix(fitCB$lods[ind2])
  par(new=T)
  plot(x, y, col="orange",  ylim=c(y0,y1), 
       xlim=c(x0,x1), main="", pch = 20, xlab="",
       ylab="",cex.lab=1.2)
  x = as.matrix(fitCB$coef[ind3])
  y = as.matrix(fitCB$lods[ind3])
  par(new=T)
  plot(x, y, col="red",  ylim=c(y0,y1), 
       xlim=c(x0,x1), main="", pch = 20, xlab="", 
       ylab="",cex.lab=1.2)
  x = as.matrix(fitCB$coef[ind4])
  y = as.matrix(fitCB$lods[ind4])
  par(new=T)
  plot(x, y, col="darkgreen", ylim=c(y0,y1), 
       xlim=c(x0,x1), main="", pch = 20, xlab="", 
       ylab="",cex.lab=1.2)
}

#Loading package
library("oligo")
celfiles<-list.celfiles("Data",full.name=TRUE)
data<-read.celfiles(celfiles)
colnames(data)<-gsub(".CEL","",colnames(data))
colores<-rep(c("red","blue"),each=2)

#Plots raw data
boxplot(data,col=colores,main="Boxplot datos crudos",las=2)
hist(data,col=colores,main="Histograma datos crudos")

#RMA
data.rma<-rma(data,background=TRUE, normalize=TRUE)
data.rma

#Plots RMA
boxplot(data.rma,col=colores,main="Boxplot RMA",las=2)
hist(data.rma,col=colores,main="Histograma RMA")


#RMA no background
data.rma.nobg<-rma(data,background=FALSE, normalize=TRUE)

#Plots RMA
boxplot(data.rma.nobg,col=colores,main="Boxplot RMA",las=2)
hist(data.rma.nobg,col=colores,main="Histograma RMA")

#Análisis diferencial no mostrado

#Revisar DEG
data.limma<-read.csv("DEG.csv")
volcanoplot2(data.limma,1,2)

#Extraer Genes DE
dge<-subset(data.limma,abs(logFC)>1 & B>2)
dge
#Genes diferencialmente expresados
nrow(dge)
