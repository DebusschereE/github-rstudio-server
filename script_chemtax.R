rm(list=ls())
require(ade4)
require(maps)
require(mapdata)
require(dplyr)

#Data
t=read.csv("Cruise_data_2017_V9_complete_Flowcam_NGS.csv",h=T,row.names=1)

t$N.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_SiO2
t$N.P<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4
t$N.P.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4/t$Nutr_SiO2

colnames(t)
#select dataframe geographic coordinates, chemtax, abiotic parameters
#
df=t[,c(3:4, 5:10, 166:168,72,75:76 ,127,107:112)]
rownames(df)
#remove NA rows for biolplankton
df<-na.omit(df)
#biol
colnames(df)
biol=df[,16:ncol(df)]
#check for high number of zeros per column
colSums(biol == 0)
biol<-biol[,2:6]
#Geographic coordinates
xy=df[,c(2,1)]
colnames(df)
#Abiotic
abio=df[,3:15]

#Mapping abiotic data
#Large white squares, low values
#Small squares, intermediate values
#Large squares, high values
#Remark: a priori, is "relative humidity"
#a pertinent variable for aquatic organisms?
par(mfrow=c(3,3),mar=c(2,2,3,1))
for(i in 1:ncol(abio)){
  plot(xy,type="n",main=colnames(abio)[i],cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,scalewt(abio)[,i],cleg=0,add.p=T)
}
par(mfrow=c(3,3),mar=c(2,2,3,1))
for(i in 1:ncol(biol)){
  plot(xy,type="n",main=colnames(biol)[i],cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,scalewt(biol)[,i],cleg=0,add.p=T)
}
#PCA of biol
#Keep 3 axes
#"scale=F", no reduction, better for
#matrices with a unique measurement unit
#"scale=T", systematically in the opposite case
#Here, keep 3 axes
pcaz=dudi.pca(log(biol+1),scale=F, scannf=FALSE, nf=3)

#Eigenvalues expressed in %
100*pcaz$eig/sum(pcaz$eig)

#biol variables
s.arrow(pcaz$co,clab=1)

#biol stations
s.label(pcaz$li)
s.label(pcaz$li,add.p=T)

#Correlation of abiotic variables?
cor(abio,pcaz$li)
## N.P; Salinity have the highest correlation with axis 1
## STIKSTOF, fosfaat, current velocity, depth to be correlated with Axis 2
## stikstof, N.P.Si--> axis3
#visualisations:
#axes 1 and 2:
par(mfrow=c(3,3))
for(i in 1:ncol(abio)){
  s.value(pcaz$li,scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#The same with axes 1 and 3
par(mfrow=c(3,3))
for(i in 1:ncol(abio)){
  s.value(pcaz$li[,c(1,3)],scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#PCA on Instrumental Variables (PCAIV)
#In this numerical context, you can say that PCAIV <=> RDA
#Keep 3 axes
iv=pcaiv(pcaz,abio, scannf=FALSE, nf=3)
iv

#Proportion of explained variance (R squared)
sum(iv$eig)/sum(pcaz$eig)
##65% explained
#Global summary (axes 1 and 2)
plot(iv)
dev.off()
scatter(iv,clab.col =0.5)
#Global summary (axes 1 and 3)
plot(iv,yax=3)

#Predictions
iv$li

#Observations
iv$ls

#Discripancies
#ls = "lignes suppl?mentaires
#They results from the passive projection
#of the raw biol data onto the axes (arrow tips)
#Arrow length = lack of fitting = residual
s.match(iv$li,iv$ls)

#Abio variables
#Salinity is the most influential variable
dev.off()
s.corcircle(iv$cor)

#Permutation test
test=randtest(iv,999)
test
#significant permutation test: siginficant relationship between chemtax phytoplankton and abiotic parameters
plot(test)

#Mapping the three axes
#This enables to check a possible spatialization
#of the three gradients
par(mfrow=c(3,1),mar=c(2,2,3,1))
for(i in 1:3){
  plot(xy,type="n",main=paste ("Axis",i,sep=" "),cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,iv$li[,i],cleg=0,add.p=T)
  box()
}

###spatial analysis ####

require(adegraphics)
require(spdep)
require(adespatial)
require(vegan)
install.packages("packfor", repos="http://R-Forge.R-project.org") 
require(packfor)


#A relatively recent development generating
#spatial predictors hierarchically
#organized from large to small scale
#They are build simply from the xy geographic coordinates
#A specifically weighted distance matrix among stations is computed
#and spatial predictors are derived
#They are called Moran's Eigenvector Maps (MEM)
#They are perfectly independent (orthogonal)
library(spdep)
#Choosing a connection network among stations
load("chooseCN.Rdata")
cn=chooseCN(xy,res="listw",type=1,plot.nb=F)
cn
#cn<-connection.network
#Triangulation is the simplest one
dev.off()
plot(cn,xy)
map("worldHires",add=T)

#Note that some distance crossing lands could be removed
#You can remove them
nb=tri2nb(xy)
#The map becomes interactive
#Remove the links crossing the lands
#by clicking on the connected stations
nb=edit(nb,xy) ##R loopt hierop vast

#Plot again
cn=nb2listw(nb,style="W")
cnchemtax<-cn
save(cnchemtax, file="cnchemtax.RData")
load("cnchemtax.RData")
cn<-cnchemtax
plot(cn,xy)
map("worldHires",add=T)
library(adespatial)
#MEM computation
#There are 38 stations, so 38-1=37 spatial predictors
#They represent all the possible scales
#that can be investigated in this data set
umem=scores.listw(cn)
#Test the significance of Moran's I for MEMs
moran.randtest <- function(x, listw, nrepet = 999, ... ){
  
  if(missing(listw) & inherits(x, "orthobasisSp"))
    
    if(!is.null(attr(x, "listw")))
      
      listw <- attr(x, "listw")
    
    
    
    x <- as.matrix(x)
    
    if(!is.numeric(x))
      
      stop("x should contain only numeric values")
    
    
    
    if(NCOL(x) == 1){
      
      res <- moran.mc(x, listw = listw, nsim = nrepet, zero.policy = TRUE)
      
      res <- as.randtest(obs = res$statistic, sim = res$res[-(nrepet + 1)], call = match.call(), ...)
      
    } else {
      
      res <- apply(x, 2, moran.mc, listw = listw, nsim = nrepet, zero.policy = TRUE)
      
      res <- as.krandtest(obs = sapply(res, function(x) x$statistic), sim = sapply(res, function(x) x$res[-(nrepet + 1)]), call = match.call(), ...)
      
    }
    
    
    
    return(res)
    
    
    
}
memtest=moran.randtest(umem,cn,999,"two-sided")
#Keeping only MEMS with significant values
names(memtest)
Umem=umem[,memtest[[7]]<0.05]



#Test of biol vs xy correlation
randtest(pcaiv(pcaz,xy,scan=F),999)#net significant 0.04
#Test of abio vs xy correlation
randtest(pcaiv(dudi.pca(abio,scan=F),xy,scan=F),999)#significant
#They are strongly correlated to xy
#suggesting that processes expressed at larger
#scale than the scale of the studied area
#affect the data
#Hence, this influence must be removed
#by detrending data from xy (removing xy effect)
#through Orthogonal PCAIV
biol.o=pcaivortho(pcaz,xy,scan=F)
abio.o=pcaivortho(dudi.pca(abio,scan=F),xy,scan=F)

#The detrended PCAIV
iv=pcaiv(biol.o,abio.o$tab, scannf=FALSE, nf=3)
randtest(iv,999)##significant
#https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
#Finding spatial predictors of iv
#24 % of iv variance is explained by space (last line, AdjR2Cum)
###spatialized of bio
rda.mem1=rda(biol.o$tab,Umem)

R2a.mem1=RsquareAdj(rda.mem1)$adj.r.squared
fw.mem1=forward.sel(biol.o$tab,Umem,adjR2thresh=R2a.mem1,99999)
fw.mem1##8%
fw.mem1=fw.mem1[order(fw.mem1$order),]##selection of MEM 
##biotic and abio ##17%!! 
rda.mem=rda(iv$tab,Umem)
summary(rda.mem)
R2a.mem=RsquareAdj(rda.mem)$adj.r.squared
fw.mem=forward.sel(iv$tab,Umem,adjR2thresh=R2a.mem,99999)
fw.mem
fw.mem=fw.mem[order(fw.mem$order),]##selection of MEM 


#Selected predictors ##for each station
pre.mem=data.frame(Umem[,colnames(Umem)%in%fw.mem$variables==T])
colnames(pre.mem)=fw.mem$variables
pre.mem=data.frame(pre.mem[,order(fw.mem$R2,decreasing=T)])

#pre.mem=Umem[,colnames(Umem)%in%fw.mem$variables==T]
#pre.mem=pre.mem[,order(fw.mem$R2,decreasing=T)]
#Represent them
#These predictors evidence processes
#(1) dominantly along waves of circa 50 km (MEM5, best predictor)
#    (alternation of white and black squares)
#(2) slightly along waves of similar length (MEM4, MEM7 and MEM8) 
##wat bedoelt hij hiermee????
detach('package:adegraphics')
detach('package:spdep')
detach('package:adespatial')
detach('package:vegan')
detach('package:packfor')

par(mfrow=c(2,2),mar=c(3,3,3,3))
for(i in 1:ncol(pre.mem)){
  plot(xy,type="n",main=colnames(pre.mem)[i],
       cex.main=1.5,bty="n",xaxt="n",yaxt="n")
  map("worldHires",add=T)
  s.value(xy,pre.mem[,i],cleg=0,csize=1.5,add.p=T)
}

#PCAIV on MEMs
iv.mem=pcaiv(iv,pre.mem)
plot(iv.mem)

#Projection of abiotic variables
w=pcaiv(dudi.pca(abio.o$tab,scan=F),pre.mem,scan=F)
#The data (iv.abio$tab) are projected onto the iv.mem system
sup=t(w$tab)%*%as.matrix(iv.mem$l1)
s.arrow(sup)

#The spatialization of the ecological pattern is mostly
#explained along a first axis by MEM5 and MEM7
#They explain the increasing concentration of
#biolpl_Appendicularia and biolpl_Noctiluca densities
#with decreasing salinity and increasing N concentration
#in coastal areas, from the center of the area to the coasts
par(mfrow=c(2,3),mar=c(3,3,3,3))
for(i in 1:2){
  plot(xy,type="n",main=colnames(pre.mem)[i],
       cex.main=1.5,bty="n",xaxt="n",yaxt="n")
  map("worldHires",add=T)
  s.value(xy,pre.mem[,i],cleg=0,csize=1.5,add.p=T)
}
s.corcircle(iv.mem$cor,clab=1.3)
s.arrow(sup,clab=1,xlim=c(-50,40))
s.arrow(iv.mem$c1,clab=1,xlim=c(-1.2,1.3))
dev.off()
s.arrow(iv.mem$c1,clab=1,xlim=c(-1.2,1.5))


