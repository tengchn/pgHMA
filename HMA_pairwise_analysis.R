library(ape)
library(phangorn)
library(phytools)
library(ggplot2)
setwd("") ##change it to your own path
#1.Phylogenetic analysis based on HMA and nucleotide sequencing data
#Build a distance matrix for 16 samples 
df16<-data.frame(matrix(0,ncol=16,nrow=16))
x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
rownames(df16)<-x
colnames(df16)<-y
df16
#Import distance matrix of nucleotide data and HMA data
pair_data<-read.csv("Data/HMA-16-ober-real.csv")
head(pair_data)
#Construct neighbour-joining(NJ) tree based on DNA data
df16[lower.tri(df16)]<-pair_data$real
df16<-t(df16)
df16[lower.tri(df16)]<-pair_data$real
df16_DNA<-df16
tr_real_16<-nj(df16_DNA)
plot(midpoint.root(tr_real_16))
add.scale.bar()
write.tree(tr_ober_16,file = "Data/Tree_file/DNA_nj.tre")
#Construct NJ tree based on HMA data
df16[lower.tri(df16)]<-pair_data$ober
df16<-t(df16)
df16[lower.tri(df16)]<-pair_data$ober
df16_HMA<-df16
tr_ober_16<-nj(df16_HMA)
tr_ober_16$edge.length <- pmax(tr_ober_16$edge.length,0) #Change the negative edge length value to 0
plot(midpoint.root(tr_ober_16))
add.scale.bar()
write.tree(tr_ober_16,file = "Data/Tree_file/HMA_nj.tre")
#Compare these two trees Robinson-Foulds (RF) distance
RF.dist(tr_ober_16,tr_real_16)


#2.Correlation analysis between HMA and nucleotide sequencing data
#Plot the data
plot(pair_data$ober, pair_data$real, xlab="Heteroduplex distace", ylab="Nucleotide distance",cex = 1, pch = 21, bg = 'gray')
#Linear regression with prediction interval
lm.out <- lm(real ~ ober,data=pair_data)
conf_interval <- predict(lm.out, interval="predict")
ggplot(pair_data, aes(x=ober, y=real))+
  theme_bw(base_size = 15)+
  xlab(expression(paste("Heteroduplex mobility distance (",d[H], ")")))+
  ylab("Nucleotide distance")+
  geom_point(size = 3.5, shape = 21, bg = 'gray')+
  annotate(label = paste('y ==',round(lm.out$coefficients[2],3),'~ d[H] + ',round(lm.out$coefficients[1],3),"*plain(\",\")",'~~ R^2 ==', round(summary(lm.out)$r.squared,3)), geom = "text", x = 0.15, y = 0.18, size = 7, parse = TRUE)+
  geom_smooth(method="lm", size=1.5)+
  geom_line(aes(x=sort(ober), y=sort(conf_interval[,2])), col = "red",size=1)+
  geom_line(aes(x=sort(ober), y=sort(conf_interval[,3])), col = "red",size=1)
#Summary of the linear regression
summary(lm.out)
summary(conf_interval)


#3.Parametric bootstrap sampling using the predict interval
library(boot)
pair_data<-read.csv("Data/HMA-16-ober-real.csv")
sd<- vector()
np<-10000
n_pober<- vector()
n_pober_pos<- vector()
set.seed(12345)
for (j in 1:np){
  df16<-data.frame(matrix(0,ncol=16,nrow=16))
  x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
  y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
  rownames(df16)<-x
  colnames(df16)<-y
  for(i in 1:length(pair_data$real))
    {
  sd[i]<-(conf_interval[i,3]-conf_interval[i,1])/1.96
  n_pober[i]<-rnorm(1,conf_interval[i,1],sd[i])
  }
  n_pober_pos<-pmax(n_pober,0)
  df16[lower.tri(df16)]<-n_pober_pos
  df16<-t(df16)
  df16[lower.tri(df16)]<-n_pober_pos
  tr_ober16<-nj(df16)
  tr_ober16$edge.length <- pmax(tr_ober16$edge.length,0)
  write.tree(tr_ober16,file="Data/Tree_file/10000bp_tree_HMA_nj.tre",append = T)
}

#Phylogenetic analysis with the expected genetic distance matrix obtained from the linear regression
df16<-data.frame(matrix(0,ncol=16,nrow=16))
x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
rownames(df16)<-x
colnames(df16)<-y
df16[lower.tri(df16)]<-conf_interval[,1]
df16<-t(df16)
df16[lower.tri(df16)]<-conf_interval[,1]
df16_fitted<-df16
tr_ober16_fitted<-nj(df16_fitted)
plot(midpoint.root(tr_ober16_fitted))
add.scale.bar()
write.tree(tr_ober16_fitted,file = "Data/Tree_file/HMA_fitted_nj.tre")


#4.Skyline plot analyses
#Parametric bootstrap sampling using UPGMA for HMA data
library(boot)
pair_data<-read.csv("Data/HMA-16-ober-real.csv")
lm.out <- lm(real ~ ober, data=pair_data)
conf_interval <- predict(lm.out, interval="predict")
sd<- vector()
np<-10000
n_pober<- vector()
n_pober_pos<- vector()
set.seed(12345)
for (j in 1:np){
  df16<-data.frame(matrix(0,ncol=16,nrow=16))
  x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
  y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
  rownames(df16)<-x
  colnames(df16)<-y
  for(i in 1:length(pair_data$real))
  {
    sd[i]<-(conf_interval[i,3]-conf_interval[i,1])/1.96
    n_pober[i]<-rnorm(1,conf_interval[i,1],sd[i])
  }
  n_pober_pos<-pmax(n_pober,0)
  df16[lower.tri(df16)]<-n_pober_pos
  df16<-t(df16)
  df16[lower.tri(df16)]<-n_pober_pos
  tr_ober16<-upgma(df16)
  tr_ober16$edge.length <- pmax(tr_ober16$edge.length,0)
  write.tree(tr_ober16,file="Data/Tree_file/10000bp_tree_HMA_upgma.tre",append = T)
}
#Load the UPGMA trees file with first 10% trees discarded
all_upgma_tre<-read.newick(file="Data/Tree_file/10000bp_tree_HMA_upgma.tre")
all_upgma_tre<-all_upgma_tre[c(floor(0.1 * length(all_upgma_tre)):length(all_upgma_tre))]
#Set the function for finding the lower 95% HPD of the matrix
library(HDInterval)
get_lower_HPD <- function(x,y) {
all_tree_height<-list()
for(i in 1:length(x)){
  all_tree_height[[i]]<-tail(skyline(x[[i]],y)$time,n=1)
}
lower_HPD<-hdi(do.call(rbind, all_tree_height))[1]
return(lower_HPD)
}
#Set the function for generating the classic skyline plot matrix
get_fake_beast_para <- function(i,x,y,lower_HPD) {
  tmp <- data.frame(skyline(x[[i]],y)$time,skyline(x[[i]],y)$population.size)
  colnames(tmp)<-c("time","population.size")
  breaks1<-data.frame()
  breaks1<-seq(0,lower_HPD,length.out=1000)
  df_total<-data.frame(matrix(ncol=1000,nrow=1))
  df_total[1,]<-data.matrix(cut(breaks1, breaks=c(-Inf,tmp$time),labels = c(tmp$population.size)))
  df_total[1,][is.na(df_total[1,])] <- tail(skyline(x[[i]],y)$population.size,n=1)
  return(df_total[1,])
}
#Apply the above function with multi-threads
library(foreach)
library(doParallel)
#Setup parallel backend to use many processors (whole cpus-2)
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
finalMatrix_upgma<-data.frame(matrix(ncol=1000,nrow=length(all_upgma_tre)))
lower_HPD_HMA<-get_lower_HPD(all_upgma_tre,1e-7) # set epsilon = 1e-7, as few trees coalescent intervals = 0 if set epsilon = 0
finalMatrix_upgma <- foreach(i=1:length(all_upgma_tre), .combine=rbind) %dopar% {
  library(ape)
  library(phangorn)
  library(phytools)
  get_fake_beast_para(i,all_upgma_tre,1e-7,lower_HPD_HMA) #calling the function
}
stopCluster(cl)
#Plot the following four results in one figure
layout(mat= matrix(1:4,2,2,byrow=TRUE))
#Plot the classic skyline plot with HMA data
library(devtools)
#devtools::install_github("laduplessis/bdskytools") #install bdskytools from github
library(bdskytools)
finalMatrix_upgma<-data.matrix(finalMatrix_upgma)
times_hma<-seq(0,lower_HPD_HMA,length.out=1000)
Re_hpd_hma <- getMatrixHPD(finalMatrix_upgma)
plotSkylinePretty(times_hma, Re_hpd_hma, type='smooth', xlab="Mutational units", ylab="Population size", main="HMA epsilon=0.0", col=4, yline=2, xline=2)
#Set the sliding window function to make the skyline plot smooth (100 by 20)
library(zoo)
get_fake_slide_matrix <- function(x,w,s) {
  n<-length(rollapply(zoo(c(1:1000)), width = w, by = s, FUN = mean, align = "left"))
  df_total<-data.frame(matrix(ncol=n,nrow=nrow(x)))
  for (i in 1:nrow(x))
  {
    df_total[i,] <- rollapply(x[i,], width = w, by = s, FUN = mean, align = "left")
  }
  df_total <- data.matrix(df_total)
  return(df_total)
}
Re_hpd_hma_slide <- getMatrixHPD(get_fake_slide_matrix(finalMatrix_upgma,100,20))
times_hma_slide<-seq(0,lower_HPD_HMA,length.out=length(rollapply(zoo(c(1:1000)), width = 100, by = 20, FUN = mean, align = "left")))
plotSkylinePretty(times_hma_slide, Re_hpd_hma_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="HMA sliding window 100 by 20", col=4, yline=2, xline=2)

#Load the Beast trees file with burnin first 10% trees
all_beast_tree<-read.nexus(file="Data/Tree_file/16KAN-CL-mafft-beast.trees")
all_beast_tree<-all_beast_tree[c(floor(0.1 * length(all_beast_tree)):length(all_beast_tree))]
#Apply the previous function for generating the classic skyline plot matrix with multi-threads 
library(foreach)
library(doParallel)
#Setup parallel backend to use many processors (whole cpus-2)
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
finalMatrix_beast<-data.frame(matrix(ncol=1000,nrow=length(all_beast_tree)))
lower_HPD_beast<-get_lower_HPD(all_beast_tree,0)
finalMatrix_beast <- foreach(i=1:length(all_beast_tree), .combine=rbind) %dopar% {
  library(ape)
  library(phangorn)
  library(phytools)
  get_fake_beast_para(i,all_beast_tree,0,lower_HPD_beast) #calling the function
}
stopCluster(cl)
#Plot the classic skyline plot with nucleotide data
finalMatrix_beast<-data.matrix(finalMatrix_beast)
Re_hpd_beast <- getMatrixHPD(finalMatrix_beast)
times_beast<-seq(0,lower_HPD_beast,length.out=1000)
plotSkylinePretty(times_beast, Re_hpd_beast, type='smooth', xlab="Mutational units", ylab="Population size", main="Bayesian epsilon=0.0", col=4, yline=2, xline=2)
#Use the slide window function to make the plot smooth (100 by 20)
Re_hpd_beast_slide <- getMatrixHPD(get_fake_slide_matrix(finalMatrix_beast,100,20))
times_beast_slide<-seq(0,lower_HPD_beast,length.out=length(rollapply(zoo(c(1:1000)), width = 100, by = 20, FUN = mean, align = "left")))
plotSkylinePretty(times_beast_slide, Re_hpd_beast_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="Bayesian sliding window 100 by 20", col=4, yline=2, xline=2)
#Plot two results without sliding window in one figure
layout(mat= matrix(1:2,1,2,byrow=TRUE))
plotSkylinePretty(times_hma, Re_hpd_hma, type='smooth', xlab="Mutational units", ylab="Population size", main="A) Skyline plot based on HMA data", ylims = c(0,0.6), col=4, yline=2, xline=2)
plotSkylinePretty(times_beast, Re_hpd_beast, type='smooth', xlab="Mutational units", ylab="Population size", main="B) Skyline plot based on nucleotide data", ylims = c(0,0.6), col=4, yline=2, xline=2)
#Plot two results with sliding window in one figure
layout(mat= matrix(1:2,1,2,byrow=TRUE))
plotSkylinePretty(times_hma_slide, Re_hpd_hma_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="A) Skyline plot based on HMA data", col=4, ylims = c(0,0.4), yline=2, xline=2)
plotSkylinePretty(times_beast_slide, Re_hpd_beast_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="B) Skyline plot based on nucleotide data", col=4, ylims = c(0,0.4), yline=2, xline=2)