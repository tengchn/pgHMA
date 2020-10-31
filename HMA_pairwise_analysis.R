library(ape)
library(phangorn)
library(phytools)
library(ggplot2)
library(devtools)
devtools::install_github("tengchn/pgHMA/pgHMAtools") #install pgHMAtools from github
library(pgHMAtools)
setwd("") ##change it to your own path

#1.Phylogenetic analysis based on HMA and nucleotide sequencing data
#Build a distance data frame for 16 samples 
df16<-data.frame(matrix(0,ncol=16,nrow=16))
x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
rownames(df16)<-x
colnames(df16)<-y
df16
#Import distance matrix of nucleotide data and HMA data
#pair_data<-read.csv("Data/HMA-16-ober-real.csv")
data(pair_data)
head(pair_data)
#Construct neighbour-joining(NJ) tree based on DNA data
df16[lower.tri(df16)]<-pair_data$real
df16<-t(df16)
df16[lower.tri(df16)]<-pair_data$real
df16_DNA<-df16
tr_real_16<-nj(df16_DNA)
plot(midpoint.root(tr_real_16))
add.scale.bar()
write.tree(tr_real_16,file = "Data/Tree_file/DNA_nj.tre")
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
plot_lm(pair_data) #call the function plot_lm from pgHMAtools
#Summary of the linear regression
summary(lm.out)
summary(conf_interval)

#3.Parametric bootstrap sampling using the predict interval
np<-10000
set.seed(12345)
df16<-data.frame(matrix(0,ncol=16,nrow=16))
x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
rownames(df16)<-x
colnames(df16)<-y
for (j in 1:np){
  tr_ober16 <- para_boot(pair_data,df16,"nj") #call the function para_boot from pgHMAtools
  tr_ober16$edge.length <- pmax(tr_ober16$edge.length,0) #Change the negative edge length value to 0
  write.tree(tr_ober16,file="Data/Tree_file/10000bp_tree_HMA_nj.tre",append = T)
}

#4.Skyline plot analyses
#Parametric bootstrap sampling using UPGMA for HMA data
np<-10000
set.seed(12345)
df16<-data.frame(matrix(0,ncol=16,nrow=16))
x<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
y<-c("KS1","KS2","KS3","KS4","KS5","KS6","KS8","KD1","KD2","KD3","KD4","KD5","KD6","KD7","KF1","KS7")
rownames(df16)<-x
colnames(df16)<-y
for (j in 1:np){
  tr_ober16 <- para_boot(pair_data,df16,"upgma") #call the function para_boot from pgHMAtools
  tr_ober16$edge.length <- pmax(tr_ober16$edge.length,0) #Change the negative edge length value to 0
  write.tree(tr_ober16,file="Data/Tree_file/10000bp_tree_HMA_upgma.tre",append = T)
}
#Load the UPGMA trees file
all_upgma_tre<-read.newick(file="Data/Tree_file/10000bp_tree_HMA_upgma.tre")
#Apply related functions with multi-threads
library(foreach)
library(doParallel)
#Setup parallel backend to use multi-threads (whole cpus-2)
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
lower_HDI_HMA<-get_lower_HDI(all_upgma_tre) #calling the function get_lower_HDI from pgHMAtools
finalMatrix_upgma <- foreach(i=1:length(all_upgma_tre), .combine=rbind) %dopar% {
  pgHMAtools::get_DH(i,all_upgma_tre,lower_HDI_HMA,1e-7) #calling the function get_DH from pgHMAtools, set epsilon = 1e-7, as few trees coalescent intervals = 0 if set epsilon = 0
}
stopCluster(cl)
#Plot the following four results in one figure
layout(mat= matrix(1:4,2,2,byrow=TRUE))
#Plot the skyline plot with HMA data
devtools::install_github("laduplessis/bdskytools") #install bdskytools from github
library(bdskytools)
finalMatrix_upgma<-data.matrix(finalMatrix_upgma)
times_hma<-seq(0,lower_HDI_HMA,length.out=1000)
Re_hdi_hma <- get_HDI_matrix(finalMatrix_upgma) #calling the function get_HDI_matrix from pgHMAtools
plotSkylinePretty(times_hma, Re_hdi_hma, type='smooth', xlab="Mutational units", ylab="Population size", main="HMA epsilon=0.0", col=4, yline=2, xline=2)
#Use the sliding window method to make the skyline plot smooth (100 by 20)
library(zoo)
Re_hdi_hma_slide <- get_slide_matrix(finalMatrix_upgma,100,20) #calling the function get_slide_matrix from pgHMAtools
times_hma_slide<-seq(0,lower_HDI_HMA,length.out=length(rollapply(zoo(c(1:1000)), width = 100, by = 20, FUN = mean, align = "left")))
plotSkylinePretty(times_hma_slide, Re_hdi_hma_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="HMA sliding window 100 by 20", col=4, yline=2, xline=2)

#Load the Beast trees file with burnin first 10% trees
all_beast_tree<-read.nexus(file="Data/Tree_file/16KAN-CL-mafft-beast.trees")
all_beast_tree<-all_beast_tree[c(floor(0.1 * length(all_beast_tree)):length(all_beast_tree))]
#Apply the function for generating the classic skyline plot matrix with multi-threads 
library(foreach)
library(doParallel)
#Setup parallel backend to use multi-threads (whole cpus-2)
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
lower_HPD_beast<-get_lower_HDI(all_beast_tree) #calling the function get_lower_HDI from pgHMAtools
finalMatrix_beast <- foreach(i=1:length(all_beast_tree), .combine=rbind) %dopar% {
  pgHMAtools::get_DH(i,all_beast_tree,lower_HPD_beast,0) #calling the function get_DH from pgHMAtools
}
stopCluster(cl)
#Plot the classic skyline plot with nucleotide data
finalMatrix_beast <- data.matrix(finalMatrix_beast)
Re_hpd_beast <- get_HDI_matrix(finalMatrix_beast) #calling the function get_HDI_matrix from pgHMAtools
times_beast<-seq(0,lower_HPD_beast,length.out=1000)
plotSkylinePretty(times_beast, Re_hpd_beast, type='smooth', xlab="Mutational units", ylab="Population size", main="Bayesian epsilon=0.0", col=4, yline=2, xline=2)
#Use the slide window function to make the plot smooth (100 by 20)
Re_hpd_beast_slide <- get_slide_matrix(finalMatrix_beast,100,20) #calling the function get_slide_matrix from pgHMAtools
times_beast_slide<-seq(0,lower_HPD_beast,length.out=length(rollapply(zoo(c(1:1000)), width = 100, by = 20, FUN = mean, align = "left")))
plotSkylinePretty(times_beast_slide, Re_hpd_beast_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="Bayesian sliding window 100 by 20", col=4, yline=2, xline=2)
#Plot two results without sliding window in one figure
layout(mat= matrix(1:2,1,2,byrow=TRUE))
plotSkylinePretty(times_hma, Re_hdi_hma, type='smooth', xlab="Mutational units", ylab="Population size", main="A) Skyline plot based on HMA data", ylims = c(0,0.6), col=4, yline=2, xline=2)
plotSkylinePretty(times_beast, Re_hpd_beast, type='smooth', xlab="Mutational units", ylab="Population size", main="B) Skyline plot based on nucleotide data", ylims = c(0,0.6), col=4, yline=2, xline=2)
#Plot two results with sliding window in one figure
layout(mat= matrix(1:2,1,2,byrow=TRUE))
plotSkylinePretty(times_hma_slide, Re_hdi_hma_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="A) Skyline plot based on HMA data", col=4, ylims = c(0,0.4), yline=2, xline=2)
plotSkylinePretty(times_beast_slide, Re_hpd_beast_slide, type='smooth', xlab="Mutational units", ylab="Population size", main="B) Skyline plot based on nucleotide data", col=4, ylims = c(0,0.4), yline=2, xline=2)
