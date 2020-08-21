setwd("") ##change it to your own path

#1.Optimization
mixture_data<-read.csv("Data/20 mixtures for optimization.csv")
head(mixture_data)
str(mixture_data)
dna_freq<-mixture_data[,1]
hma_freq<-mixture_data[,2]
HD<-mixture_data[,3]

#(1).Use linear regression model 
obj1<-function(df){
  a<-df[1]
  b<-df[2]
  opti_freq<-hma_freq*(a+b*HD)
  c1<-sqrt(sum((dna_freq-opti_freq)^2))
  #print(opti_freq)
  return(c1)
}
print("Lineal regression")
opti_lineal<-optim(par=c(1,1),fn=obj1, method = "L-BFGS-B",control = list(maxit=1e9, fnscale=1))
#(2).Use exponential distribution model
obj2<-function(df){
  a<-df[1]
  b<-df[2]
  opti_freq<-hma_freq*(a*exp(b*HD))
  c2<-sqrt(sum((dna_freq-opti_freq)^2))
  #print(opti_freq)
  return(c2)
}
print("Exponential")
opti_expon<-optim(par=c(1,1),fn=obj2, method = "L-BFGS-B",control = list(maxit=1e9, fnscale=1))
#(3).Use log-normal distribution model
obj3<-function(df){
  a<-df[1]
  b<-df[2]
  c<-df[3]
  opti_freq<-c*hma_freq*(dlnorm(x = HD, meanlog = b, sdlog = a))
  c3<-sqrt(sum((dna_freq-opti_freq)^2))
  #print(opti_freq)
  return(c3)
}
print("Log-normal")
opti_logn<-optim(par=c(1,1,1),fn=obj3, method = "L-BFGS-B",control = list(maxit=1e9, fnscale=1))
#(4).Use gamma distribution model
obj4<-function(df){
  a<-df[1]
  b<-df[2]
  c<-df[3]
  opti_freq<-c*hma_freq*(dgamma(x = HD, shape = a, scale = b))
  c4<-sqrt(sum((dna_freq-opti_freq)^2))
  #print(opti_freq)
  return(c4)
}
print("Gamma")
opti_gamma<-optim(par=c(1,1,1),fn=obj4, method = "L-BFGS-B",control = list(maxit=1e9, fnscale=1),lower = c(1,1e-10,-Inf))
#Sort by value
df_opti<-data.frame(matrix(ncol=2,nrow=4))
colnames(df_opti)<-c("Methods","Optimize_value")
df_opti[,1]<-c("Linear","Exponential","Log-normal","Gamma")
df_opti[,2]<-c(opti_lineal$value,opti_expon$value,opti_logn$value,opti_gamma$value)
df_opti[order(df_opti$Optimize_value),] #Top one is the best fit model


#2.Plot the regression figures for each group with before and after optimization
library(ggplot2)
library(gridExtra)
Mismatch_data <- read.csv("Data/Mismatch distribution/Mismatch_distribution_20 groups.csv", header = TRUE)
head(Mismatch_data)
str(Mismatch_data)
#Plots for 20 mixtures with the relationship between nucleotide frequency and HMA frequency
Mismatch_lm <- list()
for(i in 1:20){
  Mismatch_lm[[i]] <- summary(lm(Nucleotide.frequence ~ HMA.frequence, data = Mismatch_data[Mismatch_data$Group.number==i,]))
} 
par(mfrow = c(4, 5))
for(i in 1:20){
  plot(Mismatch_data$HMA.frequence[Mismatch_data$Group.number==i], Mismatch_data$Nucleotide.frequence[Mismatch_data$Group.number==i],
       col = "blue", type = "p",
       xlab = "Heteroduplex frequency",
       ylab = "Nucleotide frequency", cex.lab =1.2,
       main = paste(Mismatch_data$Group.ID[Mismatch_data$Group.number==i][1],"\n","R²=",format(Mismatch_lm[[i]][["r.squared"]],digits = 3)))
  intercept <- Mismatch_lm[[i]][["coefficients"]][[1,1]]
  slope <- Mismatch_lm[[i]][["coefficients"]][[2,1]]
  abline(a = intercept, b = slope, col = "red", lwd = 3)
}
#Plots for 20 mixtures with the relationship between nucleotide frequency and optimized heteroduplex frequency
Mismatch_optimize_lm <- list()
for(i in 1:20){
  Mismatch_optimize_lm[[i]] <- summary(lm(Nucleotide.frequence ~ Optimized.heteroduplex.frequence.after.normolization, data = Mismatch_data[Mismatch_data$Group.number==i,]))
}
par(mfrow = c(4, 5))
for(i in 1:20){
  plot(Mismatch_data$Optimized.heteroduplex.frequence.after.normolization[Mismatch_data$Group.number==i], Mismatch_data$Nucleotide.frequence[Mismatch_data$Group.number==i],
       col = "blue", type = "p",
       xlab = "Optimized HMA frequency",
       ylab = "Nucleotide frequency", cex.lab =1.2,
       main = paste(Mismatch_data$Group.ID[Mismatch_data$Group.number==i][1],"\n","R²=",format(Mismatch_optimize_lm[[i]][["r.squared"]],digits = 3)))
  intercept <- Mismatch_optimize_lm[[i]][["coefficients"]][[1,1]]
  slope <- Mismatch_optimize_lm[[i]][["coefficients"]][[2,1]]
  abline(a = intercept, b = slope, col = "red", lwd = 3)
}
# Plot all 20 mixture groups together in one figure
lm_total_1 <- summary(lm(Nucleotide.frequence ~ HMA.frequence, data = Mismatch_data)) #linear regression between nucleotide frequency and HMA frequency
lm_total_2 <- summary(lm(Nucleotide.frequence ~ Optimized.heteroduplex.frequence.after.normolization, data = Mismatch_data)) #linear regression between nucleotide frequency and optimized heteroduplex frequency
P1<-ggplot(Mismatch_data, aes(x=HMA.frequence, y=Nucleotide.frequence))+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = -0.1))+
  ggtitle("A)")+
  xlab ("Heteroduplex frequency")+
  ylab ("Nucleotide frequency")+
  geom_point(col="royalblue")+
  annotate(label = sprintf("y = %.3f x + %.3f\nR² = %.3f", lm_total_1$coefficients[2,1], lm_total_1$coefficients[1,1], lm_total_1$r.squared), geom = "text", x = 68, y = 105, size = 4)+
  geom_smooth(method="lm", formula = y ~ x, col="red")
P2<-ggplot(Mismatch_data, aes(x=Optimized.heteroduplex.frequence.after.normolization, y=Nucleotide.frequence))+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = -0.1))+
  ggtitle("B)")+
  xlab ("Optimized heteroduplex frequency")+
  ylab ("Nucleotide frequency")+
  geom_point(col="royalblue")+
  geom_smooth(method="lm", formula = y ~ x, col='red')+
  annotate(label = sprintf("y = %.3f x + %.3f\nR² = %.3f", lm_total_2$coefficients[2,1], lm_total_2$coefficients[1,1], lm_total_2$r.squared), geom = "text", x = 68, y = 105, size = 4)
grid.arrange(P1, P2, ncol=2)


#3.Cost plot
#Total cost
sanger <- function(x){x*7.5*2} #dual direction sanger sequencing cost (A$7.5*2)
hma <- function(x){factorial(x)/factorial(2)/factorial(x-2)*0.45} #LabChip per reaction cost A$0.45 
P1 <- ggplot(data.frame(x=c(2, 100)), aes(x=x)) + stat_function(fun=hma, geom="line", aes(colour = "HMA"), size=2.5) + stat_function(fun=sanger, geom="line", aes(colour = "Sanger"), size=2.5) + xlab("Number of samples") + ylab("Total cost (A$)")+ scale_colour_manual("Methods",values = c("red", "blue"))
#Cost per sample
sanger_per_indiv <- function(x){x*7.5*2/x}
hma_per_indiv <- function(x){factorial(x)/factorial(2)/factorial(x-2)*0.45/x}
P2<-ggplot(data.frame(x=c(2, 100)), aes(x=x)) + stat_function(fun=hma_per_indiv, geom="line", aes(colour = "HMA"), size=2.5) + stat_function(fun=sanger_per_indiv, geom="line", aes(colour = "Sanger"),size=2.5) + xlab("Number of samples") + ylab("Cost per sample (A$)") + scale_colour_manual("Methods",values = c("red", "blue"))
#Plot together
grid.arrange(P1+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text = element_text(size = 16),legend.title = element_text(size=16, face="bold")), P2+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text = element_text(size = 16),legend.title = element_text(size=16, face="bold")), ncol=2)
#Calculate the number of samples saving most by HMA 
obj_saving<-function(a){
  abs(a*15-factorial(a)/factorial(2)/factorial(a-2)*0.45)
}
HMA_saving<-optim(par=2,fn=obj_saving, method = "Brent",lower = 2, upper=100,control = list(maxit=1e9, fnscale=-1)) #find the maximum
cat("HMA save the most (A$",round(HMA_saving$value),") when apply", round(HMA_saving$par), "samples, compare with sanger sequencing.")
#Calculate the number of samples have similar cost between HMA and Sanger sequencing
HMA_similiar<-optim(par=2,fn=obj_saving, method = "Brent",lower = 2, upper=100,control = list(maxit=1e9, fnscale=1)) #find the minimum
cat("The cost of HMA cheaper than the commonly applied sanger sequencing strategy for less than", round(HMA_similiar$par), "samples")