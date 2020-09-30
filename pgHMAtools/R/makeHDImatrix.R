#' Get Lower HDI value of the phylogenetic tree
#'
#' This function was used to calculate the lower HDI value of the phylogenetic tree.
#' Uses hdi (Highest Density Interval) algorithm as implemented in HDInterval package.
#' Uses coalescent.intervals as implemented in ape package.
#' @param tree either an ultra-metric phylogenetic tree (i.e. an object of class "phylo") or a vector of interval lengths.
#' @param ci The value of credible interval.
#' @return lower limits of the HDI for all trees in the tree file.
#'
#' @export
get_lower_HDI <- function(tree,ci=0.95) {
  all_tree_CI<-list()
  for(i in 1:length(tree)){
    all_tree_CI[[i]]<-ape::coalescent.intervals(tree[[i]])$total.depth
  }
  lower_HDI<-HDInterval::hdi(do.call(rbind, all_tree_CI),ci)[1]
  return(lower_HDI)
}


#' Get the demographic history from each genealogy
#'
#' This function was used to obtain the demographic history from the genealogy
#' Uses classic skyline implemented in ape package to find the time at the end of each coalescent interval and effective population size of each interval.
#' @param i which number of the phylogenetic tree need to be analyzed.
#' @param tree either an ultra-metric phylogenetic tree (i.e. an object of class "phylo") or a vector of interval lengths.
#' @param epsilon collapsing parameter that controls the amount of smoothing.
#' @param lower_HDI Lower HDI value of the phylogenetic tree.
#' @param n A number of intervals used to build the matrix.
#' @return Normalized coalescent interval and effective population size for the phylogenetic tree.
#'
#' @export
get_DH <- function(i,tree,lower_HDI,epsilon=0,n=1000) {
  tmp <- data.frame(ape::skyline(tree[[i]],epsilon)$time,ape::skyline(tree[[i]],epsilon)$population.size)
  colnames(tmp)<-c("time","population.size")
  breaks1<-data.frame()
  breaks1<-seq(0,lower_HDI,length.out=n)
  df_total<-data.frame(matrix(ncol=n,nrow=1))
  df_total[1,]<-data.matrix(cut(breaks1, breaks=c(-Inf,tmp$time),labels = c(tmp$population.size)))
  df_total[1,][is.na(df_total[1,])] <- utils::tail(ape::skyline(tree[[i]],epsilon)$population.size,n=1)
  return(df_total[1,])
}


#' Get a HDI matrix for the phylogenetic tree
#'
#' This function was used to obtain a HDI matrix of the phylogenetic tree.
#' Uses hdi (Highest Density Interval) algorithm as implemented in HDInterval package.
#'
#' @param dat Data matrix used for extracting highest density intervals, and assume that each row represents a sample and each column show the coalescent interval.
#' @param margin Default by column.
#' @param ci The value of credible interval.
#' @return A HDI matrix with each column is the value of lower, median, and upper bound.
#'
#' @export
get_HDI_matrix <- function(dat, margin=2, ci=0.95) {
  if (is.matrix(dat)==FALSE) {dat<-data.matrix(dat)}
  getHDI <- function(x,ci) {c(HDInterval::hdi(x,ci)[1], stats::median(x), HDInterval::hdi(x,ci)[2])}
  return(apply(dat, margin, getHDI, ci))
}


#' Get a smooth HDI matrix by sliding window method
#'
#' This function was used to obtain a smooth HDI matrix by sliding window.
#' Uses zoo algorithm as implemented in zoo package.
#'
#' @param x Data matrix of each row represents a sample and each column show the coalescent interval.
#' @param w window size of the sliding window.
#' @param s step size of the sliding window.
#' @return A smooth HDI matrix with each column is the value of lower, median, and upper bound.
#'
#' @export
get_slide_matrix <- function(x,w=100,s=20) {
  n<-length(zoo::rollapply(zoo::zoo(c(1:ncol(x))), width = w, by = s, FUN = mean, align = "left"))
  df_total<-data.frame(matrix(ncol=n,nrow=nrow(x)))
  for (i in 1:nrow(x))
  {
    df_total[i,] <- zoo::rollapply(x[i,], width = w, by = s, FUN = mean, align = "left")
  }
  df_total <- data.matrix(df_total)
  return(get_HDI_matrix(df_total))
}
