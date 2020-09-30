#' Parametric bootstrap sampling using the predict interval
#'
#' This function was used to bootstrap sampling for HMA data using the predict interval.
#' Uses neighbour-joining (nj) algorithm as implemented in ape package.
#' Uses upgma method as implemented in phangorn package.
#' @param pair_data Assume by default that the data matrix is the pairwise distances matrix, which contain two column named ober (HMA data) and real (DNA data).
#' @param dframe A data frame for the samples.
#' @param method only allow "nj" or "upgma" method to build the phylogeny.
#' @return the phylogenetic tree of the parametric bootstrap.
#'
#' @export
para_boot<- function(pair_data, dframe, method) {
  stopifnot(is.data.frame(dframe))
  lm.out <- stats::lm(real ~ ober, data=pair_data)
  conf_interval <- suppressWarnings(stats::predict(lm.out, interval="predict"))
  sd<- vector()
  n_pober<- vector()
  n_pober_pos<- vector()
  for(i in 1:nrow(pair_data))
  {
    sd[i]<-(conf_interval[i,3]-conf_interval[i,1])/1.96
    n_pober[i] <- stats::rnorm(1,conf_interval[i,1],sd[i])
  }
  n_pober_pos<-pmax(n_pober,0)
  dframe[lower.tri(dframe)]<-n_pober_pos
  dframe<-t(dframe)
  dframe[lower.tri(dframe)]<-n_pober_pos
  if (method == "nj") return(ape::nj(dframe))
  else if (method == "upgma") return(phangorn::upgma(dframe))
  else stop("'Method' must be neighbour-joining (nj) or upgma")
  }
