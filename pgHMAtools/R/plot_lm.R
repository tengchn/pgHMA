#' Plot the linear regression with prediction interval
#'
#' This function was used to perform the correlation analysis between HMA and nucleotide sequencing data.
#' Uses ggplot2 package.
#' @param pair_data Assume by default that the data matrix is the pairwise distances matrix, which contain two column named ober (HMA data) and real (DNA data).
#'
#' @examples
#' data(pair_data)
#' plot_lm(pair_data)
#'
#' @export
#'
plot_lm <- function(pair_data) {
  lm.out <- stats::lm(real ~ ober, data=pair_data)
  conf_interval <- stats::predict(lm.out, interval="predict")
  ggplot2::ggplot(pair_data, ggplot2::aes(x=pair_data$ober, y=pair_data$real))+
    ggplot2::theme_bw(base_size = 15)+
    ggplot2::xlab(expression(paste("Heteroduplex mobility distance (",d[H], ")")))+
    ggplot2::ylab("Nucleotide distance")+
    ggplot2::geom_point(size = 3.5, shape = 21, bg = 'gray')+
    ggplot2::annotate(label = paste('y ==',round(lm.out$coefficients[2],3),'~ d[H] + ',round(lm.out$coefficients[1],3),"*plain(\",\")",'~~ R^2 ==', round(summary(lm.out)$r.squared,3)), geom = "text", x = 0.15, y = 0.18, size = 7, parse = TRUE)+
    ggplot2::geom_smooth(method="lm", size=1.5)+
    ggplot2::geom_line(ggplot2::aes(x=sort(pair_data$ober), y=sort(conf_interval[,2])), col = "red",size=1)+
    ggplot2::geom_line(ggplot2::aes(x=sort(pair_data$ober), y=sort(conf_interval[,3])), col = "red",size=1)
}
