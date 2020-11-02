#' Draw the Skyline plot with or without Y log axis.
#'
#' This function was used to draw the Skyline plot from the HDI matrix with or without Y log axis. Use ggplot2 package.
#' Assuming each column is a different time point, and the dimension between times and skyline_matrix is match.
#'
#' @param times The time points or mutational units (X-axis) used to draw the Skyline plot.
#' @param skyline_matrix The HDI matrix (including three rows of lower, median and upper bound) used to draw the Skyline plot.
#' @param log_axis Use the Y log axis or not. Default is false.
#' @param xlims scale limits for X-axis.
#' @param ylims scale limits for Y-axis.
#' @param base_size Set the base font size.
#' @param col Set the color for the median trend.
#' @param fill Set the color for the credible interval.
#' @param xlab name for X-axis.
#' @param ylab name for Y-axis.
#' @param main name for the plot title.
#'
#' @export
Skyline_plot <- function(times, skyline_matrix, log_axis=FALSE, col="blue", fill="gray85", xlims=NULL, ylims=NULL, xlab="Time", ylab="Effective population size", main="",
                         base_size = 12) {

if (length(times) == ncol(skyline_matrix) && nrow(skyline_matrix) == 3) {
    if (log_axis == FALSE) {

      P<-ggplot2::ggplot() +
        ggplot2::geom_polygon(data.frame(x=c(times, rev(times)), y=c(skyline_matrix[1,], rev(skyline_matrix[3,]))), mapping = ggplot2::aes(x=x,y=y), fill=fill)+
        ggplot2::geom_line(data.frame(x=times, y=skyline_matrix[2,]), mapping=ggplot2::aes(x=x,y=y), col=col)+
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(), expand = c(0, 0), labels = scales::format_format(zero.print = "0", drop0trailing = TRUE))+
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(), labels = scales::format_format(zero.print = "0", drop0trailing = TRUE))+
        ggplot2::labs(x=xlab, y=ylab, title = main)+
        ggplot2::theme_classic(base_size = base_size)
      if (is.null(ylims)){
        P + ggplot2::coord_cartesian(xlim = xlims)
      }
      else {
        P + ggplot2::coord_cartesian(xlim = xlims, ylim = ylims, expand = F)
      }
    }
  else
    if (log_axis == TRUE) {

      P<-ggplot2::ggplot() +
        ggplot2::geom_polygon(data.frame(x=c(times, rev(times)), y=c(skyline_matrix[1,], rev(skyline_matrix[3,]))), mapping = ggplot2::aes(x=x,y=y), fill = fill)+
        ggplot2::geom_line(data.frame(x=times, y=skyline_matrix[2,]), mapping=ggplot2::aes(x=x,y=y), col=col)+
        ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^as.integer(x)), labels = scales::format_format(zero.print = "0", drop0trailing = TRUE, digits = 1, scientific = FALSE), expand = ggplot2::expansion(mult = c(0, 0.1)))+
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(), expand = c(0, 0), labels = scales::format_format(zero.print = "0", drop0trailing = TRUE))+
        ggplot2::annotation_logticks(sides = "l", outside = TRUE, mid = ggplot2::unit(0.1, "cm"), long = ggplot2::unit(0.15, "cm"))+
        ggplot2::labs(x=xlab, y=ylab, title = main)+
        ggplot2::theme_classic(base_size = base_size)
      if (is.null(ylims)){
        P + ggplot2::coord_cartesian(clip = "off", xlim = xlims)
      }
      else {
        P + ggplot2::coord_cartesian(clip = "off", xlim = xlims, ylim = ylims, expand = F)
      }
    }}
else
  if (length(times) != ncol(skyline_matrix))
    stop("The dimension between times and skyline_matrix is not match!")
else
    stop("Please check the parameters for Skyline plot!")
}
