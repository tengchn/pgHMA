% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Skyline_plot.R
\name{Skyline_plot}
\alias{Skyline_plot}
\title{Draw the Skyline plot with or without Y log axis.}
\usage{
Skyline_plot(
  times,
  skyline_matrix,
  log_axis = FALSE,
  col = "blue",
  fill = "gray85",
  xlims = NULL,
  ylims = NULL,
  xlab = "Time",
  ylab = "Effective population size",
  main = "",
  base_size = 12
)
}
\arguments{
\item{times}{The time points or mutational units (X-axis) used to draw the Skyline plot.}

\item{skyline_matrix}{The HDI matrix (including three rows of lower, median and upper bound) used to draw the Skyline plot.}

\item{log_axis}{Use the Y log axis or not. Default is false.}

\item{col}{Set the color for the median trend.}

\item{fill}{Set the color for the credible interval.}

\item{xlims}{scale limits for X-axis.}

\item{ylims}{scale limits for Y-axis.}

\item{xlab}{name for X-axis.}

\item{ylab}{name for Y-axis.}

\item{main}{name for the plot title.}

\item{base_size}{Set the base font size.}
}
\description{
This function was used to draw the Skyline plot from the HDI matrix with or without Y log axis. Use ggplot2 package.
Assuming each column is a different time point, and the dimension between times and skyline_matrix is match.
}
