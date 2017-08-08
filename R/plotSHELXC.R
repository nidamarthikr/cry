#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Plot SHELXC log files
#'
#' @param filenameDF A data frame containing the variable from shelxc.
#' @param var the variable to be plotted vs the resolution
#' @return A graphics containg the solution founded by SHELXD.
#'
#' @examples
#' \dontrun{
#' plot <- plotSHELXC(filenameDF, var, message = TRUE)
#' plot
#' }
#' @export

plotSHELXC <-function(filenameDF, var, message = TRUE)
{
  plot <- .plotSHELXC(filenameDF, var, message = TRUE)
  return(plot)
}
