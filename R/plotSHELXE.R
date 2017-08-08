#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Plot SHELXD log files
#'
#' @param filenameDF_i A data frame containing the inverted hand.
#' @param filenameDF_o A data frame containing the original hand.
#' @return A graphics containg the solution founded by SHELXD.
#'
#' @examples
#' \dontrun{
#' plot <- plotSHELXE(filenameDF_i, filenameDF_o, message = TRUE)
#' plot
#' }
#' @export

plotSHELXE <-function(filenameDF_i, filenameDF_o, message = TRUE)
{
  plot <- .plotSHELXD(filenameDF_i, filenameDF_o, message = TRUE)
  return(plot)
}
