#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Plot SHELXD log files
#'
#' @param filenameDF A data frame containing CCall, CCweak
#' @return A graphics containg the solution founded by SHELXD.
#'
#' @examples
#' \dontrun{
#' plot <- plotSHELXD(filenameDF)
#' plot
#' }
#' @export

plotSHELXD <-function(filenameDF)
{
  ggplot(filenameDF,aes(CCall,CCweak)) + geom_point() + theme_bw() +
    xlab('CCall') + ylab('CCweak')
}
