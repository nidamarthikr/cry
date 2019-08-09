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
    gg <- ggplot(filenameDF, aes(1/(Res)^2, var))
    if(var == filenameDF$d_sig && length(var) >= 1) {
      gp <- gg + geom_point() + geom_line() + theme_bw() +
        xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Delta*F/sig(Delta*F)))
    }
    if(var == filenameDF$Chi_sq && length(var) >= 1){
      gp <- gg + geom_point() + geom_line() + theme_bw() +
        xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Chi^2))
    }
    if(var == filenameDF$I_sig && length(var) >= 1){
      gp <- gg + geom_point() + geom_line() + theme_bw() +
        xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(I/sig(I)))
    }
    if(var == filenameDF$Complete && length(var) >= 1) {
      gp <- gg + geom_point() + geom_line() + theme_bw() +
        xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Completeness))
    }
    if(var == filenameDF$CC1_2 && length(var) >= 1){
      gp <- gg + geom_point() + geom_line() + theme_bw() +
        xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(CC_(1/2)))
    }
    # gg + geom_point() + geom_line() + theme_bw() +
    #   xlab(expression(h^2 * (A^-2))) + ylab(expression(Delta*F/sig(Delta*F)))

    return(gp)


}
