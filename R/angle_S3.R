#
# This file is part of the cry package
#

#########
### S3 class "angle" and related methods
########

#' Constructor for an S3 object of class "angle"
#'
#' @param ang A real number, in degrees or radians depending on rad_flag.
#' @param rad_flag A logical flag. If FALSE, the value is meant to be in radians.
#' @return An object of class "angle" whose numerical value is always in degrees.
#' @examples
#' # Create an angle of 60 degrees
#' ang1 <- angle(60)
#' class(ang1)
#' @export
angle <- function(ang,rad_flag=TRUE) {
  if (!is.numeric(ang)) stop("ang must be numeric")
  if (!is.logical(rad_flag)) stop("rad_flag must be logical")

  # If rad_flag is FALSE the angle must be transformed in degrees
  if (!rad_flag) ang <- ang*180/pi

  # Now enrich and with the attribute "rad_flag"
  attr(ang,'rad_flag') <- TRUE

  # S3 class
  ang <- structure(ang,class = "angle")

  return(ang)
}


#' Print method for an object of class "angle".
#'
#' The value is displayed in degrees
#'
#' @param x An object of class "angle".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create an angle of 90 degrees using radians
#' ang1 <- angle(pi/2,FALSE)
#'
#' # Display its value
#' print(ang1)
#' @rdname print.angle
#' @export
print.angle <- function(x,...) {
  print(x[1])
  invisible(x)
}
