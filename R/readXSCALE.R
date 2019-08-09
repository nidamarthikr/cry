#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads and output XSCALE.ahkl file with header and reflections
#'
#' @param filename A character string. The path to a valid mtz file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the mtz header.
#' @return A named list. Each name correspond to a valid field in the mtz
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/ahkl/file.ahkl"
#' ltmp <- readXSCALE(filename)
#' print(names(ltmp))
#' }
#' @export



#filename = "/Users/ritagiordano/Documents/ECM30/Dati_ECM/cl1.ahkl"


# read reflection
readXSCALE<-function(filename, messages = TRUE, ISIGMA = 1) {
  # Read XSCALE.ahkl file from XSCALE data processing
  # Determine the number of row present in the file
  #f_xscale<-file(filename,open="rb")
  xscale <- read.table(filename, comment.char = "!")
  #ahkl<-as.data.frame(ahkl)
  # XSCALE file contains 5 columns H,K,L Iobserved and sigma
  colnames(xscale)<-c("H","K","L","Iobs","sigma")
  ISIG <- xscale$Iobs/xscale$Sigma_Iobs
  xds <- cbind(xds, ISIG)
  if (ISIGMA != 1)
  {
    xds <- xds[xds$ISIG <= ISIGMA, ]
  }
  return(xscale)
}

#read header
# readXSCALEHeader <- function(filename, message = TRUE)
# {
#   f_xscale<-file(filename,open="r")
#   header <- .extract_RE_xds(filename, message)
#   return(header)
# }

