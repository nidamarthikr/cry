#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads XDS file with header and reflections.
#'
#' @param filename A character string. The path to a valid xds file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the xds header.
#' @param ISIGMA A numerical value, which the default value is 1
#' @return A named list. Each name correspond to a valid field in the xds
#'    header.
#' @description ISIGMA will cut the data following the wished I over sigma
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/HKL/file.HKL"
#' ltmp <- readXDS(filename)
#' print(names(ltmp))
#' print(ltmp$UNIT_CELL_CONSTANTS)
#' }
#' @export

#filename = "/Users/ritagiordano/Documents/ECM30/Files/XDS.HKL"

readXDS<-function(filename, message = TRUE, ISIGMA = 1) {
  # Read XDS_ASCII.HKL file from xds data processing.
  # Determine the number of row present in the file.
  # create and open a connection to the xds file.
  xds <- read.table(filename, comment.char = "!")
  colnames(xds)<-c("H","K","L","Iobs","Sigma_Iobs", "XD", "YD", "ZD", "RLP", "PEAK", "CORR", "PSI")
  ISIG <- xds$Iobs/xds$Sigma_Iobs
  xds <- cbind(xds, ISIG)
  if (ISIGMA != 1)
  {
    xds <- xds[xds$ISIG <= ISIGMA, ]
  }
  return(xds)
}

# readXDSHeader <- function(filename, message = TRUE)
# {
#   f_xds<-file(filename,open="r")
#   header <- .extract_RE_xds(filename, message)
#   return(header)
# }

## Function to write new XDS_ASCII.HKL
#' Write a new XDS file with header and reflections.
#'
#' @param filename A character string. The path to a valid xds file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the xds header.
#' @param path A numerical value, which the default value is 1
#' @return A named list. Each name correspond to a valid field in the xds
#'    header.
#'
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/HKL/file.HKL"
#' ltmp <- readXDS(filename)
#' print(names(ltmp))
#' print(ltmp$UNIT_CELL_CONSTANTS)
#' }
#' @export

# writeXDS_ASCII <- function(filename, path, append = FALSE, message = TRUE)
# {
#   f_xds <- file(filename, open = "rw")
#   xds <- write(f_xds, path, append = FALSE)
# }


