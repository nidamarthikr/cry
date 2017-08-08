#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads and SHELXD log files
#'
#' @param filename A character string. The path to a valid log file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the log header.
#' @return A named list. Each name correspond to a valid field in the log
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/shelxd/shelxd.log"
#' ltmp <- .readSHELXDlog(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export

readSHELXlog <- function(filename)
{
  # data<-readLines(filename)
  header <- scan(filename, nlines = 3, what = character(), quiet = TRUE)
  logfile <- .extract_RE(filename, message = TRUE)
  # Convert a character data frame ti a numeric one.

  ## Specifi case for shelxe, becasue there are four datasets
  if (header[3] == grep("SHELXE", header, value = TRUE)) {
    logfile2 <- list()
    CYCLE <- logfile$CYCLE
    logfile2$CYCLE <- as.data.frame(sapply(CYCLE,as.numeric))
    FOM_mapCC <- logfile$FOM_mapCC
    logfile2$FOM_mapCC <- as.data.frame(sapply(FOM_mapCC,as.numeric))
    Site1 <- logfile$Site1
    logfile2$Site1 <- as.data.frame(sapply(Site1,as.numeric))
    Site2 <- logfile$Site2
    logfile2$Site2 <- as.data.frame(sapply(Site2,as.numeric))
  }

  else logfile2 <- as.data.frame(sapply(logfile, as.numeric))
  return(logfile2)
}

