#
# This file is part of the cry package
#



#' Constructor for an S3 object of class "CIF".
#'
#' This represents a CIF (Crystallographic Information Framework) data file.
#'
#' Crystallographic data can be collected in CIF files in the most complete and
#' structured way. An object of this class is created when the input is the name
#' of a valid CIF file. The output is a named list with number of fields and names
#' dependent on the CIF information included in the specific file.
#'
#' @param filename A string character. The name of the CIF file containing the data.
#' @return An object of class "CIF". It's a named list whose number and name of fields
#'  depends on the specific information content of the file.
#' @examples
#' # Read a CIF file
#'
#' @export
CIF <- function(filename=NULL) {
  # If one of the input parameters is NULL use default values
  #if (is.null(filename)) {
  #  #filename <- "name of pre-loaded file"

  #  # S3 class
  #  ltmp <- readCIF(filename)
  #  CIFobj <- structure(ltmp,class = "CIF")

  #  return(uc)
  ##}

  # General case
  # Check file exists
  if (is.character(filename)) {
    ltmp <- readCIF(filename)
  } else {
    warning("Input must be a valid file name.")

    return(NULL)
  }

  # S3 class
  CIFobj <- structure(ltmp,class = "CIF")

  return(CIFobj)
}


#' Print method for an object of class "CIF".
#'
#' CIF files include data related to crystallography in a popular framework supported by
#' IUCr, known as Crystallographic Information Framework. The content included varies according
#' to the type of data stored and it is organised in a list of objects, each one with a given
#' standard name. A few, commonly-used names are HEADER, SYMM, REFL, COOR, etc.
#'
#' @param x An object of class "CIF".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create a CIF object from a data file contained in the cry package
#' # cobj <- CIF("directory/file.cif")
#'
#' # Display its value
#' #print(cobj)
#' @rdname print.CIF
#' @export
print.CIF <- function(x,...) {
  cat("This is an object of class 'CIF'\n")
  cat("\nIts content includes:\n")

  # Find out how many objects are in the list
  tmp <- names(x)
  nn <- length(tmp)%/%10
  nmod <- length(tmp)%%10
  for (i in 1:nn) {
    msg <- ""
    istart <- 10*(i-1)+1
    iend <- 10*i
    cat(tmp[istart:iend],sep="  ")
  }
  if (nmod > 0) {
    cat(tmp((iend+1):length(tmp)),sep="  ")
  }

  # Cell sides
  a <- x$HEADER$CELL$A$VAL
  b <- x$HEADER$CELL$B$VAL
  c <- x$HEADER$CELL$C$VAL

  # Cell angles
  aa <- x$HEADER$CELL$ALPHA$VAL
  bb <- x$HEADER$CELL$BETA$VAL
  cc <- x$HEADER$CELL$GAMMA$VAL

  # Output unit cell info
  cat("\n\nThe unit cell of the crystal represented has sides:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    angstroms.\n",a,b,c)
  cat(msg)
  cat("and angles:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    degrees.\n",aa,bb,cc)
  cat(msg)
  sgr <- x$HEADER$HM
  msg <- paste0("\nThe space group represented is ",sgr,".\n")
  cat(msg)
  tmp <- translate_SG(sgr,"xHM","number")
  msg <- paste0("The number associated to this group in the International Tables is: ",tmp$msg,".\n")
  cat(msg)
  stmp <- crystal_system(tmp$msg)
  msg <- paste0("Its crystal system is ",stmp,".\n")
  cat(msg)

  invisible(x)
}
