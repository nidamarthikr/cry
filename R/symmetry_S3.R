#
# This file is part of the cry package
#



#' Constructor for an S3 object of class "cryst_symm".
#'
#' This represents a crystallographic space group.
#'
#' The constructor can be used with less than the full set of two input parameters,
#' to create an object of class cryst_symm corresponding to space group P 1. If the input
#' string is not recognised, a warning is raised and space group P 1 is assigned.
#'
#' @param SG A character string or an integer identifying the space group. There are 230
#'  used space group in crystallography and each one corresponds to a unique and so-called
#'  extended Hermann-Mauguin symbol. An example is space group number 19, identified by the
#'  extended Hermann-Mauguin symbol "P 21 21 21". Several formats are possible and some of them
#'  are now rarely used. Attempts are made to transform the input into a correct Hermann-Mauguin
#'  symbol, but if all fails, a warning is raised and the space group P 1 is assigned.
#' @param set An integer defining which setting of many possible for the given space group. Some
#'  crystallographic space groups can be implemented with small variants known as "settings".
#' @return An object of class "cryst_symm". It is a named list of length 4. The names are, "SG",
#'  "PG", "T" and "C".
#' \itemize{
#'  \item{1) SG. This is a string containing the correct extended Hermann-Mauguin sybol.}
#'  \item{2) PG. This is a list whose elements are all the \eqn{3\times 3} matrices forming
#'               the point-group part of the symmetry transformation.}
#'  \item{3) T. This is a list whose elements are all the \eqn{3\times 1} vectors forming
#'               the translational part of the symmetry transformation.}
#'  \item{4) T. This is a list whose elements are all the \eqn{3\times 1} vectors forming
#'               the centering of the unit cell.}
#' }
#'
#' @examples
#' # The symplest symmetry: P 1
#' crsym <- cryst_symm("P 1")
#' print(crsym)
#'
#' # The second space group: P -1
#' crsym <- cryst_symm(2)
#' print(crsym)
#'
#' # Something more complicated
#' crsym <- cryst_symm("P 21 21 21")
#' print(crsym)
#'
#' @export
cryst_symm <- function(SG=NULL,set=NULL) {
  # If the input parameters is NULL use default value (P1)
  if (is.null(SG) & is.null(set)) {
    SG <- "P 1"

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }

  # If SG is not NULL, set is needed
  if (is.null(set)) set <- 1

  # Try commonly-used names, if SG is a character string
  if (is.character(SG)) {
    SG <- findHM(SG)
    if (is.null(SG)) {
     warning("This space group does not exist. Use P 1.")

     SG <- "P 1"

     # Extract all operators
     ltmp <- syminfo_to_matrix_list(SG)

     # S3 class
     crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

     return(crsym)
    }
  }

  # Whatever the input, turn it into an xHM symbol
  tmp <- translate_SG(SG,"number","xHM",set)
  if (tmp$ans) {
    SG <- tmp$msg

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }
  tmp <- translate_SG(SG,"ccp4","xHM",set)
  if (tmp$ans) {
    SG <- tmp$msg

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }
  tmp <- translate_SG(SG,"xHM","xHM",set)
  if (tmp$ans) {
    SG <- tmp$msg

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }
  tmp <- translate_SG(SG,"Hall","xHM",set)
  if (tmp$ans) {
    SG <- tmp$msg

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }
  tmp <- translate_SG(SG,"old","xHM",set)
  if (tmp$ans) {
    SG <- tmp$msg

    # Extract all operators
    ltmp <- syminfo_to_matrix_list(SG)

    # S3 class
    crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

    return(crsym)
  }

  # The input does not correspond to anything known
  warning("This space group does not exist. Use P 1.")

  SG <- "P 1"

  # Extract all operators
  ltmp <- syminfo_to_matrix_list(SG)

  # S3 class
  crsym <- structure(list(SG=SG,PG=ltmp$PG,T=ltmp$T,C=ltmp$C),class = "cryst_symm")

  return(crsym)
}


#' Print method for an object of class "cryst_symm".
#'
#' xxx
#'
#' @param x An object of class "cryst_symm".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create an object of P 2 symmetry
#' crsym <- cryst_symm("P 2")
#'
#' # Display its value
#' print(crsym)
#'
#' @rdname print.cryst_symm
#' @export
print.cryst_symm <- function(x,...) {
  cat("This is an object of class 'cryst_symm'\n")
  msg <- paste0("The space group represented is ",x$SG,".\n")
  cat(msg)
  tmp <- translate_SG(x$SG,"xHM","number")
  msg <- paste0("The number associated to this group in the International Tables is: ",tmp$msg,".\n")
  cat(msg)
  stmp <- crystal_system(tmp$msg)
  msg <- paste0("Its crystal system is ",stmp,".\n")
  cat(msg)

  # Crystal symmetry
  # Symmetry operations
  tmp <- full_symm_strings(x$SG)
  cat("\n")
  cat("List of symmetry operations:\n")
  for (sline in tmp) {
    lne <- trimws(sline,"right")
    msg <- substr(lne,6,nchar(lne))
    cat(paste0("  ",msg,"\n"))
  }

  invisible(x)
}
