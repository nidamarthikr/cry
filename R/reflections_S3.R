#
# This file is part of the cry package
#



#' Constructor for an S3 object of class "merged_reflections".
#'
#' This represents scaled and merged x-ray data from one crystal.
#'
#' If the constructor is used without arguments, the default
#' object created will be create reflections for a cubic crystal
#' with cell of side 10 angstroms, and symmetry P 2 3, up to 5
#' angstroms resolution. The only available columns
#' will be of dtype "H", named H, K, L (the Miller indices), and
#' of dtype "S", inverse resoluton, named S.
#'
#' The possible dtypes are:
#' \describe{
#'   \item{H}{Miller index}
#'   \item{S}{Inverse resolution (1/angstroms)}
#'   \item{J}{Reflectionn intensity}
#'   \item{F}{Amplitude of a structure factor}
#'   \item{D}{Anomalous difference}
#'   \item{Q}{Standard deviation of J, F, D}
#'   \item{G}{Amplitude associated with a Friedel pair (F(+), F(-))}
#'   \item{L}{Standard deviation of G}
#'   \item{K}{Intensity associated with G (I(+), I(-))}
#'   \item{M}{Standard deviation of K}
#'   \item{E}{Amplitude of the normalised structure factors}
#'   \item{P}{Phase angle (in degrees)}
#'   \item{W}{Weight of some sort}
#'   \item{A}{Phase probability coefficients (Hendrickson-Lattman)}
#'   \item{B}{Batch number (from raw data)}
#'   \item{I}{Any other integer}
#'   \item{R}{Any other real}
#' }
#' More vaues can become available in a future release.
#'
#' @param ruc An object of class "rec_unit_cell" (which represents
#'        a reciprocal unit cell).
#' @param csym An object of class "cryst_symm" (which represents
#'        a crystallographic symmetry group).
#' @param records A data frame containing all reflections coming
#'        from the x-ray data collection on the crystal. This
#'        data frame must include at least the three Miller
#'        indices, H, K, L (of dtype "H").
#' @param dtypes A character vector whose length is the same as
#'        the number of columns in 'records'. One character (a
#'        capital letter) is associated with each type of data.
#'        For example, a Miller index is of dtype "H"; a structure
#'        amplitude is of dtype "F"; an anomalous difference is of
#'        dtype "D"; etc (see details later).
#' @return An object of class "merged_reflections". It is a named
#'         list of length 4 whose names are:
#'         \describe{
#'           \item{ruc}{An object of class "rec_unit_cell".}
#'           \item{csym}{An object of class "cryst_symm".}
#'           \item{records}{A data frame containing the data.}
#'           \item{dtypes}{A character vector containing the
#'                 type of data (Miller indices, structure
#'                 factors, etc).}
#'         }
#'
#' @examples
#' # Create an orthorombic (default) cell
#' uc <- unit_cell(10,30,15)
#'
#' # Create the related reciprocal cell
#' ruc <- create_rec_unit_cell(uc)
#'
#' # Create symmetry (P n c 2)
#' csym <- cryst_symm(30)
#'
#' # Create a few records (these include syst. absences)
#' records <- expand.grid(H=-2:2,K=-2:2,L=-2:2)
#' print(length(records[,1]))
#'
#' # dtypes are all H
#' dtypes <- c("H","H","H")
#'
#' # Create merged_reflections object with H, K, L
#' # Systematic absences have been eliminated
#' mrefs <- merged_reflections(ruc,csym,records,dtypes)
#' print(length(mrefs$records[,1]))
#'
#' @export
merged_reflections <- function(ruc=NULL,csym=NULL,
                               records=NULL,dtypes=NULL) {
  # If all input parameters are NULL use default values
  # Cubic cell of side 10. Data to 5 angstroms resolution
  if (is.null(ruc) & is.null(csym) &
      is.null(records) & is.null(dtypes)) {
    uc <- unit_cell(10)
    ruc <- create_rec_unit_cell(uc)
    csym <- cryst_symm("P 2 3")
    records <- generate_miller(uc=uc,SG=csym$SG,reso=5)
    dtypes <- c("H","H","H","S")

    # S3 class
    mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
              dtypes=dtypes),class = "merged_reflections")

    return(mrefs)
  }

  # Input is only 'rec_unit_cell'. Symmetry is chosen as the
  # lowest compatible with the corresponding unit cell. Records
  # are generated up to 5 angstroms resolution
  if (!is.null(ruc) & is.null(csym) &
      is.null(records) & is.null(dtypes)) {
    # Check input is 'rec_unit_cell'
    ans <- check_rec_unit_cell_validity(ruc)
    if (!ans) {
      msg <- paste("Input is not a valid object",
                   "of class 'rec_unit_cell'.\n")
      cat(msg)

      return(NULL)
    }

    # Unit cell needed for the rest
    uc <- create_unit_cell(ruc)
    csym <- lowest_uc_compatible_SG(uc)
    records <- generate_miller(uc,csym$SG,5)
    dtypes <- c("H","H","H","S")

    # S3 class
    mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
                  dtypes=dtypes),class = "merged_reflections")

    return(mrefs)
  }

  # Input is only 'cryst_symm'. The unit cell is chosen as to
  # be compatible with the symmetry. Records are generated up
  # to 5 angstroms resolution
  if (is.null(ruc) & !is.null(csym) &
      is.null(records) & is.null(dtypes)) {
    # Check input is 'cryst_symm'
    ans <- check_cryst_symm_validity(csym)
    if (!ans) {
      msg <- paste("Input is not a valid object",
                   "of class 'cryst_symm'.\n")
      cat(msg)

      return(NULL)
    }

    # Carry on
    ruc <- create_rec_unit_cell(csym)
    uc <- create_unit_cell(ruc)
    records <- generate_miller(uc,csym$SG,5)
    dtypes <- c("H","H","H","S")

    # S3 class
    mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
                  dtypes=dtypes),class = "merged_reflections")

    return(mrefs)
  }

  # Input is 'ruc' and 'cryst_symm'. Check compatibility and
  # generate data
  if (!is.null(ruc) & !is.null(csym) &
      is.null(records) & is.null(dtypes)) {
    # Check both inputs
    ans1 <- check_rec_unit_cell_validity(ruc)
    ans2 <- check_cryst_symm_validity(csym)
    if (!ans1 & ans2) {
      msg <- paste("First input is not a valid object",
                   "of class 'rec_unit_cell'.\n")
      cat(msg)

      return(NULL)
    }
    if (ans1 & !ans2) {
      msg <- paste("Second input is not a valid object",
                   "of class 'cryst_symm'.\n")
      cat(msg)

      return(NULL)
    }
    if (!ans1 & !ans2) {
      msg <- paste("First input is not a valid object",
                   "of class 'rec_unit_cell' and second",
                   "input is not a valid object of",
                   "class 'cryst_symm'.\n")
      cat(msg)

      return(NULL)
    }

    # Carry on. Check unit cell and symmetry compatibility
    uc <- create_unit_cell.rec_unit_cell(ruc)
    ans <- check_validity(uc,csym)
    if (!ans) {
      msg <- paste("'ruc' object is not compatible with",
                   "'csym' object because cell parameters",
                   "are not compatible with symmetry.\n")
      cat(msg)

      return(NULL)
    }

    # All checks passed. Produce data
    records <- generate_miller(uc,csym$SG,5)
    dtypes <- c("H","H","H","S")

    # S3 class
    mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
                  dtypes=dtypes),class = "merged_reflections")

    return(mrefs)
  }

  # None of the parameters is NULL
  if (!is.null(ruc) & !is.null(csym) &
      !is.null(records) & !is.null(dtypes)) {
    # Check validity of ruc
    ans <- check_rec_unit_cell_validity(ruc)
    if (!ans) {
      msg <- paste("'ruc' is not a valid object of class",
                   "'rec_unit_cell'.\n")
      cat(msg)

      return(NULL)
    }
    # Check validity of csym
    ans <- check_cryst_symm_validity(csym)
    if (!ans) {
      msg <- paste("'csym' is not a valid object of class",
                   "'cryst_symm'.\n")
      cat(msg)

      return(NULL)
    }
    # Check 'records' is a data frame
    ans <- is(records,"data.frame")
    if (!ans) {
      msg <- paste("'records' is not a valid object",
                   "of class 'data.frame'.\n")
      cat(msg)

      return(NULL)
    }
    # Check 'records' has at least 3 columns
    ans <- length(records[1,]) >= 3
    if (!ans) {
      msg <- paste("'records' must have at least 3 columns,",
                   "the 3 Miller indices, H, K, L.\n")
      cat(msg)

      return(NULL)
    }
    # Check 'dtypes' is a character vector
    ans <- is(dtypes,"character")
    if (!ans) {
      msg <- paste("'dtypes' is not of a valid object",
                   "of class 'character'.\n")
      cat(msg)

      return(NULL)
    }
    # Allowed dtypes
    Vdtypes <- c("H","S","J","F","D","Q","G","L",
                 "K","M","E","P","W","A","B","I","R")
    for (dt in dtypes) {
      if (!(dt %in% Vdtypes)) {
        msg <- "One or more 'dtypes' are not recognised.\n"
        cat(msg)

        return(NULL)
      }
    }
    # Number of columns = length(dtypes)
    if (ncol(records) != length(dtypes)) {
      ans <- FALSE
      msg <- paste("Number of columns of 'records' must be",
                   "equal to length of 'dtypes'.\n")
      cat(msg)

      return(NULL)
    }
    # First three columns have to be of dtype "H"
    ans <- (dtypes[1] == "H" & dtypes[2] == "H" &
      dtypes[3] == "H")
    if (!ans) {
      msg <- paste("The first three columns of 'records' has",
                   "to be of dtypes 'H'.\n")
      cat(msg)

      return(NULL)
    }
    # Unit cell compatible with symmetry?
    uc <- create_unit_cell(ruc)
    ans <- check_validity(uc,csym)
    if (!ans) {
      msg <- paste("The reciprocal unit cell of this object",
                   "is not compatible with its symmetry.\n")
      cat(msg)

      return(NULL)
    }
    # Check systematic absences
    idx <- sysabs(records[,1:3],csym$SG)
    if (length(idx) != length(records[,1])) {
      msg <- paste("Systematic absences have been",
                   "eliminated from 'records'.\n")
      cat(msg)
      records <- records[idx,]
    }
  }

  # All clear: proceed
  # S3 class
  mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
           dtypes=dtypes),class = "merged_reflections")

  return(mrefs)
}




#' Print method for an object of class "merged_reflections".
#'
#' The output includes details on the unit cell, the crystal
#' symmetry and the first 10 records (data).
#'
#' @param x An object of class "merged_reflections".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create a default 'merged_reflections' object
#' mrefs <- merged_reflections()
#'
#' # Display its value
#' print(mrefs)
#'
#' @rdname print.unit_cell
#' @export
print.merged_reflections <- function(x,...) {
  cat("This is an object of class 'merged_reflections'\n")

  # Unit cell
  uc <- create_unit_cell(x$ruc)
  cat("\n")
  cat("Its unit cell has sides:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    angstroms.\n",
                 uc$a,uc$b,uc$c)
  cat(msg)
  cat("And angles:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    degrees.\n",
                 uc$alpha[1],uc$beta[1],uc$gamma[1])
  cat(msg)

  # Symmetry
  cat("\n")
  msg <- "Its symmetry is represented by space group: %s.\n"
  msg <- sprintf(msg,x$csym$SG)
  cat(msg)

  # Records
  cat("These are the first few data lines:\n\n")
  print(x$records[1:10,])

  invisible(x)
}



#' Default method for generic "create_merged_reflections"
#'
#' This method is an alternative call to 'merged_reflections'.
#'.
#' @param ruc An object of class 'rec_unit_cell'.
#' @param csym An object of class 'cryst_symm'.
#' @param records A data frame containing all reflections coming
#'        from the x-ray data collection on the crystal. This
#'        data frame must include at least the three Miller
#'        indices, H, K, L (of dtype "H").
#' @param dtypes A character vector whose length is the same as
#'        the number of columns in 'records'. One character (a
#'        capital letter) is associated with each type of data.
#'        For example, a Miller index is of dtype "H"; a structure
#'        amplitude is of dtype "F"; an anomalous difference is of
#'        dtype "D"; etc (see details later).
#' @params ... Additional arguments passed to the
#'             create_merged_reflections methods.
#' @return An object of class "merged_reflections". It is a named
#'         list of length 4 whose names are:
#'         \describe{
#'           \item{ruc}{An object of class "rec_unit_cell".}
#'           \item{csym}{An object of class "cryst_symm".}
#'           \item{records}{A data frame containing the data.}
#'           \item{dtypes}{A character vector containing the
#'                 type of data (Miller indices, structure
#'                 factors, etc).}
#'         }
#'
#' @examples
#' # Create merged data for a cubic (10 angstrom) unit cell
#' # of space group P 2 3. Data up to 5 angstroms resolution
#' mrefs <- create_merged_reflections()
#' print(mrefs)
#'
#' @seealso
#' \code{\link{merged_reflections}}
#'
#' @rdname create_merged_reflections.default
#' @export
create_merged_reflections.default <- function(ruc=NULL,csym=NULL,
                      records=NULL,dtypes=NULL,...) {
  # Makes use of all constructions in "merged_reflections"
  mrefs <- merged_reflections(ruc=ruc,csym=csym,records=records,
                              dtypes=dtypes)

  return(mrefs)
}
