#
# This file is part of the cry package
#



#' Constructor for an S3 object of class "merged_reflections".
#'
#' This represents scaled and merged x-ray data from one crystal.
#'
#' If the constructor is used without arguments, the default
#' object created will be composed of 3 reflections, (1,0,0),
#' (0,1,0), (0,0,1) from a cubic crystal with cell of side 10
#' angstroms, and symmetry P 2 3. The only available columns
#' will be of dtype "H", and named H, K, L (the Miller indices).
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
#' @param set An integer number corresponding to the specific
#'            setting for the given space group. Default is 1.
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
#' ruc <- create_rec_unt_cell(uc)
#'
#' # Create symmetry (P n c 2)
#' csym <- cryst_symm(30)
#'
#' # Create a few records
#' records <- expand.grid(H=-2:2,K=-2:2,L=-2:2)
#'
#' # Create merged_reflections object with H, K, L
#' mrefs <- merged_reflections(ruc,csym,records,set=1)
#'
#' @export
merged_reflections <- function(ruc=NULL,csym=NULL,
                               records=NULL,dtypes=NULL,set=1) {
  # If all input parameters are NULL use default values
  if (is.null(ruc) & is.null(csym) &
      is.null(records) & is.null(dtypes)) {
    uc <- unit_cell(10)
    ruc <- create_rec_unit_cell(uc)
    csym <- cryst_symm(195,set=set)
    records <- expand.grid(H=-2:2,K=-2:2,L=-2:2)
    dtypes <- c("H","H","H")
  }

    # S3 class
    mrefs <- structure(list(ruc=ruc,csym=csym,records=records,
                       dtypes=dtypes),class = "merged_reflections")

    return(mrefs)
}
