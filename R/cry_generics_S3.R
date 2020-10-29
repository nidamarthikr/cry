# This file is part of the cry package
#

#########
### S3 generics specific for cry objects
########

#' S3 generic to create unit_cell objects
#'
#' The unit_cell object can be created starting from specific objects, files, etc.
#'
#' @param x An object used to select a method.
#' @param ... Further arguments passed to or from other methods.
#' @examples
#' # Create a unit_cell in default (no arguments)
#' uc <- create_unit_cell()
#' print(uc)
#' @export
create_unit_cell <- function(x,...) UseMethod("create_unit_cell")


#' S3 generic to create rec_unit_cell objects
#'
#' The rec_unit_cell object can be created starting from specific objects, files, etc.
#'
#' @param x An object used to select a method.
#' @param ... Further arguments passed to or from other methods.
#' @examples
#' # Create a rec_unit_cell in default (no arguments)
#' ruc <- create_rec_unit_cell()
#' print(ruc)
#' @export
create_rec_unit_cell <- function(x,...) UseMethod("create_rec_unit_cell")
