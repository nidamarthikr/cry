#
# This file is part of the cry package
#

## Validity and compatibility of cry objects
##
## Compatibility of cry objects The objects that can be created in cry are subject to issues of compatibility of the
## parameters forming them. For example, the unit cell of a cubic system cannot have the
## a, b, c sides different from each other. The present function returns TRUE if all parts
## composing the object are compatible with each other.
##
## @param x A first object of one of the new cry classes.
## @param y A second  object of one of the new cry classes.
## @return ans A logical value. ans = TRUE means that either the parameters of the single
##  object (if only one input is provided) are valid parameters, or that the two objects are
##  compatible with each other.
##
## @examples
## # Create a cubic cell with side 50
## uc <- create_unit_cell(50)
##
## # Create an object of class "cryst_symm"
## crsym <- cryst_symm("P m -3")
##
## # Are they compatible with each other?
## check_validity(uc,crsym)
##
## @export
#check_validity <- function(x,y=NULL) {
#  # First check if validity (one input) or compatibility (two inputs) is required
#  if (is.null(y)) {
#    #
#    ## Validity
#    #
#    ans <- TRUE
#
#    # Angles (class "angle")
#
#    }
#
#  } else {
#    #
#    ## Compatibility
#    #
#
#  }
#}


#' Validity and compatibility of a cry object of class 'angle'
#'
#' An object of class 'angle' is a numeric with logical attribure "rad_flag".
#'
#' @param x Object of class angle.
#' @return ans A logical value. TRUE means that the input is a valid object of class'angle'.
#'
#' @examples
#' # Create an object of class angle
#' x <- angle(80)
#'
#' # Check its validity
#' check_angle_validity(x)
#'
#' # Modify the 'rad_flag' attribute
#' attr(x,"rad_flag") <- 12.5
#'
#' # Check its validity
#' check_angle_validity(x)
#'
#' @export
check_angle_validity <- function(x) {
  ans <- TRUE
  fans <- is(x,"angle")
  if (fans) {
    ltmp <- attributes(x)
    ans1 <- "rad_flag" %in% names(ltmp)
    if (ans1) {
      ans2 <- is.logical(attr(x,"rad_flag"))
      if (!ans2) {
        warning("This object of class 'angle' has not a valid 'rad_flag' attribute.\n")
        ans <- FALSE
      }
    } else {
      warning("This object of class 'angle' has not 'rad_flag' attribute.\n")
      ans <- FALSE
    }
    if (!is.numeric(x)) {
      warning("This object of class 'angle' has no numeric value.\n")
      ans <- FALSE
    }

    return(ans)
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'bravais'
#'
#' An object of class 'bravais' is a named list of length 4. The first slot, "bt", is
#' the universally-used two-letter symbol. The second, third and fourth slots are, respectively,
#' "cr_fam" (the corresponding crystal family), "cr_sys" (the corresponding crystal system) and
#' "lt_sys" (the corresponding lattice system).
#'
#' @param x Object of class 'bravais'.
#' @return ans A logical value. TRUE means that the input is a valid object of class'bravais'.
#'
#' @examples
#' # Create an object of class 'bravais'
#' x <- bravais("mP")
#'
#' # Check its validity
#' check_bravais_validity(x)
#'
#' @export
check_bravais_validity <- function(x) {
  ans <- TRUE
  fans <- is(x,"bravais")
  if (!fans) {
    ans <- FALSE
  } else {
    ltmp <- names(x)
    ans1 <- "bt" %in% ltmp
    if (ans1) {
      # More checks
      if (x$bt != "aP" & x$bt != "mP" & x$bt != "mS" & x$bt != "oP" & x$bt != "oS" &
          x$bt != "oF" & x$bt != "oI" & x$bt != "tP" & x$bt != "tI" & x$bt != "hP" &
          x$bt != "hR" & x$bt != "cP" & x$bt != "cF" & x$bt != "cI") {
        warning("The Bravais type 'bt' for this object of class 'bravais' cannot be recognised.")
        ans <- FALSE
      }
    } else {
      warning("'bt' field not included in this object of class 'bravais'.")
      ans <- FALSE
    }
    ans1 <- "cr_fam" %in% ltmp
    if (ans1) {
      # More checks
      if (x$cr_fam != "triclinic" & x$cr_fam != "monoclinic" & x$cr_fam != "orthorombic" &
          x$cr_fam != "tetragonal" & x$cr_fam != "cubic" & x$cr_fam != "cubic" &
          x$cr_fam != "hexagonal") {
        warning("The crystal family for this object of class 'bravais' cannot be recognised.")
        ans <- FALSE
      }
    } else {
      warning("'cr_fam' field not included in this object of class 'bravais'.")
      ans <- FALSE
    }
    ans1 <- "cr_sys" %in% ltmp
    if (ans1) {
      # More checks
      if (x$cr_sys != "triclinic" & x$cr_sys != "monoclinic" & x$cr_sys != "orthorombic" &
          x$cr_sys != "tetragonal" & x$cr_sys != "cubic" & x$cr_sys != "trigonal" &
          x$cr_sys != "hexagonal") {
        warning("The crystal system for this object of class 'bravais' cannot be recognised.")
      }
    } else {
      warning("'cr_sys' field not included in this object of class 'bravais'.")
      ans <- FALSE
    }
    ans1 <- "lt_sys" %in% ltmp
    if (ans1) {
      # More checks
      if (x$lt_sys != "triclinic" & x$lt_sys != "monoclinic" & x$lt_sys != "orthorombic" &
          x$lt_sys != "tetragonal" & x$lt_sys != "cubic" & x$lt_sys != "hexagonal" &
          x$lt_sys != "rhombohedral or hexagonal (centred)") {
        warning("The lattice system for this object of class 'bravais' cannot be recognised.")
        ans <- FALSE
      }
    } else {
      warning("'lt_sys' field not included in this object of class 'bravais'.")
      ans <- FALSE
    }
  }

  # Now check compatibilities

  # Crystal Family
  if (x$bt == "aP" & x$cr_fam != "triclinic") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "mP" | x$bt == "mS") & x$cr_fam != "monoclinic") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" | x$bt == "oI" ) & x$cr_fam != "orthorombic") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "tP" | x$bt == "tI") & x$cr_fam != "tetragonal") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") & x$cr_fam != "cubic") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "hP" | x$bt == "hR") & x$cr_fam != "hexagonal") {
    warning("Bravais type and crystal family are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }

  # Crystal System
  if (x$bt == "aP" & x$cr_sys != "triclinic") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "mP" | x$bt == "mS") & x$cr_sys != "monoclinic") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" | x$bt == "oI" ) & x$cr_sys != "orthorombic") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "tP" | x$bt == "tI") & x$cr_sys != "tetragonal") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") & x$cr_sys != "cubic") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "hP" | x$bt == "hR") & (x$cr_sys != "hexagonal" & x$cr_sys != "trigonal")) {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }

  # Lattice System
  if (x$bt == "aP" & x$cr_sys != "triclinic") {
    warning("Bravais type and crystal system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "mP" | x$bt == "mS") & x$lt_sys != "monoclinic") {
    warning("Bravais type and lattice system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" | x$bt == "oI" ) & x$lt_sys != "orthorombic") {
    warning("Bravais type and lattice system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "tP" | x$bt == "tI") & x$lt_sys != "tetragonal") {
    warning("Bravais type and lattice system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") & x$lt_sys != "cubic") {
    warning("Bravais type and lattice system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }
  if ((x$bt == "hP" | x$bt == "hR") & (x$lt_sys != "hexagonal" &
                                       x$lt_sys != "rhombohedral or hexagonal (centred)")) {
    warning("Bravais type and lattice system are incompatile in this object of class 'bravais'.")
    ans <- FALSE
  }

  # Compatibilities among cr_fam, cr_sys and lt_sys
  log1 <- x$cr_fam == "triclinic" & x$cr_sys == "triclinic" & x$lt_sys == "triclinic"
  log2 <- x$cr_fam == "monoclinic" & x$cr_sys == "monoclinic" & x$lt_sys == "monoclinic"
  log3 <- x$cr_fam == "orthorombic" & x$cr_sys == "orthorombic" & x$lt_sys == "orthorombic"
  log4 <- x$cr_fam == "tetragonal" & x$cr_sys == "tetragonal" & x$lt_sys == "tetragonal"
  log5 <- x$cr_fam == "cubic" & x$cr_sys == "cubic" & x$lt_sys == "cubic"
  log6 <- x$cr_fam == "hexagonal" & x$cr_sys == "trigonal" &
    x$lt_sys == "rhombohedral or hexagonal (centred)"
  log7 <- x$cr_fam == "hexagonal" & x$cr_sys == "hexagonal" & x$lt_sys == "hexagonal"
  if (!log1 & !log2 & !log3 & !log4 & !log5 & !log6 & !log7) {
    msg <- paste0("The crystal family, crystal system and lattice system of this object of class",
                  " 'bravais' are not compatible.")
    warning(msg)
    ans <- FALSE
  }

  return(ans)
}
