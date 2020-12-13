#
# This file is part of the cry package
#

#' Validity and compatibility of cry objects
#'
#' Compatibility of cry objects The objects that can be created in cry are subject to issues of
#' compatibility of the parameters forming them. For example, the unit cell of a cubic system
#' cannot have the a, b, c sides different from each other. The present function returns TRUE if
#' all parts composing the object are compatible with each other. Otherwise it returns FALSE and
#' one or more warnings with details.
#'
#' @param x A first object of one of the new cry classes.
#' @param y A second  object of one of the new cry classes.
#' @return ans A logical value. ans = TRUE means that either the parameters of the single
#'  object (if only one input is provided) are valid parameters, or that the two objects are
#'  compatible with each other.
#'
#' @examples
#' # Create a cubic cell with side 50
#' uc <- create_unit_cell(50)
#'
#' # Create an object of class "cryst_symm"
#' crsym <- cryst_symm("P m -3")
#'
#' # Are they compatible with each other?
#' check_validity(uc,crsym)
#'
#' @export
check_validity <- function(x,y=NULL) {
  # Default value is true
  ans <- TRUE

  # Class of first object
  classX <- class(x)
  classY <- class(y)

  # One or two objects?
  if (classY == "NULL") {
    # What class?
    if (classX == "angle") {
      ans <- check_angle_validity(x)
    } else if (classX == "bravais") {
      ans <- check_bravais_validity(x)
    } else if (classX == "unit_cell") {
      ans <- check_unit_cell_validity(x)
    } else if (classX == "rec_unit_cell") {
      ans <- check_rec_unit_cell_validity(x)
    } else if (classX == "cryst_symm") {
      ans <- check_cryst_symm_validity(x)
    } else {
      warning("No checks are implemented for this input class.")
    }
  } else {
    # Some pairs are checked, others are not
    if ((classX == "bravais" & classY == "unit_cell") |
        (classX == "unit_cell" & classY == "bravais")) {
      #
      ### Unit Cell and Bravais Lattice ###
      #

      # Individual objects OK
      ans1 <- check_bravais_validity(x)
      ans2 <- check_unit_cell_validity(y)
      if (ans1 & ans2) {
        bt <- x
        uc <- y
      } else {
        bt <- y
        uc <- x
      }

      # Sides and angles
      a <- uc$a
      b <- uc$b
      c <- uc$c
      aa <- uc$alpha[1]
      bb <- uc$beta[1]
      cc <- uc$gamma[1]

      # Extract lattice and centering
      latt <- substr(bt$bt,1,1)
      centring <- substr(bt$bt,2,2)

      # Triclinic (no checks)
      #
      # Monoclinic (unique axis can be b or c)
      if (latt == "m") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.00001 | dc > 0.000001) {
          if (da > 0.000001 | db > 0.000001) {
            warning("Angles of 'unit_cell' object are not compatible with 'bravais' object.")
            ans <- FALSE
          }
        }
      }
      # Orthorombic
      if (latt == "o") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
      }
      # Tetragonal
      if (latt == "t") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
        diffab <- abs(a-b)
        if (diffab > 0.000001) {
          warning("Sides of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
      }
      # Cubic
      if (latt == "c") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
        diffab <- abs(a-b)
        diffbc <- abs(b-c)
        if (diffab > 0.000001 | diffbc > 0.000001) {
          warning("Sides of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
      }
      # Hexagonal
      if (latt == "h") {
        # First check angles
        # Hexagonal
        ansHe <- TRUE
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(120-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) ansHe <- FALSE
        # Rombohedral
        ansRo <- TRUE
        diffab <- abs(aa-bb)
        diffbc <- abs(bb-cc)
        da90 <- abs(90-aa)
        db90 <- abs(90-bb)
        dc90 <- abs(90-cc)
        if (diffab > 0.000001 | diffbc > 0.000001 |
            da90 < 0.000001 | db90 < 0.000001 | dc90 < 0.000001) ansRo <- FALSE
        if (!ansHe & !ansRo) {
          warning("Angles of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
        # Next check sides
        # Hexagonal
        ansHe <- TRUE
        diffab <- abs(a-b)
        if (diffab < 0.000001) ansHe <- FALSE
        # Rombohedral
        ansRo <- TRUE
        diffab <- abs(a-b)
        diffbc <- abs(b-c)
        if (diffab > 0.000001 | diffbc > 0.000001) ansRo <- FALSE
        if (!ansHe & !ansRo) {
          warning("Sides of 'unit_cell' object are not compatible with 'bravais' object.")
          ans <- FALSE
        }
      }
    } else if ((classX == "cryst_symm" & classY == "unit_cell") |
               (classX == "unit_cell" & classY == "cryst_symm")) {
      #
      ### Unit Cell and Symmetry ###
      #

      # Individual objects OK
      ans1 <- check_cryst_symm_validity(x)
      ans2 <- check_unit_cell_validity(y)
      if (ans1 & ans2) {
        crsym <- x
        uc <- y
      } else {
        crsym <- y
        uc <- x
      }

      # From cryst_symm to bravais (to crystal family, really!)
      ltmp <- extract_symmetry_info(crsym$SG)
      xxx <- strsplit(ltmp$NUMBER," ")[[1]]
      gn <- as.integer(xxx[length(xxx)])
      tmp <- crystal_family(gn)
      if (tmp == "TRICLINIC") latt <- "a"
      if (tmp == "MONOCLINIC") latt <- "m"
      if (tmp == "ORTHOROMBIC") latt <- "o"
      if (tmp == "TETRAGONAL") latt <- "t"
      if (tmp == "HEXAGONAL") latt <- "h"
      if (tmp == "CUBIC") latt <- "c"

      # Now we can proceed as in the unit_cell - bravais case

      # Sides and angles
      a <- uc$a
      b <- uc$b
      c <- uc$c
      aa <- uc$alpha[1]
      bb <- uc$beta[1]
      cc <- uc$gamma[1]

      # Triclinic (no checks)
      #
      # Monoclinic (unique axis can be b or c)
      if (latt == "m") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.00001 | dc > 0.000001) {
          if (da > 0.000001 | db > 0.000001) {
            warning("Angles of 'unit_cell' object are not compatible with 'cryst_symm' object.")
            ans <- FALSE
          }
        }
      }
      # Orthorombic
      if (latt == "o") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
      }
      # Tetragonal
      if (latt == "t") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
        diffab <- abs(a-b)
        if (diffab > 0.000001) {
          warning("Sides of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
      }
      # Cubic
      if (latt == "c") {
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(90-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) {
          warning("Angles of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
        diffab <- abs(a-b)
        diffbc <- abs(b-c)
        if (diffab > 0.000001 | diffbc > 0.000001) {
          warning("Sides of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
      }
      # Hexagonal
      if (latt == "h") {
        # First check angles
        # Hexagonal
        ansHe <- TRUE
        da <- abs(90-aa)
        db <- abs(90-bb)
        dc <- abs(120-cc)
        if (da > 0.000001 | db > 0.000001 | dc > 0.000001) ansHe <- FALSE
        # Rombohedral
        ansRo <- TRUE
        diffab <- abs(aa-bb)
        diffbc <- abs(bb-cc)
        da90 <- abs(90-aa)
        db90 <- abs(90-bb)
        dc90 <- abs(90-cc)
        if (diffab > 0.000001 | diffbc > 0.000001 |
            da90 < 0.000001 | db90 < 0.000001 | dc90 < 0.000001) ansRo <- FALSE
        if (!ansHe & !ansRo) {
          warning("Angles of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
        # Next check sides
        # Hexagonal
        ansHe <- TRUE
        diffab <- abs(a-b)
        if (diffab < 0.000001) ansHe <- FALSE
        # Rombohedral
        ansRo <- TRUE
        diffab <- abs(a-b)
        diffbc <- abs(b-c)
        if (diffab > 0.000001 | diffbc > 0.000001) ansRo <- FALSE
        if (!ansHe & !ansRo) {
          warning("Sides of 'unit_cell' object are not compatible with 'cryst_symm' object.")
          ans <- FALSE
        }
      }
    } else {
      warning("No checks are implemented for this pair of input classes.")
    }
  }

  return(ans)
}


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


#' Validity and compatibility of a cry object of class 'unit_cell'
#'
#' An object of class 'unit_cell' is a named list of length 6. The first three fields are
#' numeric, the last three of class 'angle'.
#'
#' @param x Object of class 'unit_cell'.
#' @return ans A logical value. TRUE means that the input is a valid object of class'unit_cell'.
#'
#' @examples
#' # Create an object of class 'unit_cell'
#' x <- create_unit_cell()
#'
#' # Check its validity
#' check_unit_cell_validity(x)
#'
#' # Now change two fields
#' x$a <- "a"
#' x$alpha <- 123
#'
#' # Check validity again
#' check_unit_cell_validity(x)
#'
#' @export
check_unit_cell_validity <- function(x) {
  ans <- TRUE
  fans <- is(x,"unit_cell")
  if (!fans) {
    ans <- FALSE
  } else {
    # Check a, b, c are numeric
    if (!is.numeric(x$a) | !is.numeric(x$b) | !is.numeric(x$c)) {
      warning("One or more of the a, b, c of this object of class 'unit_cell' are not numeric.")
      ans <- FALSE
    }
    if (!is(x$alpha,"angle") | !is(x$beta,"angle") | !is(x$gamma,"angle")) {
      msg <- paste0("One or more of the alpha, beta, gamma of this object of class ",
                    "'unit_cell' are not of class 'angle'.")
      warning(msg)
      ans <- FALSE
    }
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'rec_unit_cell'
#'
#' An object of class 'rec_unit_cell' is a named list of length 6. The first three fields are
#' numeric, the last three of class 'angle'.
#'
#' @param x Object of class 'rec_unit_cell'.
#' @return ans A logical value. TRUE means that the input is a valid object of class'rec_unit_cell'.
#'
#' @examples
#' # Create an object of class 'rec_unit_cell'
#' x <- create_rec_unit_cell()
#'
#' # Check its validity
#' check_rec_unit_cell_validity(x)
#'
#' # Now change two fields
#' x$ar <- "ar"
#' x$alphar <- 123
#'
#' # Check validity again
#' check_rec_unit_cell_validity(x)
#'
#' @export
check_rec_unit_cell_validity <- function(x) {
  ans <- TRUE
  fans <- is(x,"rec_unit_cell")
  if (!fans) {
    ans <- FALSE
  } else {
    # Check ar, br, cr are numeric
    if (!is.numeric(x$ar) | !is.numeric(x$br) | !is.numeric(x$cr)) {
      warning("One or more of the ar, br, cr of this object of class 'rec_unit_cell' are not numeric.")
      ans <- FALSE
    }
    if (!is(x$alphar,"angle") | !is(x$betar,"angle") | !is(x$gammar,"angle")) {
      msg <- paste0("One or more of the alphar, betar, gammar of this object of class ",
                    "'rec_unit_cell' are not of class 'angle'.")
      warning(msg)
      ans <- FALSE
    }
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'cryst_symm'
#'
#' An object of class 'cryst_symm' is a named list of length 4. The first field is a
#' character string, the second field is a \eqn{3\times 3} array and the third and fourth
#' field are \eqn{3 \times 1} arrays.
#'
#' @param x Object of class 'cryst_symm'.
#' @return ans A logical value. TRUE means that the input is a valid object of class'cryst_symm'.
#'
#' @examples
#' # Create an object of class 'cryst_symm'
#' x <- cryst_symm(15)
#'
#' # Check its validity
#' check_cryst_symm_validity(x)
#'
#' # Now change two fields
#' x$SG <- "PPP"
#' x$PG[[1]] <- matrix(rep(0,times=9),ncol=3)
#'
#' # Check validity again
#' check_cryst_symm_validity(x)
#'
#' @export
check_cryst_symm_validity <- function(x) {
  ans <- TRUE
  fans <- is(x,"cryst_symm")
  if (!fans) {
    warning("This is not an object of class 'cryst_symm'.")
    ans <- FALSE

    return(ans)
  } else {
    # Names must be SG, PG, T and C
    tmpnames <- names(x)
    if (tmpnames[1] != "SG" | tmpnames[2] != "PG" | tmpnames[3] != "T" | tmpnames[4] != "C") {
      msg <- paste0("Names for the fields of this object of class 'cryst_symm' must be ",
                    "'SG', 'PG', 'T' and 'C'.")
      warning(msg)
      ans <- FALSE

      return(ans)
    }
    # First field must be a string character
    if (!is.character(x$SG)) {
      warning("SG of the object of class 'cryst_symm' must be a string character.")
      ans <- FALSE
    }
    # First field is a valid space group symbol
    ans1 <- extract_symmetry_info(x$SG)
    if (is.null(ans1)) {
      msg <- paste0("SG of the object of class 'cryst_symm' must be a valid space group ",
                    "extended Hermann-Mauguin symbol.")
      warning(msg)
      ans <- FALSE
    }
    # Second field must be a list of 3 X 3 matrices
    ans1 <- is.list(x$PG)
    if (!ans1) {
      warning("'PG' of this object of class 'cryst_symm' must be a list of 3X3 matrices.")
      ans <- FALSE

      return(FALSE)
    }
    ans2 <- TRUE
    for (i in 1:length(x$PG)) {
      if (!is.matrix(x$PG[[i]])) ans2 <- FALSE
    }
    if (!ans2) {
      warning("'PG' of this object of class 'cryst_symm' must be a list of 3X3 matrices.")
      ans <- FALSE

      return(ans)
    }
    # Third field must be a 3 X 1 vector
    ans1 <- is.list(x$T)
    if (!ans1) {
      warning("'T' of this object of class 'cryst_symm' must be a list of 3X1 vectors.")
      ans <- FALSE

      return(FALSE)
    }
    ans2 <- TRUE
    for (i in 1:length(x$T)) {
      if (!is.matrix(x$T[[i]])) ans2 <- FALSE
    }
    if (!ans2) {
      warning("'T' of this object of class 'cryst_symm' must be a list of 3X1 vectors.")
      ans <- FALSE

      return(ans)
    }
    # Fourth field must be a list of 3 X 1 vectors
    ans1 <- is.list(x$C)
    if (!ans1) {
      warning("'C' of this object of class 'cryst_symm' must be a list of 3X1 vectors.")
      ans <- FALSE

      return(FALSE)
    }
    ans2 <- TRUE
    for (i in 1:length(x$C)) {
      if (!is.matrix(x$C[[i]])) ans2 <- FALSE
    }
    if (!ans2) {
      warning("'C' of this object of class 'cryst_symm' must be a list of 3X1 vectors.")
      ans <- FALSE

      return(ans)
    }
  }

  return(ans)
}
