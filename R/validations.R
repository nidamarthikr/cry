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
#' The available checks for individual objects are:
#' \itemize{
#'   \item \strong{angle} (unit cell angle)
#'   \item \strong{bravais} (Bravais system)
#'   \item \strong{unit_cell} (unit cell)
#'   \item \strong{rec_unit_cell} (reciprocal unit cell)
#'   \item \strong{cryst_symm} (crystallographic symmetry)
#'   \item \strong{merged_reflections} (scaled and merged data)
#' }
#' The available checks for couple of objects are:
#' \itemize{
#'   \item \strong{bravais} and \strong{unit_cell}
#'   \item \strong{unit_cell} and \strong{cryst_symm}
#' }
#'
#' @param x A first object of one of the new cry classes.
#' @param y A second  object of one of the new cry classes.
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
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
#' check_validity(uc,crsym,TRUE)
#'
#' @export
check_validity <- function(x,y=NULL,message=FALSE) {
  # Flag to indicate covered cases
  iscove <- FALSE

  # Default value is true
  ans <- TRUE

  # Class of first and second object
  classX <- class(x)
  classY <- class(y)

  # One or two objects?
  if (classY == "NULL") {
    # What class?
    if (classX == "angle") {
      iscove <- TRUE
      ans <- check_angle_validity(x,message)
    } else if (classX == "bravais") {
      iscove <- TRUE
      ans <- check_bravais_validity(x,message)
    } else if (classX == "unit_cell") {
      iscove <- TRUE
      ans <- check_unit_cell_validity(x,message)
    } else if (classX == "rec_unit_cell") {
      iscove <- TRUE
      ans <- check_rec_unit_cell_validity(x,message)
    } else if (classX == "cryst_symm") {
      iscove <- TRUE
      ans <- check_cryst_symm_validity(x,message)
    } else if (classX =="merged_reflections") {
      iscove <- TRUE
      ans <- check_merged_reflections_validity(x,message)
    } else {
      msg <- paste("No checks are implemented for this",
                  "input class.\n")
      if (message) cat(msg)

      return(FALSE)
    }
  }

  # Next, some pairs are checked, others are not.
  if ((classX == "bravais" & classY == "unit_cell") |
      (classX == "unit_cell" & classY == "bravais")) {

    ### Bravais Lattice and Unit Cell ###
    iscove <- TRUE

    # bravais is x and unit_cell is y
    if (is(y,"bravais")) {
      tmp <- x
      x <- y
      y <- tmp
    }
    bt <- x
    uc <- y

    # Check individual objects
    ans <- check_bravais_validity(bt)
    if (!ans) {
      msg <- paste("Input does not include a valid object",
                   "of class 'bravais'.\n")
      if (message) cat(msg)
      ans <- FALSE

      return(ans)
    }
    ans <- check_unit_cell_validity(uc)
    if (!ans) {
      msg <- paste("Input does not include a valid object",
                   "of class 'unit_cell'.\n")
      if (message) cat(msg)
      ans <- FALSE

      return(ans)
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
    # Monoclinic (unique axis can be a or b or c)
    if (latt == "m") {
      da <- abs(90-aa)
      db <- abs(90-bb)
      dc <- abs(90-cc)
      ntmp <- 3
      if (da > 0.000001) ntmp <- ntmp-1
      if (db > 0.000001) ntmp <- ntmp-1
      if (dc > 0.000001) ntmp <- ntmp-1
      if (ntmp < 2) {
        msg <- paste("Angles of 'unit_cell' object are not",
                     "compatible with 'bravais'","object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
    }

    # Orthorombic
    if (latt == "o") {
      da <- abs(90-aa)
      db <- abs(90-bb)
      dc <- abs(90-cc)
      ntmp <- 3
      if (da > 0.000001) ntmp <- ntmp-1
      if (db > 0.000001) ntmp <- ntmp-1
      if (dc > 0.000001) ntmp <- ntmp-1
      if (ntmp < 3) {
        msg <- paste("Angles of 'unit_cell' object are not",
                     "compatible with 'bravais'","object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
    }

    # Tetragonal
    if (latt == "t") {
      da <- abs(90-aa)
      db <- abs(90-bb)
      dc <- abs(90-cc)
      ntmp <- 3
      if (da > 0.000001) ntmp <- ntmp-1
      if (db > 0.000001) ntmp <- ntmp-1
      if (dc > 0.000001) ntmp <- ntmp-1
      if (ntmp < 3) {
        msg <- paste("Angles of 'unit_cell' object are not",
                     "compatible with 'bravais'","object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
      diffab <- abs(a-b)
      if (diffab > 0.000001) {
        msg <- paste("Sides of 'unit_cell' object are not",
                     "compatible with 'bravais'","object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
    }

    # Cubic
    if (latt == "c") {
      da <- abs(90-aa)
      db <- abs(90-bb)
      dc <- abs(90-cc)
      ntmp <- 3
      if (da > 0.000001) ntmp <- ntmp-1
      if (db > 0.000001) ntmp <- ntmp-1
      if (dc > 0.000001) ntmp <- ntmp-1
      if (ntmp < 3) {
        msg <- paste("Angles of 'unit_cell' object are not",
                     "compatible with 'bravais'","object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
      diffab <- abs(a-b)
      diffbc <- abs(b-c)
      if (diffab > 0.000001 | diffbc > 0.000001) {
        msg <- paste("Angles of 'unit_cell' object are",
                     "not compatible with 'bravais'",
                     "object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
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
      if (da > 0.000001 | db > 0.000001 | dc > 0.000001)
        ansHe <- FALSE
      # Rombohedral
      ansRo <- TRUE
      diffab <- abs(aa-bb)
      diffbc <- abs(bb-cc)
      da90 <- abs(90-aa)
      db90 <- abs(90-bb)
      dc90 <- abs(90-cc)
      if (diffab > 0.000001 | diffbc > 0.000001 |
          da90 < 0.000001 | db90 < 0.000001 | dc90 < 0.000001)
        ansRo <- FALSE
      if (!ansHe & !ansRo) {
        msg <- paste("Angles of 'unit_cell' object are",
                     "not compatible with 'bravais'",
                     "object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
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
      if (diffab > 0.000001 | diffbc > 0.000001)
        ansRo <- FALSE
      if (!ansHe & !ansRo) {
        msg <- paste("Sides of 'unit_cell' object are",
                     "not compatible with 'bravais'",
                     "object.\n")
        if (message) cat(msg)
        ans <- FALSE

        return(ans)
      }
    }
  }

  if ((classX == "cryst_symm" & classY == "unit_cell") |
               (classX == "unit_cell" & classY == "cryst_symm")) {

    ### Unit Cell and Symmetry ###
    iscove <- TRUE

    # cryst_symm is x and unit_cell is y
    if (is(y,"cryst_symm")) {
      tmp <- x
      x <- y
      y <- tmp
    }
    csym <- x
    uc <- y

    # Individual objects OK
    ans1 <- check_cryst_symm_validity(csym)
    ans2 <- check_unit_cell_validity(uc)
    if (!ans1) {
      msg <- paste("Input does not include  a valid object",
                   "of class 'cryst_symm'.\n")
      if (message) cat(msg)
    }
    if (!ans2) {
      msg <- paste("Input does not include a valid object",
                   "of class 'unit_cell'.\n")
      if (message) cat(msg)
    }
    if (!ans1 | !ans2) {
      return(FALSE)
    }

    # Proceed with comparison

    # Derive constrains from symmetry
    vcons <- symm_to_cell_const(csym$SG)

    # Sides and angles
    a <- uc$a
    b <- uc$b
    c <- uc$c
    aa <- uc$alpha[1]
    bb <- uc$beta[1]
    cc <- uc$gamma[1]

    ans <- .ucCcsym(a,b,c,aa,bb,cc,vcons)
    if (!ans) {
      msg <- paste("Unit cell parameters not compatible with",
                   "given symmetry.\n")
      if (message) cat(msg)

      return(ans)
    }
  }

  # Pair not covered by present cases
  if (!iscove) {
    msg <- paste("No checks are implemented for this pair",
                 "of input classes.\n")
    ans <- FALSE
    if (message) cat(msg)
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'angle'
#'
#' An object of class 'angle' is a numeric with logical attribure "rad_flag".
#'
#' @param x Object of class angle.
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
#' @return ans A logical value. TRUE means that the input is
#'             a valid object of class 'angle'.
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
#' check_angle_validity(x,TRUE)
#'
#' @export
check_angle_validity <- function(x,message=FALSE) {
  ans <- is(x,"angle")
  if (!ans) {
    msg <- paste("This is not a valid object of class",
                 "'angle'.\n")
    if (message) cat(msg)

    return(ans)
  }
  # Check attributes
  ltmp <- attributes(x)
  ans <- "rad_flag" %in% names(ltmp)
  if (!ans) {
    msg <- "This object has no 'rad_flag' attribute.\n"
    if (message) cat(msg)

    return(ans)
  }
  ans <- is.logical(attr(x,"rad_flag"))
  if (!ans) {
    msg <- "'rad_flag' is not a logical attribute.\n"
    if (message) cat(msg)

    return(ans)
  }

  # The value of angle has to be numeric
  ans <- is.numeric(x[1])
  if (!ans) {
    msg <- "This object has no numeric value.\n"
    if (message) cat(msg)

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
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
#' @return ans A logical value. TRUE means that the input is
#'             a valid object of class 'bravais'.
#'
#' @examples
#' # Create an object of class 'bravais'
#' x <- bravais("mP")
#'
#' # Check its validity
#' check_bravais_validity(x,TRUE)
#'
#' @export
check_bravais_validity <- function(x,message=FALSE) {
  ans <- is(x,"bravais")
  if (!ans) {
    msg <- paste("This is not a valid object of class",
                 "'bravais'.\n")
    if (message) cat(msg)

    return(ans)
  }
  # Check on Bravais type
  ltmp <- names(x)
  ans <- "bt" %in% ltmp
  if (!ans) {
    msg <- "This object has no 'bt' attribute.\n"
    if (message) cat(msg)

    return(ans)
  }
  ## More checks on 'bt'
  ans2 <- TRUE
  if (x$bt != "aP" & x$bt != "mP" & x$bt != "mS" &
      x$bt != "oP" & x$bt != "oS" & x$bt != "oF" &
      x$bt != "oI" & x$bt != "tP" & x$bt != "tI" &
      x$bt != "hP" & x$bt != "hR" & x$bt != "cP" &
      x$bt != "cF" & x$bt != "cI") {
    ans <- FALSE
    msg <- paste("The Bravais type in 'bt'",
                 "cannot be recognised ")
    if (message) cat(msg)

    return(ans)
  }
  # Check on crystal family
  ans <- "cr_fam" %in% ltmp
  if (!ans) {
    msg <- "This object has no 'cr_fam' attribute.\n"
    if (message) cat(msg)

    return(ans)
  }
  ## More checks on 'cr_fam'
  if (x$cr_fam != "triclinic" & x$cr_fam != "monoclinic" &
      x$cr_fam != "orthorombic" & x$cr_fam != "tetragonal" &
      x$cr_fam != "cubic" & x$cr_fam != "cubic" &
      x$cr_fam != "hexagonal") {
    ans <- FALSE
    msg <- paste("The crystal family in 'cr_fam'",
                 "cannot be recognised.\n")

    return(ans)
  }
  # Check on crystal system
  ans <- "cr_sys" %in% ltmp
  if (!ans) {
    msg <- "This object has no 'cr_sys' attribute.\n"
    if (message) cat(msg)

    return(ans)
  }
  ## More checks on 'cr_sys'
  if (x$cr_sys != "triclinic" & x$cr_sys != "monoclinic" &
      x$cr_sys != "orthorombic" & x$cr_sys != "tetragonal" &
      x$cr_sys != "cubic" & x$cr_sys != "trigonal" &
      x$cr_sys != "hexagonal") {
    ans <- FALSE
    msg <- paste("The crystal system in 'cr_fam'",
                 "cannot be recognised.\n")

    return(ans)
  }
  # Check on lattice system
  ans <- "lt_sys" %in% ltmp
  if (!ans) {
    msg <- "This object has no 'lt_sys' attribute.\n"
    if (message) cat(msg)

    return(ans)
  }
  ## More checks on 'lt_sys'
  if (x$lt_sys != "triclinic" & x$lt_sys != "monoclinic" &
      x$lt_sys != "orthorombic" & x$lt_sys != "tetragonal" &
      x$lt_sys != "cubic" & x$lt_sys != "hexagonal" &
      x$lt_sys != "rhombohedral or hexagonal (centred)") {
    ans <- FALSE
    msg <- paste("The lattice system in 'cr_fam'",
                 "cannot be recognised.\n")

    return(ans)
  }

  # Now check compatibilities

  # Crystal Family
  if (x$bt == "aP" & x$cr_fam != "triclinic") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "mP" | x$bt == "mS") &
       x$cr_fam != "monoclinic") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" |
       x$bt == "oI" ) & x$cr_fam != "orthorombic") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "tP" | x$bt == "tI") &
       x$cr_fam != "tetragonal") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") &
      x$cr_fam != "cubic") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "hP" | x$bt == "hR") &
       x$cr_fam != "hexagonal") {
    msg <- paste("Bravais type and crystal family are",
                 "incompatile in this object of",
                 "class 'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }

  # Crystal System
  if (x$bt == "aP" & x$cr_sys != "triclinic") {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "mP" | x$bt == "mS") &
       x$cr_sys != "monoclinic") {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" |
       x$bt == "oI" ) & x$cr_sys != "orthorombic") {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "tP" | x$bt == "tI") &
       x$cr_sys != "tetragonal") {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") &
       x$cr_sys != "cubic") {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  if ((x$bt == "hP" | x$bt == "hR") &
      (x$cr_sys != "hexagonal" & x$cr_sys != "trigonal")) {
    msg <- paste("Bravais type and crystal system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }

  # Lattice System
  if (x$bt == "aP" & x$cr_sys != "triclinic") {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }
  if ((x$bt == "mP" | x$bt == "mS") &
       x$lt_sys != "monoclinic") {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }
  if ((x$bt == "OP" | x$bt == "oS" | x$bt == "oF" |
       x$bt == "oI" ) & x$lt_sys != "orthorombic") {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }
  if ((x$bt == "tP" | x$bt == "tI") &
       x$lt_sys != "tetragonal") {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }
  if ((x$bt == "cP" | x$bt == "cI" | x$bt == "cF") &
       x$lt_sys != "cubic") {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }
  if ((x$bt == "hP" | x$bt == "hR") &
      (x$lt_sys != "hexagonal" &
       x$lt_sys != "rhombohedral or hexagonal (centred)")) {
    msg <- paste("Bravais type and lattice system are",
                 "incompatile in this object of class",
                 "'bravais'.\n")
    ans <- FALSE

    return(ans)
  }

  # Compatibilities among cr_fam, cr_sys and lt_sys
  log1 <- x$cr_fam == "triclinic" & x$cr_sys == "triclinic" &
          x$lt_sys == "triclinic"
  log2 <- x$cr_fam == "monoclinic" & x$cr_sys == "monoclinic" &
          x$lt_sys == "monoclinic"
  log3 <- x$cr_fam == "orthorombic" &
          x$cr_sys == "orthorombic" & x$lt_sys == "orthorombic"
  log4 <- x$cr_fam == "tetragonal" & x$cr_sys == "tetragonal" &
          x$lt_sys == "tetragonal"
  log5 <- x$cr_fam == "cubic" & x$cr_sys == "cubic" &
          x$lt_sys == "cubic"
  log6 <- x$cr_fam == "hexagonal" & x$cr_sys == "trigonal" &
          x$lt_sys == "rhombohedral or hexagonal (centred)"
  log7 <- x$cr_fam == "hexagonal" & x$cr_sys == "hexagonal" &
          x$lt_sys == "hexagonal"
  if (!log1 & !log2 & !log3 & !log4 & !log5 & !log6 & !log7) {
    msg <- paste("The crystal family, crystal system and",
                 "lattice system of this object of class",
                "'bravais' are not compatible.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'unit_cell'
#'
#' An object of class 'unit_cell' is a named list of length 6. The first three fields are
#' numeric, the last three of class 'angle'.
#'
#' @param x Object of class 'unit_cell'.
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
#' @return ans A logical value. TRUE means that the input is a
#'             valid object of class'unit_cell'.
#'
#' @examples
#' # Create an object of class 'unit_cell'
#' x <- create_unit_cell()
#'
#' # Check its validity
#' check_unit_cell_validity(x)
#'
#' # Now change a field
#' x$alpha <- 123
#'
#' # Check validity again
#' check_unit_cell_validity(x,TRUE)
#'
#' @export
check_unit_cell_validity <- function(x,message=FALSE) {
  ans <- is(x,"unit_cell")
  if (!ans) {
    msg <- paste("This is not a valid object of class",
                 "'unit_cell'.\n")
    if (message) cat(msg)

    return(ans)
  }
  # Check a, b, c are numeric
  if (!is.numeric(x$a) | !is.numeric(x$b) |
      !is.numeric(x$c)) {
    msg <- paste("One or more of the a, b, c of this object",
                 "of class 'unit_cell' are not numeric.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # Check alpha, beta, gamma are valid angles
  ans1 <- check_angle_validity(x$alpha)
  ans2 <- check_angle_validity(x$beta)
  ans3 <- check_angle_validity(x$gamma)
  if (!ans1 | !ans2 | !ans3) {
    msg <- paste("One or more of the alpha, beta, gamma",
                 "of this object of class 'unit_cell'",
                 "are not of class 'angle'.\n")
    if (message) cat(msg)
    ans <- FALSE

    return(ans)
  }

  # Check sides
  a <- x$a
  b <- x$b
  c <- x$c
  if (a <= 0 | b <= 0 | c <= 0) {
    ans <- FALSE
    msg <- "A unit cell must have positive sides.\n"
    if (message) cat(msg)

    return(ans)
  }

  # Check angle values
  aa <- x$alpha[1]
  bb <- x$beta[1]
  cc <- x$gamma[1]
  if (aa < 0 | bb < 0 | cc < 0 |
      aa > 180 | bb > 180 | cc > 180) {
    msg <- paste("A unit cell must have angles between",
                 "0 and 180 degrees.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aa+bb+cc
  if (ss >= 360) {
    msg <- "A unit cell with these angles cannot exist.\n"
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aa+bb-cc
  if (ss <= 0 | ss >= 360) {
    msg <- "A unit cell with these angles cannot exist.\n"
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aa-bb+cc
  if (ss <= 0 | ss >= 360) {
    msg <- "A unit cell with these angles cannot exist.\n"
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- -aa+bb+cc
  if (ss <= 0 | ss >= 360) {
    msg <- "A unit cell with these angles cannot exist.\n"
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }

  return(ans)
}


#' Validity and compatibility of a cry object of class 'rec_unit_cell'
#'
#' An object of class 'rec_unit_cell' is a named list of length 6. The first three fields are
#' numeric, the last three of class 'angle'.
#'
#' @param x Object of class 'rec_unit_cell'.
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
#' @return ans A logical value. TRUE means that the input is a
#'             valid object of class'rec_unit_cell'.
#'
#' @examples
#' # Create an object of class 'rec_unit_cell'
#' x <- create_rec_unit_cell()
#'
#' # Check its validity
#' check_rec_unit_cell_validity(x)
#'
#' # Now change a field
#' x$alphar <- 123
#'
#' # Check validity again
#' check_rec_unit_cell_validity(x,TRUE)
#'
#' @export
check_rec_unit_cell_validity <- function(x,message=FALSE) {
  ans <- is(x,"rec_unit_cell")
  if (!ans) {
    msg <- paste("This is not a valid object of class",
                 "'rec_unit_cell'.\n")
    if (message) cat(msg)

    return(ans)
  }
  # Check ar, br, cr are numeric
  if (!is.numeric(x$ar) | !is.numeric(x$br) |
      !is.numeric(x$cr)) {
    msg <- paste("One or more of the ar, br, cr of this object",
                 "of class 'rec_unit_cell' are not numeric.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # Check alphar, betar, gammar are valid angles
  ans1 <- check_angle_validity(x$alphar)
  ans2 <- check_angle_validity(x$betar)
  ans3 <- check_angle_validity(x$gammar)
  if (!ans1 | !ans2 | !ans3) {
    msg <- paste("One or more of the alphar, betar, gammar",
                 "of this object of class 'rec_unit_cell'",
                 "are not of class 'angle'.\n")
    if (message) cat(msg)
    ans <- FALSE

    return(ans)
  }

  # Check sides
  ar <- x$ar
  br <- x$br
  cr <- x$cr
  if (ar <= 0 | br <= 0 | cr <= 0) {
    ans <- FALSE
    msg <- "A reciprocal unit cell must have positive sides.\n"
    if (message) cat(msg)

    return(ans)
  }

  # Check angle values
  aar <- x$alphar[1]
  bbr <- x$betar[1]
  ccr <- x$gammar[1]
  if (aar < 0 | bbr < 0 | ccr < 0 |
      aar > 180 | bbr > 180 | ccr > 180) {
    msg <- paste("A reciprocal unit cell must have angles",
                 "between 0 and 180 degrees.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aar+bbr+ccr
  if (ss >= 360) {
    msg <- paste("A reciprocal unit cell with these angles",
                 "cannot exist.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aar+bbr-ccr
  if (ss <= 0 | ss >= 360) {
    msg <- paste("A reciprocal unit cell with these angles",
                 "cannot exist.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- aar-bbr+ccr
  if (ss <= 0 | ss >= 360) {
    msg <- paste("A reciprocal unit cell with these angles",
                 "cannot exist.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ss <- -aar+bbr+ccr
  if (ss <= 0 | ss >= 360) {
    msg <- paste("A reciprocal unit cell with these angles",
                 "cannot exist.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
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
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any
#'                (default is FALSE, i.e. no message printed).
#' @return ans A logical value. TRUE means that the input is a valid object of class'cryst_symm'.
#'
#' @examples
#' # Create an object of class 'cryst_symm'
#' x <- cryst_symm(15)
#'
#' # Check its validity
#' check_cryst_symm_validity(x)
#'
#' # Now change a field
#' x$PG[[1]] <- matrix(rep(0,times=9),ncol=3)
#'
#' # Check validity again
#' check_cryst_symm_validity(x,TRUE)
#'
#' @export
check_cryst_symm_validity <- function(x,message=FALSE) {
  ans <- TRUE
  fans <- is(x,"cryst_symm")
  if (!fans) {
    msg <- paste("This is not a valid object of class",
                  "'cryst_symm'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # Names must be SG, PG, T and C
  tmpnames <- names(x)
  if (tmpnames[1] != "SG" | tmpnames[2] != "PG" |
      tmpnames[3] != "T" | tmpnames[4] != "C") {
    msg <- paste("Names for the fields of this object of",
                 "class 'cryst_symm' must be 'SG', 'PG',",
                  "'T' and 'C'.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # SG must be a string character
  if (!is.character(x$SG)) {
    msg <- paste("SG of the object of class 'cryst_symm'",
                  "must be a string character.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ## It must be a valid space group symbol
  ans1 <- extract_symmetry_info(x$SG)
  if (is.null(ans1)) {
    msg <- paste("SG of the object of class 'cryst_symm'",
                 "must be a valid space group",
                 "extended Hermann-Mauguin symbol.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # PG must be a list of 3 X 3 matrices
  ans1 <- is.list(x$PG)
  if (!ans1) {
    msg <- paste("'PG' of this object of class 'cryst_symm'",
                 "must be a list of 3X3 matrices.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:length(x$PG)) {
    if (!is.matrix(x$PG[[i]])) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("'PG' of this object of class 'cryst_symm'",
                 "must be a list of 3X3 matrices.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ## They must be the correct 3X3 matrices
  ltmp <- syminfo_to_matrix_list(x$SG)
  ntrue <- length(ltmp$PG)
  ncur <- length(x$PG)
  ans2 <- ntrue == ncur
  if (!ans2) {
    msg <- paste("'PG' of this object of class 'cryst_symm'",
                 "has an incorrect number of 3X3 matrices.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:ntrue) {
    tans <- all.equal(ltmp$PG[[i]],x$PG[[i]])
    if (tans != TRUE) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("One or more matrices of the 'PG' part",
                 "tof his object of class 'cryst_symm'",
                 "are incorrect.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # T must be a list of 3 X 1 vectors
  ans1 <- is.list(x$T)
  if (!ans1) {
    msg <- paste("'T' of this object of class 'cryst_symm'",
                 "must be a list of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:length(x$T)) {
    if (!is.matrix(x$T[[i]])) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("'T' of this object of class 'cryst_symm'",
                 "must be a list of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ## They must be the correct 3X1 vectors
  ntrue <- length(ltmp$T)
  ncur <- length(x$T)
  ans2 <- ntrue == ncur
  if (!ans2) {
    msg <- paste("'T' of this object of class 'cryst_symm'",
                 "has an incorrect number of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:ntrue) {
    tans <- all.equal(ltmp$T[[i]],x$T[[i]])
    if (tans != TRUE) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("One or more vectors of the'T' part",
                 "of this object of class 'cryst_symm'",
                 "are incorrect.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  # C must be a list of 3 X 1 vectors
  ans1 <- is.list(x$C)
  if (!ans1) {
    msg <- paste("'C' of this object of class 'cryst_symm'",
                 "must be a list of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:length(x$C)) {
    if (!is.matrix(x$C[[i]])) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("'C' of this object of class 'cryst_symm'",
                 "must be a list of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ## They must be the correct 3X1 vectors
  ntrue <- length(ltmp$C)
  ncur <- length(x$C)
  ans2 <- ntrue == ncur
  if (!ans2) {
    msg <- paste("'C' of this object of class 'cryst_symm'",
                 "has an incorrect number of 3X1 vectors.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }
  ans2 <- TRUE
  for (i in 1:ntrue) {
    tans <- all.equal(ltmp$C[[i]],x$C[[i]])
    if (tans != TRUE) ans2 <- FALSE
  }
  if (!ans2) {
    msg <- paste("One or more vectors of the'C' part",
                 "of this object of class 'cryst_symm'",
                 "are incorrect.\n")
    ans <- FALSE
    if (message) cat(msg)

    return(ans)
  }

  return(ans)
}

#' Validity and compatibility of a cry object of
#' class 'merged_reflections'
#'
#' An object of class 'merged_reflections' is a named list of
#' length 4:
#' \describe{
#'           \item{ruc}{An object of class "rec_unit_cell".}
#'           \item{csym}{An object of class "cryst_symm".}
#'           \item{records}{A data frame containing the data.}
#'           \item{dtypes}{A character vector containing the
#'                 type of data (Miller indices, structure
#'                 factors, etc).}
#'         }
#' Internal consistency must be displayed between the object
#' 'ruc' and the object 'csym' because groups of crystallographic
#' symmetries are compatible only with certain unit cells (and,
#' accordingly, certain reciprocal cells).
#' It is not possible to check consistency between dtypes and
#' the nature of data in each column of the data frame 'records',
#' but a check about length of 'dtypes' and number of columns is
#' possible. Therefore, the user should pay attention to the
#' nature of his/her data. Also, merged reflection data, having
#' to be compatible with crystal symmetry, have to display the
#' appropriate systematic absences. Users interested in keeping
#' systematic absences in the object, might want to look at the
#' object of class "raw_reflections".
#'
#' @param x Object of class 'merged_reflections'.
#' @param message A logical variable. If TRUE, the function
#'                prints a message on the errors, if any.
#' @return ans A logical value. TRUE means that the input is
#'             a valid object of class'merged_reflections'.
#'
#' @examples
#' # Create an object of class 'merged_reflections'
#' # (default ara data associated with a cubic cell)
#' x <- merged_reflections()
#'
#' # Check its validity
#' check_merged_reflections_validity(x)
#'
#' # Now change reciprocal unit cell (to triclinic)
#' uc <- unit_cell(10,20,30,30,50,70)
#' ruc <- create_rec_unit_cell(uc)
#' x$ruc <- ruc
#'
#' # Check validity again
#' check_merged_reflections_validity(x)
#'
#' @export
check_merged_reflections_validity <- function(x,message=FALSE) {
  ans <- is(x,"merged_reflections")
  if (!ans) {
    msg <- paste("This is not a valid object of class",
                 "'merged_reflections'.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Individual objects
  ans1 <- check_rec_unit_cell_validity(x$ruc)
  ans2 <- check_cryst_symm_validity(x$csym)
  if (!ans1) {
    msg <- paste("Input does not include  a valid object",
                 "of class 'rec_unit_cell'.\n")
    if (message) cat(msg)
  }
  if (!ans2) {
    msg <- paste("Input does not include a valid object",
                 "of class 'cryst_symm'.\n")
    if (message) cat(msg)
  }
  if (!ans1 | !ans2) {
    return(FALSE)
  }

  # Check records
  ans <- is(x$records,"data.frame")
  if (!ans) {
    msg <- paste("'records' is not a valid object",
                 "of class 'data.frame'.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Check 'records' has at least 3 columns
  ans <- length(x$records[1,]) >= 3
  if (!ans) {
    msg <- paste("'records' must have at least 3 columns,",
                 "the 3 Miller indices, H, K, L.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Check dtypes
  ans <- is(x$dtypes,"character")
  if (!ans) {
    msg <- paste("'dtypes' is not of a valid object",
                   "of class 'character'.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Allowed dtypes
  Vdtypes <- c("H","S","J","F","D","Q","G","L",
               "K","M","E","P","W","A","B","I","R")
  for (dt in x$dtypes) {
    if (!(dt %in% Vdtypes)) {
      msg <- "One or more 'dtypes' are not recognised.\n"
      if (message) cat(msg)
      ans <- FALSE

      return(ans)
    }
  }

  # Number of columns = length(dtypes)
  if (ncol(x$records) != length(x$dtypes)) {
    ans <- FALSE
    msg <- paste("Number of columns of 'records' must be",
                 "equal to length of 'dtypes'.\n")
    if (message) cat(msg)

    return(ans)
  }

  # First three columns have to be of dtype "H"
  ans <- (x$dtypes[1] == "H" & x$dtypes[2] == "H" &
            x$dtypes[3] == "H")
  if (!ans) {
    msg <- paste("The first three columns of 'records' has",
                 "to be of dtypes 'H'.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Unit cell compatible with symmetry?
  uc <- create_unit_cell(x$ruc)
  ans <- check_validity(uc,x$csym)
  if (!ans) {
    msg <- paste("The reciprocal unit cell of this object",
                 "is not compatible with its symmetry.\n")
    if (message) cat(msg)

    return(ans)
  }

  # Check systematic absences
  aidx <- sysabs(x$records[,1:3],x$csym$SG)
  if (length(aidx) != length(x$records[,1])) {
    ans <- FALSE
    msg <- paste("Systematic absences are presents in",
                 "'records'.\n")
    if (message) cat(msg)

    return(ans)
  }

  return(ans)
}


#####----------------------------------------------------------
#####
##### Auxiliary functions not available to users
#####
#####----------------------------------------------------------

# Compare unit cell parameters with constrains and
# return TRUE or FALSEaccording to compatibility
.ucCcsym <- function(a,b,c,aa,bb,cc,vcons) {
  # No constrains
  if (length(vcons) == 1 & vcons[1] == "No constrains") {
    return(TRUE)
  }

  # 'alpha=90'
  if ("alpha=90" %in% vcons) {
    da90 <- abs(90-aa)
    if (da90 > 0.000001) return(FALSE)
  }

  # 'beta=90'
  if ("beta=90" %in% vcons) {
    db90 <- abs(90-bb)
    if (db90 > 0.000001) return(FALSE)
  }

  # 'gamma=90'
  if ("gamma=90" %in% vcons) {
    dc90 <- abs(90-cc)
    if (dc90 > 0.000001) return(FALSE)
  }

  # 'gamma=120'
  if ("gamma=120" %in% vcons) {
    dc120 <- abs(120-cc)
    if (dc120 > 0.000001) return(FALSE)
  }

  # 'alpha=beta=gamma'
  if ("alpha=beta=gamma" %in% vcons) {
    diffab <- abs(aa-bb)
    diffbc <- abs(bb-cc)
    diffac <- abs(aa-cc)
    da90 <- abs(aa-90)
    db90 <- abs(bb-90)
    dc90 <- abs(cc-90)
    if (diffab > 0.000001 | diffac > 0.000001 | diffac > 0.000001)
      return(FALSE)
    if (da90 < 0.000001 & db90 < 0.000001 & dc90 < 0.000001)
      return(FALSE)
  }

  # 'a=b'
  if ("a=b" %in% vcons) {
    dab <- abs(a-b)
    if (dab > 0.000001) return(FALSE)
  }

  # 'a=b=c'
  if ("a=b=c" %in% vcons) {
    dab <- abs(a-b)
    dac <- abs(a-c)
    dbc <- abs(b-c)
    if (dab > 0.000001 | dac > 0.000001 | dbc > 0.000001)
      return(FALSE)
  }

  # It must be one of the above cases. Or ...
  return(TRUE)
}
