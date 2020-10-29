#
# This file is part of the cry package
#



#' Constructor for an S3 object of class "unit_cell.
#'
#' This represents a crystal unit cell.
#'
#' The constructor can be used with less than the full set of six input parameters,
#' to create unit cells for the cubic, tetragonal and orthogonal systems. Objects
#' of "unit_cell" class can also be created with no parameters (default to a cubic
#' cell of side length 10 angstroms).
#'
#' @param a A real number. One of the unit cell's side lengths, in angstroms.
#' @param b A real number. One of the unit cell's side lengths, in angstroms.
#' @param c A real number. One of the unit cell's side lengths, in angstroms.
#' @param aa A real number. One of the unit cell's angles, in degrees.
#' @param bb A real number. One of the unit cell's angles, in degrees.
#' @param cc A real number. One of the unit cell's angles, in degrees.
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a monoclinic unit cell
#' uc <- unit_cell(10,30,15,90,60,90)
#' print(uc)
#'
#' # Create a cubic cell (default)
#' uc <- unit_cell()
#' print(uc)
#'
#' # Create a cubic cell with side 20
#' uc <- unit_cell(20)
#' print(uc)
#'
#' # Create a tetragonal unit cell with sides 20 and 60
#' uc <- unit_cell(20,60)
#' print(uc)
#'
#' # Create an orthogonal unit cell
#' uc <- unit_cell(40,15,30)
#' print(uc)
#'
#' @export
unit_cell <- function(a=NULL,b=NULL,c=NULL,aa=NULL,bb=NULL,cc=NULL) {
  # If one of the input parameters is NULL use default values
  if (is.null(a) & is.null(b) & is.null(c) &
      is.null(aa) & is.null(bb) & is.null(cc)) {
    a <- 10
    b <- 10
    c <- 10
    aa <- 90
    bb <- 90
    cc <- 90

    # S3 class
    uc <- structure(list(a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc),class = "unit_cell")

    return(uc)
  }

  # If the first parameter only is provided: cubic cell with that side length
  if (!is.null(a) & is.null(b) & is.null(c) &
      is.null(aa) & is.null(bb) & is.null(cc)) {
    if (!is.numeric(a) | (is.numeric(a) & a <= 0)) {
      stop('The parameter provided must be a positive number.')
    }
    b <- a
    c <- a
    aa <- 90
    bb <- 90
    cc <- 90

    # S3 class
    uc <- structure(list(a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc),class = "unit_cell")

    return(uc)
  }

  # If the first and second parameter only are provided: tetragonal cell
  if (!is.null(a) & !is.null(b) & is.null(c) &
      is.null(aa) & is.null(bb) & is.null(cc)) {
    if (!is.numeric(a) | !is.numeric(b)) {
      stop('The two parameters provided must be positive numbers.')
    }
    if (a <= 0 | b <= 0) {
      stop('The two parameters provided must be positive numbers.')
    }
    c <- b
    b <- a
    aa <- 90
    bb <- 90
    cc <- 90

    # S3 class
    uc <- structure(list(a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc),class = "unit_cell")

    return(uc)
  }

  # If the first, second and third parameter only are provided: orthogonal cell
  if (!is.null(a) & !is.null(b) & !is.null(c) &
      is.null(aa) & is.null(bb) & is.null(cc)) {
    if (!is.numeric(a) | !is.numeric(b) | !is.numeric(c)) {
      stop('The three parameters provided must be positive numbers.')
    }
    if (a <= 0 | b <= 0 | c <= 0) {
      stop('The three parameters provided must be positive numbers.')
    }
    aa <- 90
    bb <- 90
    cc <- 90

    # S3 class
    uc <- structure(list(a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc),class = "unit_cell")

    return(uc)
  }

  # General case
  # Check input types
  if (!is.numeric(a) | !is.numeric(b) | !is.numeric(c) |
      !is.numeric(aa) | !is.numeric(bb) | !is.numeric(cc)) {
    stop("All six input values must be numerics.")
  }

  # Check sides
  if (a <= 0 | b <= 0 | c <= 0) {
    stop("A unit cell must have positive sides.")
  }

  # Check angles. # Angles have values between 0 and 180 degrees. But combinations
  # of alpha beta and gamma need to obey certain rules (see J. Foadi & G. Evans (2011))
  if (aa < 0 | bb < 0 | cc < 0 | aa > 180 | bb > 180 | cc > 180) {
    stop("A unit cell must have angles between 0 and 180 degrees.")
  }
  ss <- aa+bb+cc
  if (ss >= 360) stop("A unit cell with these angles cannot exist.")
  ss <- aa+bb-cc
  if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist.")
  ss <- aa-bb+cc
  if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist.")
  ss <- -aa+bb+cc
  if (ss <= 0 | ss >= 360) stop("A unit cell with these angles cannot exist.")

  # Turn angles to objects of class "angle"
  aa <- angle(aa)
  bb <- angle(bb)
  cc <- angle(cc)

  # S3 class
  uc <- structure(list(a=a,b=b,c=c,alpha=aa,beta=bb,gamma=cc),class = "unit_cell")

  return(uc)
}


#' Constructor for an S3 object of class "rec_unit_cell".
#'
#' This represents a crystal reciprocal unit cell.
#'
#' @param ar A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param br A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param cr A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param aar A real number. One of the reciprocal unit cell's angles, in degrees.
#' @param bbr A real number. One of the reciprocal unit cell's angles, in degrees.
#' @param ccr A real number. One of the reciprocal unit cell's angles, in degrees.
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a cubic reciprocal unit cell
#' ruc1 <- rec_unit_cell(1/10,1/10,1/10,90,90,90)
#' class(ruc1)
#' @export
rec_unit_cell <- function(ar,br,cr,aar,bbr,ccr) {
  # Check input types
  if (!is.numeric(ar) | !is.numeric(br) | !is.numeric(cr) |
      !is.numeric(aar) | !is.numeric(bbr) | !is.numeric(ccr)) {
    stop("All six input values must be numerics.")
  }

  # Check sides
  if (ar <= 0 | br <= 0 | cr <= 0) {
    stop("A reciprocal unit cell must have positive sides.")
  }

  # Check angles. # Angles have values between 0 and 180 degrees. But combinations
  # of alpha beta and gamma need to obey certain rules (see J. Foadi & G. Evans (2011))
  if (aar < 0 | bbr < 0 | ccr < 0 | aar > 180 | bbr > 180 | ccr > 180) {
    stop("A reciprocal unit cell must have angles between 0 and 180 degrees.")
  }
  ss <- aar+bbr+ccr
  if (ss >= 360) stop("A reciprocal unit cell with these angles cannot exist.")
  ss <- aar+bbr-ccr
  if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist.")
  ss <- aar-bbr+ccr
  if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist.")
  ss <- -aar+bbr+ccr
  if (ss <= 0 | ss >= 360) stop("A reciprocal unit cell with these angles cannot exist.")

  # Turn angles to objects of class "angle"
  aar <- angle(aar)
  bbr <- angle(bbr)
  ccr <- angle(ccr)

  # S3 class
  ruc <- structure(list(ar=ar,br=br,cr=cr,alphar=aar,betar=bbr,gammar=ccr),class = "rec_unit_cell")

  return(ruc)
}


#' Print method for an object of class "unit_cell".
#'
#' The values are displayed in angstroms and degrees
#'
#' @param x An object of class "unit_cell"
#' @examples
#' # Create a cubic unit cell
#' uc <- unit_cell(10,10,10,90,90,,90)
#'
#' # Display its value
#' print(uc)
#' @export
print.unit_cell <- function(x) {
  cat("This is an object of class 'unit_cell'\n")
  cat("Its sides are:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    angstroms.\n",x$a,x$b,x$c)
  cat(msg)
  cat("Its angles are:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    degrees.\n",x$alpha[1],x$beta[1],x$gamma[1])
  cat(msg)
  invisible(x)
}


#' Print method for an object of class "rec_unit_cell".
#'
#' The values are displayed in 1/angstroms and degrees
#'
#' @param x An object of class "rec_unit_cell"
#' @examples
#' # Create a cubic reciprocal unit cell
#' ruc <- rec_unit_cell(1/10,1/10,1/10,90,90,,90)
#'
#' # Display its value
#' print(ruc)
#' @export
print.rec_unit_cell <- function(x) {
  cat("This is an object of class 'rec_unit_cell'\n")
  cat("Its sides are:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    1/angstroms.\n",x$ar,x$br,x$cr)
  cat(msg)
  cat("Its angles are:\n")
  msg <- sprintf("  %10.3f , %10.3f , %10.3f    degrees.\n",x$alphar[1],x$betar[1],x$gammar[1])
  cat(msg)
  invisible(x)
}


#' Default method for generic "create_unit_cell"
#'
#' This method is an alternative call to "unit_cell"
#'.
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a cubic cell with side 50
#' uc <- create_unit_cell(50)
#' print(uc)
#'
#' @seealso
#' \code{\link{unit_cell}}
#'
#' @export
create_unit_cell.default <- function(a=NULL,b=NULL,c=NULL,aa=NULL,bb=NULL,cc=NULL,...) {
  # Makes use of all constructions in "unit_cell"
  uc <- unit_cell(a,b,c,aa,bb,cc)

  return(uc)
}


#' Method to create a unit_cell object starting from a Bravais symbol
#'
#' The Bravais symbols indicate the 14 possible Bravais lattices. A few
#' examples are "aP", "oF", etc. The cell parameters assigned are assigned
#' randomly, but are compatible with the Bravais lattice.
#'
#' @param bt An object of class "bravais".
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "unit_cell" object from a monoclinic primitive Bravais lattice
#' # Cell parameters generated automatically.
#' bt <- bravais("mP")
#' uc <- create_unit_cell(bt)
#' print(uc)
#'
#' # Create a "unit_cell" object from a monoclinic primitive Bravais lattice
#' # Cell parameters fixed by user.
#' bt <- bravais("mP")
#' uc <- create_unit_cell(bt,c(10,30,20,90,105,90))
#' print(uc)
#' @export
create_unit_cell.bravais <- function(bt,cpar=NULL) {
  # Check bt is an object of class "bravais"
  if(!is(bt,"bravais")) {
    stop("Input must be a valid object of class 'bravais'.")
  }

  # Extract lattice and centering (also check on symbol validity)
  latt <- substr(bt$bt,1,1)
  centring <- substr(bt$bt,2,2)

  if (!is.null(cpar)) {
    if (!is.numeric(cpar)) stop('Second argument is a 6 cell-parameters vector.')
    if (length(cpar) != 6) stop('Second argument is a 6 cell-parameters vector.')
    a <- cpar[1]
    b <- cpar[2]
    c <- cpar[3]
    aa <- cpar[4]
    bb <- cpar[5]
    cc <- cpar[6]

    # Are these parameters compatible with the Bravais lattice?
    if (latt == "m") {
      erang <- 0.000001      # Finite accuracy of binary representation of numbers means we have to
      diff_a <- abs(aa-90)   # test number == 90 in this way
      diff_b <- abs(bb-90)
      diff_c <- abs(cc-90)
      if (diff_a > erang & (diff_b > erang | diff_c > erang))
        stop ("Unit cell angles non compatible with lattice type")
      if (diff_b > erang & (diff_a > erang | diff_c > erang))
        stop ("Unit cell angles non compatible with lattice type")
      if (diff_c > erang & (diff_a > erang | diff_b > erang))
        stop ("Unit cell angles non compatible with lattice type")
      if (diff_a > erang & centring == "A")
        stop ("Unit cell parameters non compatible with lattice type")
      if (diff_b > erang & centring == "B")
        stop ("Unit cell parameters non compatible with lattice type")
      if (diff_c > erang & centring == "C")
        stop ("Unit cell parameters non compatible with lattice type")
    }
    if (latt == "o") {
      if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
        stop ("Unit cell parameters non compatible with lattice type")
    }
    if (latt == "t") {
      if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
        stop ("Unit cell parameters non compatible with lattice type")
      if (abs(a-b) > 0.000001 & abs(a-c) > 0.000001 & abs(b-c) > 0.000001)
        stop ("Unit cell parameters non compatible with lattice type")
    }
    if (latt == "c") {
      if (abs(a-b) > 0.000001 | abs(a-c) > 0.000001 | abs(b-c) > 0.000001)
        stop ("Unit cell parameters non compatible with lattice type")
      if (abs(aa-90) > 0.000001 | abs(bb-90) > 0.000001 | abs(cc-90) > 0.000001)
        stop ("Unit cell parameters non compatible with lattice type")
    }
    if (latt == "h") {
      if (abs(a-b) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
      if (abs(a-b) < 0.000001) {
        if (abs(cc-120) > 0.000001) {
          if (abs(b-c) > 0.000001) stop ("Unit cell parameters non compatible with lattice type")
        }
      }
    }
    # All check passed. Create (further down) unit_cell object.
  } else {
    # Selection of arbitrary unit cell parameters
    sidebox <- seq(10,100,by=10)
    anglebox <- seq(5,175,by=5)
    anglebox2 <- seq(5,85,by=5)

    # Next, assign cell parameters according to Bravais lattice
    # Triclinic
    if (latt == "a") {
      a <- sample(sidebox,1)
      b <- sample(sidebox,1)
      c <- sample(sidebox,1)
      ans <- FALSE
      while (!ans) {
        aa <- sample(anglebox,1)
        bb <- sample(anglebox,1)
        cc <- sample(anglebox,1)
        ans <- TRUE
        if (aa < 0 | bb < 0 | cc < 0 | aa > 180 | bb > 180 | cc > 180) ans <- FALSE
        ss <- aa+bb+cc
        if (ss >= 360) ans <- FALSE
        ss <- aa+bb-cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
        ss <- aa-bb+cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
        ss <- -aa+bb+cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
      }
    }

    # Monoclinic
    if (latt == "m") {
      a <- sample(sidebox,1)
      b <- sample(sidebox,1)
      c <- sample(sidebox,1)
      ans <- FALSE
      while (!ans) {
        aa <- 90
        bb <- sample(anglebox,1)
        cc <- 90
        ans <- TRUE
        if (aa < 0 | bb < 0 | cc < 0 | aa > 180 | bb > 180 | cc > 180) ans <- FALSE
        ss <- aa+bb+cc
        if (ss >= 360) ans <- FALSE
        ss <- aa+bb-cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
        ss <- aa-bb+cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
        ss <- -aa+bb+cc
        if (ss <= 0 | ss >= 360) ans <- FALSE
      }
    }

    # Orthorombic
    if (latt == "o") {
      a <- sample(sidebox,1)
      b <- sample(sidebox,1)
      c <- sample(sidebox,1)
      aa <- 90
      bb <- 90
      cc <- 90
    }

    # Tetragonal
    if (latt == "t") {
      a <- sample(sidebox,1)
      b <- a
      c <- sample(sidebox,1)
      aa <- 90
      bb <- 90
      cc <- 90
    }

    # Hexagonal
    if (latt == "h" & centring == "P") {
      a <- sample(sidebox,1)
      b <- a
      c <- sample(sidebox,1)
      aa <- 90
      bb <- 90
      cc <- 120
    }

    # Rombohedral
    if (latt == "h" & centring == "R") {
      a <- sample(sidebox,1)
      b <- a
      c <- a
      aa <- sample(anglebox2,1)
      bb <- aa
      cc <- aa
    }

    # Cubic
    if (latt == "c") {
      a <- sample(sidebox,1)
      b <- a
      c <- a
      aa <- 90
      bb <- 90
      cc <- 90
    }
  }

  # S3 object
  uc <- unit_cell(a,b,c,aa,bb,cc)

  return(uc)
}


#' Reciprocal unit cell starting from a unit cell
#'
#' Method to create an object of class "rec_unit_cell" starting from an object of
#' class "unit_cell".
#'
#' @param uc An object of class "unit_cell".
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "rec_unit_cell" object starting from a cubic cell object
#' uc <- unit_cell()
#' print(uc)
#' ruc <- create_rec_unit_cell(uc)
#' print(ruc)
#'
#' @export
create_rec_unit_cell.unit_cell <- function(uc) {
  print(uc)
}
