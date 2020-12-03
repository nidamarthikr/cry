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
    aa <- angle(90)
    bb <- angle(90)
    cc <- angle(90)

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
    aa <- angle(90)
    bb <- angle(90)
    cc <- angle(90)

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
    aa <- angle(90)
    bb <- angle(90)
    cc <- angle(90)

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
    aa <- angle(90)
    bb <- angle(90)
    cc <- angle(90)

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


#' Constructor for an S3 object of class "rec_unit_cell.
#'
#' This represents a crystal reciprocal unit cell.
#'
#' The constructor can be used with less than the full set of six input parameters,
#' to create reciprocal unit cells for the cubic, tetragonal and orthogonal systems. Objects
#' of "rec_unit_cell" class can also be created with no parameters (default to a reciprocal cubic
#' cell of side length 0.1 1/angstroms).
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
#' # Create a monoclinic reciprocal unit cell
#' ruc <- unit_cell(0.115,0.033,0.077,90,120,90)
#' print(ruc)
#'
#' # Create a cubic cell (default)
#' ruc <- rec_unit_cell()
#' print(ruc)
#'
#' # Create a reciprocal cubic cell with side 1/20
#' ruc <- rec_unit_cell(1/20)
#' print(ruc)
#'
#' # Create a reciprocal tetragonal unit cell with sides 1/20 and 1/60
#' ruc <- rec_unit_cell(1/20,1/60)
#' print(ruc)
#'
#' # Create a reciprocal orthogonal unit cell
#' ruc <- rec_unit_cell(1/40,1/15,1/30)
#' print(ruc)
#'
#' @export
rec_unit_cell <- function(ar=NULL,br=NULL,cr=NULL,aar=NULL,bbr=NULL,ccr=NULL) {
   # If one of the input parameters is NULL use default values
   if (is.null(ar) & is.null(br) & is.null(cr) &
      is.null(aar) & is.null(bbr) & is.null(ccr)) {
     ar <- 0.1
     br <- 0.1
     cr <- 0.1
     aar <- angle(90)
     bbr <- angle(90)
     ccr <- angle(90)

     # S3 class
     ruc <- structure(list(ar=ar,br=br,cr=cr,alphar=aar,betar=bbr,gammar=ccr),class = "rec_unit_cell")

     return(ruc)
   }

  # If the first parameter only is provided: cubic cell with that side length
  if (!is.null(ar) & is.null(br) & is.null(cr) &
      is.null(aar) & is.null(bbr) & is.null(ccr)) {
    if (!is.numeric(ar) | (is.numeric(ar) & ar <= 0)) {
      stop('The parameter provided must be a positive number.')
    }
    br <- ar
    cr <- ar
    aar <- angle(90)
    bbr <- angle(90)
    ccr <- angle(90)

    # S3 class
    ruc <- structure(list(ar=ar,br=br,cr=cr,alphar=aar,betar=bbr,gammar=ccr),class = "rec_unit_cell")

    return(ruc)
  }

  # If the first and second parameter only are provided: tetragonal cell
  if (!is.null(ar) & !is.null(br) & is.null(cr) &
      is.null(aar) & is.null(bbr) & is.null(ccr)) {
    if (!is.numeric(ar) | !is.numeric(br)) {
      stop('The two parameters provided must be positive numbers.')
    }
    if (ar <= 0 | br <= 0) {
      stop('The two parameters provided must be positive numbers.')
    }
    cr <- br
    br <- ar
    aar <- angle(90)
    bbr <- angle(90)
    ccr <- angle(90)

    # S3 class
    ruc <- structure(list(ar=ar,br=br,cr=cr,alphar=aar,betar=bbr,gammar=ccr),class = "rec_unit_cell")

    return(ruc)
  }

  # If the first, second and third parameter only are provided: orthogonal cell
  if (!is.null(ar) & !is.null(br) & !is.null(cr) &
      is.null(aar) & is.null(bbr) & is.null(ccr)) {
    if (!is.numeric(ar) | !is.numeric(br) | !is.numeric(cr)) {
      stop('The three parameters provided must be positive numbers.')
    }
    if (ar <= 0 | br <= 0 | cr <= 0) {
      stop('The three parameters provided must be positive numbers.')
    }
    aar <- angle(90)
    bbr <- angle(90)
    ccr <- angle(90)

    # S3 class
    ruc <- structure(list(ar=ar,br=br,cr=cr,alphar=aar,betar=bbr,gammar=ccr),class = "rec_unit_cell")

    return(ruc)
  }

  # General case
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
#' @param x An object of class "unit_cell".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create a cubic unit cell
#' uc <- unit_cell(10,10,10,90,90,90)
#'
#' # Display its value
#' print(uc)
#' @rdname print.unit_cell
#' @export
print.unit_cell <- function(x,...) {
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
#' @param x An object of class "rec_unit_cell".
#' @param ... Additional arguments passed to the print methods
#' @examples
#' # Create a cubic reciprocal unit cell
#' ruc <- rec_unit_cell(1/10,1/10,1/10,90,90,90)
#'
#' # Display its value
#' print(ruc)
#' @rdname print.rec_unit_cell
#' @export
print.rec_unit_cell <- function(x,...) {
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
#' @param a A real number. One of the unit cell's side lengths, in angstroms.
#' @param b A real number. One of the unit cell's side lengths, in angstroms.
#' @param c A real number. One of the unit cell's side lengths, in angstroms.
#' @param aa A real number. One of the unit cell's angles, in degrees.
#' @param bb A real number. One of the unit cell's angles, in degrees.
#' @param cc A real number. One of the unit cell's angles, in degrees.
#' @param ... Additional arguments passed to the create_unit_cell methods
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
#' @rdname create_unit_cell.default
#' @export
create_unit_cell.default <- function(a=NULL,b=NULL,c=NULL,aa=NULL,bb=NULL,cc=NULL,...) {
  # Makes use of all constructions in "unit_cell"
  uc <- unit_cell(a,b,c,aa,bb,cc)

  return(uc)
}


#' Unit cell starting from a Bravais symbol
#'
#' Method to create a "unit_cell" object starting from a "bravais" object.
#' The Bravais symbols indicate the 14 possible Bravais lattices. A few
#' examples are "aP", "oF", etc. The cell parameters assigned are assigned
#' randomly, but are compatible with the Bravais lattice.
#'
#' @param a An object of class "bravais".
#' @param ... Additional arguments passed to the create_unit_cell methods
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "unit_cell" object from a monoclinic primitive Bravais lattice
#' # Cell parameters generated automatically.
#' bt <- bravais("mP")
#' uc <- create_unit_cell(bt)
#' print(uc)
#'
#' @rdname create_unit_cell.bravais
#' @export
create_unit_cell.bravais <- function(a,...) {
  bt <- a
  # Check bt is an object of class "bravais"
  if(!is(bt,"bravais")) {
    stop("Input must be a valid object of class 'bravais'.")
  }

  # Extract lattice and centering (also check on symbol validity)
  latt <- substr(bt$bt,1,1)
  centring <- substr(bt$bt,2,2)

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

  # S3 object
  uc <- unit_cell(a,b,c,aa,bb,cc)

  return(uc)
}


#' Unit cell starting from a reciprocal unit cell
#'
#' Method to create an object of class "unit_cell" starting from an object of
#' class "rec_unit_cell".
#'
#' @param a An object of class "rec_unit_cell".
#' @param ... Additional arguments passed to the create_unit_cell methods
#' @return An object of class "unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "unit_cell" object starting from a reciprocal cubic cell object
#' ruc <- rec_unit_cell()
#' print(ruc)
#' uc <- create_unit_cell(ruc)
#' print(uc)
#'
#' @rdname create_unit_cell.rec_unit_cell
#' @export
create_unit_cell.rec_unit_cell <- function(a,...) {
  ruc <- a
  # Check ruc is an object of class "rec_unit_cell"
  if(!is(ruc,"rec_unit_cell")) {
    stop("Input must be a valid object of class 'rec_unit_cell'.")
  }

  # Extract reciprocal unit cell parameters
  ar <- ruc$ar
  br <- ruc$br
  cr <- ruc$cr
  aar <- ruc$alphar[1]
  bbr <- ruc$betar[1]
  ccr <- ruc$gammar[1]

  # Lattice calculations
  ltmp <- lattice_stuff(ar,br,cr,aar,bbr,ccr)

  # Cell parameters
  a <- ltmp[[7]]
  b <- ltmp[[8]]
  c <- ltmp[[9]]
  aa <- atan2(ltmp[[10]],ltmp[[13]])*180/pi
  bb  <- atan2(ltmp[[11]],ltmp[[14]])*180/pi
  cc <- atan2(ltmp[[12]],ltmp[[15]])*180/pi

  # unit_cell_object
  uc <- unit_cell(a,b,c,aa,bb,cc)

  return(uc)
}


#' Volume of a unit cell (in angstroms^3)
#'
#' Method of the S3 generic class "calculate_cell_volume", to calculate the
#' volume, in cubic angstroms, of the unit cell corresponding to the input
#' object of class "unit_cell".
#'
#' @param x An object of class "unit_cell".
#' @param ... Additional arguments passed to the calculate_cell_volume methods
#' @return A positive numeric, the volume in cubic angstroms of the unit cell
#'  corresponding to the input.
#' @examples
#' # Create a monoclinic cell
#' bt <- bravais("mP")
#' uc <- create_unit_cell(bt)
#' print(uc)
#'
#' # Calculate cell volume
#' calculate_cell_volume(uc)
#'
#' @rdname calculate_cell_volume.unit_cell
#' @export
calculate_cell_volume.unit_cell <- function(x,...) {
  uc <- x
  # Check uc is an object of class "unit_cell"
  if(!is(uc,"unit_cell")) {
    stop("Input must be a valid object of class 'unit_cell'.")
  }

  # Extract unit cell parameters
  a <- uc$a
  b <- uc$b
  c <- uc$c
  aa <- uc$alpha[1]
  bb <- uc$beta[1]
  cc <- uc$gamma[1]

  # Lattice calculations
  ltmp <- lattice_stuff(a,b,c,aa,bb,cc)

  # Volume
  V <- ltmp[[16]]

  return(V)
}


#' Default method for generic "create_rec_unit_cell"
#'
#' This method is an alternative call to "rec_unit_cell"
#'.
#' @param ar A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param br A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param cr A real number. One of the reciprocal unit cell's side lengths, in 1/angstroms.
#' @param aar A real number. One of the reciprocal unit cell's angles, in degrees.
#' @param bbr A real number. One of the reciprocal unit cell's angles, in degrees.
#' @param ccr A real number. One of the reciprocal unit cell's angles, in degrees.
#' @param ... Additional arguments passed to the create_rec_unit_cell methods
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a reciprocal cubic cell with side 1/50
#' ruc <- create_rec_unit_cell(1/50)
#' print(ruc)
#'
#' @seealso
#' \code{\link{rec_unit_cell}}
#'
#' @rdname create_rec_unit_cell.default
#' @export
create_rec_unit_cell.default <- function(ar=NULL,br=NULL,cr=NULL,aar=NULL,bbr=NULL,ccr=NULL,...) {
  # Makes use of all constructions in "rec_unit_cell"
  ruc <- rec_unit_cell(ar,br,cr,aar,bbr,ccr)

  return(ruc)
}


#' Reciprocal unit cell starting from a Bravais symbol
#'
#' Method to create a "rec_unit_cell" object starting from a "bravais" object.
#' The Bravais symbols indicate the 14 possible Bravais lattices. A few
#' examples are "aP", "oF", etc. The cell parameters assigned are assigned
#' randomly, but are compatible with the Bravais lattice.
#'
#' @param ar An object of class "bravais".
#' @param ... Additional arguments passed to the create_rec_unit_cell methods
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "rec_unit_cell" object from a monoclinic primitive Bravais lattice
#' # Cell parameters generated automatically.
#' bt <- bravais("mP")
#' ruc <- create_rec_unit_cell(bt)
#' print(ruc)
#'
#' @rdname create_rec_unit_cell.bravais
#' @export
create_rec_unit_cell.bravais <- function(ar,...) {
  bt <- ar
  # Check bt is an object of class "bravais"
  if(!is(bt,"bravais")) {
    stop("Input must be a valid object of class 'bravais'.")
  }

  # Extract lattice and centering (also check on symbol validity)
  latt <- substr(bt$bt,1,1)
  centring <- substr(bt$bt,2,2)

  # Selection of arbitrary reciprocal unit cell parameters
  sidebox <- seq(10,100,by=10)
  sidebox <- 1/sidebox
  anglebox <- seq(5,175,by=5)
  anglebox2 <- seq(5,85,by=5)

  # Next, assign cell parameters according to Bravais lattice
  # Triclinic
  if (latt == "a") {
    ar <- sample(sidebox,1)
    br <- sample(sidebox,1)
    cr <- sample(sidebox,1)
    ans <- FALSE
    while (!ans) {
      aar <- sample(anglebox,1)
      bbr <- sample(anglebox,1)
      ccr <- sample(anglebox,1)
      ans <- TRUE
      if (aar < 0 | bbr < 0 | ccr < 0 | aar > 180 | bbr > 180 | ccr > 180) ans <- FALSE
      ss <- aar+bbr+ccr
      if (ss >= 360) ans <- FALSE
      ss <- aar+bbr-ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
      ss <- aar-bbr+ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
      ss <- -aar+bbr+ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
    }
  }

  # Monoclinic
  if (latt == "m") {
    ar <- sample(sidebox,1)
    br <- sample(sidebox,1)
    cr <- sample(sidebox,1)
    ans <- FALSE
    while (!ans) {
      aar <- 90
      bbr <- sample(anglebox,1)
      ccr <- 90
      ans <- TRUE
      if (aar < 0 | bbr < 0 | ccr < 0 | aar > 180 | bbr > 180 | ccr > 180) ans <- FALSE
      ss <- aar+bbr+ccr
      if (ss >= 360) ans <- FALSE
      ss <- aar+bbr-ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
      ss <- aar-bbr+ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
      ss <- -aar+bbr+ccr
      if (ss <= 0 | ss >= 360) ans <- FALSE
    }
  }

  # Orthorombic
  if (latt == "o") {
    ar <- sample(sidebox,1)
    br <- sample(sidebox,1)
    cr <- sample(sidebox,1)
    aar <- 90
    bbr <- 90
    ccr <- 90
  }

  # Tetragonal
  if (latt == "t") {
    ar <- sample(sidebox,1)
    br <- ar
    cr <- sample(sidebox,1)
    aar <- 90
    bbr <- 90
    ccr <- 90
  }

  # Hexagonal
  if (latt == "h" & centring == "P") {
    ar <- sample(sidebox,1)
    br <- ar
    cr <- sample(sidebox,1)
    aar <- 90
    bbr <- 90
    ccr <- 120
  }

  # Rombohedral
  if (latt == "h" & centring == "R") {
    ar <- sample(sidebox,1)
    br <- ar
    cr <- ar
    aar <- sample(anglebox2,1)
    bbr <- aar
    ccr <- aar
  }

  # Cubic
  if (latt == "c") {
    ar <- sample(sidebox,1)
    br <- ar
    cr <- ar
    aar <- 90
    bbr <- 90
    ccr <- 90
  }

  # S3 object
  ruc <- rec_unit_cell(ar,br,cr,aar,bbr,ccr)

  return(ruc)
}


#' Reciprocal unit cell starting from a unit cell
#'
#' Method to create an object of class "rec_unit_cell" starting from an object of
#' class "unit_cell".
#'
#' @param ar An object of class "unit_cell".
#' @param ... Additional arguments passed to the create_rec_unit_cell methods
#' @return An object of class "rec_unit_cell". It is a named list of length 6 whose
#'         last three slots are of "angle" class.
#' @examples
#' # Create a "rec_unit_cell" object starting from a cubic cell object
#' uc <- unit_cell()
#' print(uc)
#' ruc <- create_rec_unit_cell(uc)
#' print(ruc)
#'
#' @rdname create_rec_unit_cell.unit_cell
#' @export
create_rec_unit_cell.unit_cell <- function(ar,...) {
  uc <- ar
  # Check uc is an object of class "unit_cell"
  if(!is(uc,"unit_cell")) {
    stop("Input must be a valid object of class 'unit_cell'.")
  }

  # Extract unit cell parameters
  a <- uc$a
  b <- uc$b
  c <- uc$c
  aa <- uc$alpha[1]
  bb <- uc$beta[1]
  cc <- uc$gamma[1]

  # Lattice calculations
  ltmp <- lattice_stuff(a,b,c,aa,bb,cc)

  # Reciprocal cell parameters
  ar <- ltmp[[7]]
  br <- ltmp[[8]]
  cr <- ltmp[[9]]
  aar <- atan2(ltmp[[10]],ltmp[[13]])*180/pi
  bbr  <- atan2(ltmp[[11]],ltmp[[14]])*180/pi
  ccr <- atan2(ltmp[[12]],ltmp[[15]])*180/pi

  # rec_unit_cell_object
  ruc <- rec_unit_cell(ar,br,cr,aar,bbr,ccr)

  return(ruc)
}


#' Volume of a reciprocal unit cell (in angstroms^(-3))
#'
#' Method of the S3 generic class "calculate_cell_volume", to calculate the
#' volume, in reciprocal cubic angstroms, of the reciprocal unit cell corresponding
#' to the input object of class "rec_unit_cell".
#'
#' @param x An object of class "rec_unit_cell".
#' @param ... Additional arguments passed to the calculate_cell_volume methods
#' @return A positive numeric, the volume in reciprocal cubic angstroms of the
#'  reciprocal unit cell corresponding to the input.
#' @examples
#' # Create a monoclinic cell and the corresponding reciprocal
#' bt <- bravais("mP")
#' uc <- create_unit_cell(bt)
#' ruc <- create_rec_unit_cell(uc)
#'
#' # Calculate reciprocall cell volume
#' calculate_cell_volume(ruc)
#'
#' @rdname calculate_cell_volume.rec_unit_cell
#' @export
calculate_cell_volume.rec_unit_cell <- function(x,...) {
  ruc <- x
  # Check ruc is an object of class "rec_unit_cell"
  if(!is(ruc,"rec_unit_cell")) {
    stop("Input must be a valid object of class 'rec_unit_cell'.")
  }

  # Extract reciprocal cell parameters
  ar <- ruc$ar
  br <- ruc$br
  cr <- ruc$cr
  aar <- ruc$alphar[1]
  bbr <- ruc$betar[1]
  ccr <- ruc$gammar[1]

  # Lattice calculations
  ltmp <- lattice_stuff(ar,br,cr,aar,bbr,ccr)

  # Volume
  Vr <- ltmp[[16]]

  return(Vr)
}
