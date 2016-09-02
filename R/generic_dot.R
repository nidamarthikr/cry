#
# This file is part of the cry package
#

#
# For use in conjunction with S4-classes based modules in cRy.
# The word "dot" has been adopted because all functions here start with a ".",
# historically used for them to be invisible in the workspace.


#' Eliminates blanks from character strings
#'
#' @param wordin A character string, potentially containing blank space at
#' the beginning or the end that has to be removed.
#' @param side A character string. It can be either "both" or "left" or
#' "right", to indicate where the blank part in wordin needs to be removed.
#' @return A character string. This is the input character string deprived
#' of blank spaces at the beginning, end or on both sides of the string.
#' @examples
#' win <- "  word  "
#' out <- .strip_blanks(win)
#' print(out)
#' @export
.strip_blanks <- function(wordin,side="both")
{
 if (side != "left" & side != "right" & side != "both") stop("Input variable side can only be 'left', 'right' or 'both'")

 if (side == "left" | side == "both")
 {
  nb <- 1
  while(substr(wordin,nb,nb) == " ") nb <- nb+1
  wordin <- substr(wordin,nb,nchar(wordin))
 }

 if (side == "right" | side == "both")
 {
  nb <- nchar(wordin)
  while(substr(wordin,nb,nb) == " ") nb <- nb-1
  wordin <- substr(wordin,1,nb)
 }

 return(wordin)
}


#' Given 3 Euler angles returns 3X3 rotation matrix
#'
#' @param alpha A real numeric. The alpha Euler angle in degrees.
#' @param beta A real numeric. The beta Euler angle in degrees.
#' @param gamma A real numeric. The gamma Euler angle in degrees.
#' @return A 3X3 rotation matrix.
#' @examples
#' M <- .euler_to_matrix(0,90,0)
#' print(M)
#' M <- .euler_to_matrix(90,0,0)
#' print(M)
#' @export
.euler_to_matrix <- function(alpha,beta,gamma)
{
 # Input angles are read in degrees
 aa <- alpha*pi/180
 bb <- beta*pi/180
 cc <- gamma*pi/180
 R1 <- matrix(c(cos(cc),sin(cc),0,-sin(cc),cos(cc),0,0,0,1),nrow=3,ncol=3)
 R2 <- matrix(c(1,0,0,0,cos(bb),sin(bb),0,-sin(bb),cos(bb)),nrow=3,ncol=3)
 R3 <- matrix(c(cos(aa),sin(aa),0,-sin(aa),cos(aa),0,0,0,1),nrow=3,ncol=3)
 Rfinal <- R3%*%R2%*%R1

 return(Rfinal)
}


#' Given 3 polar angles returns 3X3 rotation matrix
#'
#' @param chi A real numeric. The chi spherical polar angle in degrees.
#' @param psi A real numeric. The psi spherical polar angle in degrees.
#' @param fi A real numeric. The fi spherical polar angle in degrees.
#' @return A 3X3 rotation matrix.
#' @examples
#' M <- .polar_to_matrix(0,90,90)
#' print(M)
#' M <- .polar_to_matrix(90,90,90)
#' print(M)
#' @export
.polar_to_matrix <- function(chi,psi,fi)
{
 # Input angles are in degrees
 chi <- chi*pi/180
 psi <- psi*pi/180
 fi <- fi*pi/180

 # Matrix components
 o11 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(cos(fi))^2
 o12 <- -sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o13 <- -sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o21 <- sin(chi)*sin(psi)*sin(fi)+(1-cos(chi))*sin(psi)*cos(psi)*cos(fi)
 o22 <- cos(chi)+(1-cos(chi))*(cos(psi))^2
 o23 <- sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o31 <- sin(chi)*cos(psi)-(1-cos(chi))*(sin(psi))^2*sin(fi)*cos(fi)
 o32 <- -sin(chi)*sin(psi)*cos(fi)-(1-cos(chi))*sin(psi)*cos(psi)*sin(fi)
 o33 <- cos(chi)+(1-cos(chi))*(sin(psi))^2*(sin(fi))^2

 # Matrix
 oline <- c(o11,o12,o13,o21,o22,o23,o31,o32,o33)
 Om <- matrix(oline,nrow=3,ncol=3,byrow=TRUE)

 return(Om)
}


#' Cross product u X v in 3D. Returns a 3 components vector
#'
#' @param u A real numeric vector of length 3.
#' @param v A real numeric vector of length 3.
#' @return A real numeric vector of length 3. This is the cross product of
#'          vectors u and v.
#' @examples
#' u <- c(1,0,0)
#' v <- c(0,1,0)
#' w <- u%X%v
#' print(w)
#' @export
"%X%" <- function(u,v)
{
 x <- u[2]*v[3]-u[3]*v[2]
 y <- u[3]*v[1]-u[1]*v[3]
 z <- u[1]*v[2]-u[2]*v[1]
 return(c(x,y,z))
}
