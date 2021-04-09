#
# This file is part of the cry package
#


#' Translation of space group symbols, numbers, etc.
#'
#' Function to find out space group symbol given number and vice-versa.
#'
#' This function returns either a number of a specific symbol corresponding to a
#' crystallographic space group. The input is an integer number or a character symbol
#' identifying a specific space group. The output is, similarly, the corresponding
#' character symbol or number, according to what is specified in the input.
#' Possible formats are:
#' \itemize{
#'   \item{1) Space group number}
#'   \item{2) Hall symbol (e.g. ' P 2yb (z,x,y)')}
#'   \item{3) Extended Hermann-Maguin symbol (e.g. 'P 1 1 21')}
#' }
#' If more than one setting is implied in an ambiguous way in the input value,
#' then the first setting will be selected by default for the output value, unless
#' argument "set" is set to another value.
#'
#' @param value A string or an integer number corresponding to the space group being
#'  investigated.
#' @param SG_in A string representing the space group format for the input. Possible values are:
#'  \itemize{
#'    \item{1) "number"}
#'    \item{2) "ccp4"}
#'    \item{3) "Hall"}
#'    \item{4) "xHM"}
#'    \item{5) "old"}
#'  }
#' @param SG_out A string representing the space group format for the output. Possible values are:
#'  \itemize{
#'    \item{1) "number"}
#'    \item{2) "ccp4"}
#'    \item{3) "Hall"}
#'    \item{4) "xHM"}
#'    \item{5) "old"}
#'  }
#' @param set Specific setting for the given space group. A number like 1,2,...
#'   It is used if for a same symbol there are more than one choice.
#' @return list_SG A named list with two fields. The first field, "msg", is a character string
#'   representing the space group format needed as output. Possible values are the same as
#'   those for SG_in. The second field, "ans", is TRUE only if a valid symbol for "msg" is
#'   found.
#'
#' @examples
#' # Space Group P1 corresponds to number 1
#' translate_SG(value=1,SG_in="number",SG_out="xHM")
#'
#' @export
translate_SG <- function(value,SG_in="number",SG_out="xHM",set=1) {
  # Special case in which SG_in is not a number
  #if (SG_in != "number") {
  #  tmp <- .translate_SG(value,SG_in,SG_out="number",set=set)
  #  if (!tmp$ans) {
  #    list_SG <- list(msg=tmp$msg,ans=FALSE)

  #    return(list_SG)
  #  }
  #  tmp <- .translate_SG(tmp$msg,SG_in="number",SG_out=SG_out,
  #                       set=set)
  #  tmp <- .translate_SG(value,SG_in,SG_out,set)
  #  list_SG <- list(msg=tmp$msg,ans=tmp$ans)

  #  return(list_SG)
  #}
  list_SG <- .translate_SG(value,SG_in,SG_out,set)

  return(list_SG)
}



#' Information on a specific space group
#'
#' Returns human-readable information on a specific input space group.
#'
#' Crystallographic symmetry is fundamental in crystallography. It affects the way
#' atoms are arranged in a unit cell, the pattern of reflections in reciprocal space
#' and many other common occurrences in crystallography. This function returns a named
#' list with human-readable character strings which detail key symmetry information.
#'
#' @param SG A character string. The extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
#' @return infostring A named list with fields corresponding to those in the CCP4 symmetry
#'  library. The fields' name are:
#'  \itemize{
#'    \item{\strong{NUMBER} standard spacegroup number}
#'    \item{\strong{BASISOP} change of basis operator}
#'    \item{\strong{CCP4} CCP4 spacegroup number e.g. 1003 (0 if not a CCP4 group)}
#'    \item{\strong{HALL} Hall symbol}
#'    \item{\strong{xHM} extended Hermann Mauguin symbol}
#'    \item{\strong{OLD} CCP4 spacegroup name (blank if not a CCP4 group)}
#'    \item{\strong{LAUE} Laue group symbol}
#'    \item{\strong{PATT} Patterson group symbol}
#'    \item{\strong{PGRP} Point group symbol}
#'    \item{\strong{HKLASU} reciprocal space asymmetric unit (with respect to standard setting)}
#'    \item{\strong{MAPASU_CCP4} CCP4 real space asymmetric unit (with respect to standard setting.
#'     Negative ranges if not a CCP4 group)}
#'    \item{\strong{MAPASU_ZERO} origin based real space asymmetric unit (with respect to current
#'     setting)}
#'    \item{\strong{MAPASU_NONZ} non-origin based real space asymmetric uni (with respect to
#'     current setting)}
#'    \item{\strong{CESHIRE} Cheshire cell (with respect to standard setting)}
#'    \item{\strong{SYMOP} list of primitive symmetry operators}
#'    \item{\strong{CENOP} list of centering operators}
#'  }
#'
#' @examples
#' # This is the full information for space group number 19, P 21 21 21
#' SG <- translate_SG(19)
#' ltmp <- extract_symmetry_info(SG)
#' ltmp
#'
#' @export
extract_symmetry_info <- function(SG) {
  # Input is spacegroup symbol in xHM format.
  # Output is a list with all symmetry information
  # for the specific SG group, as contained in syminfo.lib.

  # The whole content of "syminfo.lib" is already in object syminfo
  bsg <- grep(SG,syminfo)
  bini <- grep("begin_spacegroup",syminfo)
  bend <- grep("end_spacegroup",syminfo)
  prima <- bini[length(bini[bini < bsg[1]])]
  dopo  <- bend[bend > bsg[1]][1]

  # If symbol used is not valid, return NULL
  if (is.na(prima) | is.na(dopo)) {
    infostring <- NULL
    return(infostring)
  }

  # Extract info in string format
  infostring <- list()
  infostring$NUMBER       <- syminfo[prima+ 1]
  infostring$BASISOP      <- syminfo[prima+ 2]
  infostring$CCP4         <- syminfo[prima+ 3]
  infostring$HALL         <- syminfo[prima+ 4]
  infostring$XHM          <- syminfo[prima+ 5]
  infostring$OLD          <- syminfo[prima+ 6]
  infostring$LAUE         <- syminfo[prima+ 7]
  infostring$PATT         <- syminfo[prima+ 8]
  infostring$PGRP         <- syminfo[prima+ 9]
  infostring$HKLASU       <- syminfo[prima+10]
  infostring$MAPASU_CCP4  <- syminfo[prima+11]
  infostring$MAPASU_ZERO  <- syminfo[prima+12]
  infostring$MAPASU_NONZ  <- syminfo[prima+13]
  infostring$CHESHIRE     <- syminfo[prima+14]
  infostring$SYMOP        <- syminfo[prima:dopo][grep("symop",syminfo[prima:dopo])]
  infostring$CENOP        <- syminfo[prima:dopo][grep("cenop",syminfo[prima:dopo])]

  return(infostring)
}


#' Operators of a specific space group
#'
#' Returns human-readable symmetry operators corresponding to a specific input space group.
#'
#' A crystallographic space group includes a set of symmetry operators that can be expressed
#' like operations on the (x,y,z) fractional coordinates of atoms in a unit cell. So, for example,
#' The only operator associated with the space group P 1 is "x,y,z", while the four operators
#' associated with P 21 21 21 are "symop x,y,z", "symop -x+1/2,-y,z+1/2", "symop x+1/2,-y+1/2,-z",
#' "symop -x,y+1/2,-z+1/2".
#'
#' @param SG A character string. The extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
#' @return op_xyz_list A named list made of two vectors. The first vector, SYMOP, contains strings
#'   describing the symmetry operators. The second vector, CENOP, contains strings describing the
#'   centring of the unit cell.
#'
#' @examples
#' # Symmetry operators for space group number 3, P 1 2 1
#' SG <- "P 1 2 1"
#' ltmp <- syminfo_to_op_xyz_list(SG)
#' ltmp
#'
#' @export
syminfo_to_op_xyz_list <- function(SG) {
  # Extract full symmetry information
  data <- extract_symmetry_info(SG)

  # Return NULL if SG not a valid symbol
  if (is.null(data)) {
    return(NULL)
  }

  # Extract "symop" bit
  tmp2 <- strsplit(data$SYMOP," ")
  symop_xyz <- c()
  for (i in 1:length(tmp2))
  {
    symop_xyz <- c(symop_xyz,tmp2[[i]][2])
  }

  # Extract "cenop" bit
  tmp2 <- strsplit(data$CENOP," ")
  cenop_xyz  <- c()
  for (i in 1:length(tmp2))
  {
    cenop_xyz <- c(cenop_xyz,tmp2[[i]][2])
  }
  op_xyz_list <- list(SYMOP=symop_xyz,CENOP=cenop_xyz)

  return(op_xyz_list)
}


#' Human-readable symmetry operator into matrix and vector
#'
#' Returns a \eqn{3\times 3} matrix and \eqn{3\times 1} vector corresponding to either a
#' symmetry operator or a centering operator in human-readable, string form.
#'
#' A string describing a symmetry or a centering operation has a format similar to, for instance,
#' '-x+1/2,-y,z+1/2'. Such a string corresponds to the symmetry operation represented
#' mathematically by the following matrix and vector:
#' \deqn{
#'   \left(\begin{array}{rrr}
#'   -1 & 0 & 0 \\
#'   0 & -1 & 0 \\
#'   0 & 0 & 1
#'   \end{array}\right)\quad,\quad
#'   \left(\begin{array}{r}
#'   1/2 \\
#'   0 \\
#'   1/2
#'   \end{array}\right)
#' }
#' Where symmetry operations in human-readable form are useful for the subjective reasoning
#' in crystallography, their mathematical counterpart is needed for all practical calculations.
#'
#' @param op_xyz A symmetry or centering operation in the form of a human-readable string, e.g.
#'  '-x+1/2,-y,z+1/2'.
#' @return mat_ops A named list including a \eqn{3\times 3} matrix 'R' and a \eqn{3\times 1}
#' vector 'T'.
#'
#' @examples
#' # Reflection and translation
#' sop <- '-x,y+1/2,z'
#' mat_ops <- op_xyz_to_matrix(sop)
#' print(mat_ops)
#'
#' @export
op_xyz_to_matrix <- function(op_xyz)
{
  # Reads in a symmetry or centering operator in character form and output it in
  # matrix or vector form. If input is a symmetry operator, output is a list of
  # a 3X3 matrix and a 3X1 vector; if input is a centering operator, output is
  # still a matrix and a vector, but the matrix is always the identity matrix.

  ltmp <- strsplit(op_xyz,",")

  # If the split is not exactly of three fields return NULL
  if (length(ltmp[[1]]) != 3) {
    warning("Input is not a valid symmetry string!")
    return(NULL)
  }

  # If the split does not contain x, y, or z, return NULL
  flg <- FALSE
  for (a in ltmp[[1]]) {
    if (!grepl("x",a) & !grepl("y",a) & !grepl("z",a)) flg <- TRUE
  }
  if (flg) {
    warning("Input is not a valid symmetry string!")
    return(NULL)
  }

  if (substr(ltmp[[1]][1],1,1) != "-" & substr(ltmp[[1]][1],1,1) != "+")
  {
    stmp <- paste("+",ltmp[[1]][1],sep="")
    ltmp[[1]][1] <- stmp
  }
  if (substr(ltmp[[1]][2],1,1) != "-" & substr(ltmp[[1]][2],1,1) != "+")
  {
    stmp <- paste("+",ltmp[[1]][2],sep="")
    ltmp[[1]][2] <- stmp
  }
  if (substr(ltmp[[1]][3],1,1) != "-" & substr(ltmp[[1]][3],1,1) != "+")
  {
    stmp <- paste("+",ltmp[[1]][3],sep="")
    ltmp[[1]][3] <- stmp
  }

  # Get 3X3 point group matrix
  m <- matrix(rep(0,times=9),nrow=3,ncol=3)
  for (j in 1:length(ltmp[[1]]))
  {
    if (grepl("\\+x",ltmp[[1]][j])) m[j,1] <-  1
    if (grepl("\\-x",ltmp[[1]][j])) m[j,1] <- -1
    if (grepl("\\+y",ltmp[[1]][j])) m[j,2] <-  1
    if (grepl("\\-y",ltmp[[1]][j])) m[j,2] <- -1
    if (grepl("\\+z",ltmp[[1]][j])) m[j,3] <-  1
    if (grepl("\\-z",ltmp[[1]][j])) m[j,3] <- -1
  }

  v <- rep(0,times=3)
  for (j in 1:length(ltmp[[1]]))
  {
    bb <- strsplit(ltmp[[1]][j],"/")
    if (length(bb[[1]]) > 1)
    {
      v[j] <- as.integer(substr(bb[[1]][1],nchar(bb[[1]][1])-1,nchar(bb[[1]][1])))/
        as.integer(bb[[1]][2])
    }
    if (length(bb[[1]]) <= 1)
    {
      v[j] <- 0
    }

    # Turn v into 3X1 vector
    v <- matrix(v,ncol=1)
  }

  # Final list
  mat_ops <- list(R=m,T=v)

  return(mat_ops)
}


#' List of matrices and vectors of a specific space group
#'
#' Returns \eqn{3\times 3} matrices and \eqn{3\times 1} vectors corresponding to point group
#' operations, group translations and cell centring of a given space group.
#'
#' A crystallographic space group consists of a series of transformations on a point
#' \eqn{(x_f,y_f,z_f)} in space that are mathematically implemented as the product of
#' a \eqn{3\times 3} point-group matrix and the point fractional coordinates, \eqn{(x_f,y_f,z_f)},
#' followed by a sum with a \eqn{3\times 1} translation vector. The complete set of points thus
#' produced can be cloned into a new and shifted set translated of an amount represented by a
#' \eqn{3\times 1} centring vector.
#'
#' @param op_xyz_list A named list made of two vectors. The first vector, SYMOP, contains strings
#'   describing the symmetry operators. The second vector, CENOP, contains strings describing the
#'   centring of the unit cell.
#' @return mat_ops_list A named list consisting of 3 lists. The first list, PG, contains
#' \eqn{3\times 3} point group matrices; the second list, T, contains the same number of
#' \eqn{3\times 1} translation vectors. The first matrix is always the identity matrix, the first
#' translation vector is always the null vector. The third list, C, consists of centering vectors;
#' the first centering vector is always the null vector. To summarize, the output looks like the
#' following:
#'
#'            [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
#' where:
#'            I                = identity    3X3 matrix
#'            0                = null        3X1 vector
#'            M2,M3,...,Mn     = point group 3X3 matrices
#'            V2,V3,...,Cn     = translation 3X1 vectors
#'            C2,C3,...,Cm     = centering   3X1 vectors
#'
#' @examples
#' # Symmetry operators for space group number 3, P 1 2 1
#' SG <- "P 1 2 1"
#' op_xyz_list <- syminfo_to_op_xyz_list(SG)
#' mat_ops_list <- op_xyz_list_to_matrix_list(op_xyz_list)
#' names(mat_ops_list)
#'
#' @export
op_xyz_list_to_matrix_list <- function(op_xyz_list)
{
  # Input: "syminfo_to_op_xyz_list" output, i.e. a list of 2 character vectors.
  # The first one contains symmetry operators in x,y,z format; the second one centering operators
  # in x,y,z format.
  # Returns a list consisting of 3 lists. The first list contains 3X3 point group matrices;
  # the second list contains the same number of 3X1 translation vectors. First matrix is always the
  # identity matrix, the first translation vector is always the null vector. the third
  # list consists of centering vectors; the first centering vector is always the null
  # vector. To summarize, the output looks like the following:
  # [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
  # where:
  # I            = identity    3X3 matrix
  # 0            = null        3X1 vector
  # M2,M3,...,Mn = point group 3X3 matrices
  # V2,V3,...,Cn = translation 3X1 vectors
  # C2,C3,...,Cm = centering   3X1 vectors

  # Stop straight away if input is NULL
  if (is.null(op_xyz_list)) {
    return(NULL)
  }

  # Create empty lists
  matrix_list <- list()
  vector_list <- list()
  centering_list <- list()

  # Point group matrices and translation vectors
  for (i in 1:length(op_xyz_list[[1]]))
  {
    data <- op_xyz_to_matrix(op_xyz_list[[1]][i])
    matrix_list[i] <- list(data[[1]])
    vector_list[i] <- list(data[[2]])
  }

  # Centering vectors
  for (i in 1:length(op_xyz_list[[2]]))
  {
    data <- op_xyz_to_matrix(op_xyz_list[[2]][i])
    centering_list[i] <- list(data[[2]])
  }

  mat_ops_list <- list(PG=matrix_list,T=vector_list,C=centering_list)

  return(mat_ops_list)
}



#' Operators of a specific space group in matrix form
#'
#' Returns \eqn{3\times 3} matrices and \eqn{3\times 1} vectors corresponding to point group
#' operations, group translations and cell centring of a given space group.
#'
#' A crystallographic space group consists of a series of transformations on a point
#' \eqn{(x_f,y_f,z_f)} in space that are mathematically implemented as the product of
#' a \eqn{3\times 3} point-group matrix and the point fractional coordinates, \eqn{(x_f,y_f,z_f)},
#' followed by a sum with a \eqn{3\times 1} translation vector. The complete set of points thus
#' produced can be cloned into a new and shifted set translated of an amount represented by a
#' \eqn{3\times 1} centring vector.
#'
#' @param SG A character string. The extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
#' @return mat_ops_list A named list consisting of 3 lists. The first list, PG, contains
#' \eqn{3\times 3} point group matrices; the second list, T, contains the same number of
#' \eqn{3\times 1} translation vectors. The first matrix is always the identity matrix, the first
#' translation vector is always the null vector. The third list, C, consists of centering vectors;
#' the first centering vector is always the null vector. To summarize, the output looks like the
#' following:
#'
#'            [[ [[I,M2,M3,...,Mn]] , [[O,V2,V3,...,Vn]] , [[O,C2,C3,...,Cm]] ]]
#' where:
#'            I                = identity    3X3 matrix
#'            0                = null        3X1 vector
#'            M2,M3,...,Mn     = point group 3X3 matrices
#'            V2,V3,...,Cn     = translation 3X1 vectors
#'            C2,C3,...,Cm     = centering   3X1 vectors
#'
#' @examples
#' # Symmetry operators for space group number 4, P 1 21 1
#' SG <- "P 1 21 1"
#' mat_ops <- syminfo_to_matrix_list(SG)
#' print(mat_ops)
#'
#' @export
syminfo_to_matrix_list <- function(SG)
{
  #
  # This function is simply a wrapper for op_xyz_list_to_matrix_list.

  # 3 functions called in one line (cool, isn't it?)
  return(op_xyz_list_to_matrix_list(syminfo_to_op_xyz_list(SG)))
}


#' Symmetry operations in human readable form
#'
#' This function returns the full set of symmetry operations in human-readable form,
#' each one as a character string starting with 'SYMM'. These are the common crystallographic
#' symmetry operations.
#'
#' @param SG A character string. The extended Hermann-Maguin symbol (e.g. 'P 1 1 21')
#' @return Symm_string A character vector whose components are strings starting by 'SYMM'
#'  and containing the symmetry operations of the given group in human-readable form.
#'
#' @examples
#' # P1 has only one symmetry operation
#' SG <- "P 1"
#' symm_string <- full_symm_strings(SG)
#' print(symm_string)
#'
#' # P 21 21 21 is has many more operations
#' SG <- "P 21 21 21"
#' symm_string <- full_symm_strings(SG)
#' print(symm_string)
#'
#' @export
full_symm_strings <- function(SG)
{
  # Returns all symmetry operators in string format (the kind SYMM)

  # Extract all symmetry information
  sinfo <- extract_symmetry_info(SG)

  if (is.null(sinfo)) {
    return(NULL)
  }

  # 1. Point symmetry
  pgstring <- sinfo$SYMOP
  npgops <- length(pgstring)

  # 2. Centring
  cstring <- sinfo$CENOP
  ncops <- length(cstring)

  # Produce npgops X ncops vectors
  symm_string <- c()

  # Double loop
  for (ssym in pgstring)
  {
    # 3 character vectors for non-centring part
    pgvec <- c("","","")
    i <- 1
    while (substr(ssym,i,i) != ",") i <- i+1
    pgvec[1]<- paste(pgvec[1],substr(ssym,7,i-1),sep="")
    ssym <- substr(ssym,i+1,nchar(ssym))
    i <- 1
    while (substr(ssym,i,i) != ",") i <- i+1
    pgvec[2] <- paste(pgvec[2],substr(ssym,1,i-1),sep="")
    ssym <- substr(ssym,i+1,nchar(ssym))
    pgvec[3] <- paste(pgvec[3],ssym,sep="")
    for (iscen in 1:length(cstring))
    {
      cnvec <- c("","","")
      scen <- cstring[iscen]
      if (iscen == 1)
      {
        symm_string <- c(symm_string,paste("symm ",pgvec[1],",  ",pgvec[2],",  ",pgvec[3],sep=""))
      }
      if (iscen > 1)
      {
        i <- 1
        while (substr(scen,i,i) != ",") i <- i+1
        sa <- substr(scen,7,i-1)
        scen <- substr(scen,i+1,nchar(scen))
        i <- 1
        while (substr(scen,i,i) != ",") i <- i+1
        sb <- substr(scen,1,i-1)
        scen <- substr(scen,i+1,nchar(scen))
        sc <- scen
        cnvec <- c("","","")
        if (nchar(sa) > 2 & substr(sa,2,2) == "+" | substr(sa,2,2) == "-") cnvec[1] <- substr(sa,2,nchar(sa))
        if (nchar(sa) > 2 & substr(sa,3,3) == "+" | substr(sa,3,3) == "-") cnvec[1] <- substr(sa,3,nchar(sa))
        if (nchar(sb) > 2 & substr(sb,2,2) == "+" | substr(sb,2,2) == "-") cnvec[2] <- substr(sb,2,nchar(sb))
        if (nchar(sb) > 2 & substr(sb,3,3) == "+" | substr(sb,3,3) == "-") cnvec[2] <- substr(sb,3,nchar(sb))
        if (nchar(sc) > 2 & substr(sc,2,2) == "+" | substr(sc,2,2) == "-") cnvec[3] <- substr(sc,2,nchar(sc))
        if (nchar(sc) > 2 & substr(sc,3,3) == "+" | substr(sc,3,3) == "-") cnvec[3] <- substr(sc,3,nchar(sc))
        symm_string <- c(symm_string,paste("symm ",pgvec[1],cnvec[1],",  ",pgvec[2],cnvec[2],",  ",pgvec[3],cnvec[3],sep=""))
      }
    }
  }
  symm_string <- toupper(symm_string)
  for (i in 1:length(symm_string)) while(nchar(symm_string[i]) != 80) symm_string[i] <- paste(symm_string[i]," ",sep="")

  return(symm_string)
}


#' Correct spelling for Herman-Maguin space groups symbols
#'
#' The commonly-used spelling of a crystallographic space group does not match the correct definition
#' given by the Herman-Maguin symbols which define all space groups in a unique and precise way. This
#' function attempt to translate a tentative string into a possible Herman-Maguin symbol, if it finds
#' one. If the input string is already in the extended Herman-Maguin form, the same string is returned
#' as output.
#'
#' @param sym_xHM A character string. The space group symbol in its commonly-used spelling.
#' @return SG A character string. The extended Hermann-Maguin symbol (e.g. 'P 1 1 21').
#'
#' @examples
#' # P21
#' print(find("P 21"))
#'
#' # P -1
#' print(find("P-1"))
#'
#' @export
findHM <- function(sym_xHM)
{
  # Eliminate trailing white spaces, if any
  sym_xHM <- trimws(sym_xHM,"both")

  # Create a space after first character, if it does not exist
  if (substr(sym_xHM,2,2) != " ") {
    tmp <- substr(sym_xHM,2,nchar(sym_xHM))
    sym_xHM <- paste0(substr(sym_xHM,1,1)," ",tmp)
  }

  # Try with translate_SG first
  tmp <- translate_SG(sym_xHM,"Hall","xHM")
  if (tmp$ans) {
    SG <- tmp$msg
    return(SG)
  }
  tmp <- translate_SG(sym_xHM,"old","xHM")
  if (tmp$ans) {
    SG <- tmp$msg
    return(SG)
  }
  # For an xHM input is trickier
  tmp <- translate_SG(sym_xHM,"xHM","number")
  if (tmp$ans) {
    sym_num <- tmp$msg
    for (i in 1:18) {
      tmp2 <- translate_SG(sym_num,set=i)
      if (tmp2$ans) {
        if (sym_xHM == tmp2$msg) return(sym_xHM)
      }
    }
  }

  #tmp <- translate_SG(sym_xHM,"xHM","xHM")
  #if (tmp$ans) {
  #  SG <- tmp$msg
  #  return(SG)
  #}

  # Possible inputs
  if (sym_xHM == "P 2") {
    sym_xHM <- "P 1 2 1"
  } else if (sym_xHM == "P 21") {
    sym_xHM <- "P 1 21 1"
  } else if (sym_xHM == "C 2") {
    sym_xHM <- "C 1 2 1"
  } else if (sym_xHM == "P m") {
    sym_xHM <- "P 1 m 1"
  } else if (sym_xHM == "P c") {
    sym_xHM <- "P 1 c 1"
  } else if (sym_xHM == "C m") {
    sym_xHM <- "C 1 m 1"
  } else if (sym_xHM == "C c") {
    sym_xHM <- "C 1 c 1"
  } else if (sym_xHM == "P 2/m") {
    sym_xHM <- "P 1 2/m 1"
  } else if (sym_xHM == "P 21/m") {
    sym_xHM <- "P 1 21/m 1"
  } else if (sym_xHM == "C 2/m") {
    sym_xHM <- "C 1 2/m 1"
  } else if (sym_xHM == "P 2/c") {
    sym_xHM <- "P 1 2/c 1"
  } else if (sym_xHM == "P 21/c") {
    sym_xHM <- "P 1 21/c 1"
  } else if (sym_xHM == "C 2/c") {
    sym_xHM <- "C 1 2/c 1"
  } else if (sym_xHM == "P n n n") {
    sym_xHM <- "P n n n :1"
  } else if (sym_xHM == "P b a n") {
    sym_xHM <- "P b a n :1"
  } else if (sym_xHM == "P m m n") {
    sym_xHM <- "P m m n :1"
  } else if (sym_xHM == "C c c a") {
    sym_xHM <- "C c c a :1"
  } else if (sym_xHM == "F d d d") {
    sym_xHM <- "F d d d :1"
  } else if (sym_xHM == "P 4/n") {
    sym_xHM <- "P 4/n :1"
  } else if (sym_xHM == "P 42/n") {
    sym_xHM <- "P 42/n :1"
  } else if (sym_xHM == "P 41/a") {
    sym_xHM <- "P 41/a :1"
  } else if (sym_xHM == "P 4/n b m") {
    sym_xHM <- "P 4/n b m :1"
  } else if (sym_xHM == "P 4/n n c") {
    sym_xHM <- "P 4/n n c :1"
  } else if (sym_xHM == "P 4/n m m") {
    sym_xHM <- "P 4/n m m :1"
  } else if (sym_xHM == "P 4/n c c") {
    sym_xHM <- "P 4/n c c :1"
  } else if (sym_xHM == "P 42/n b c") {
    sym_xHM <- "P 42/n b c :1"
  } else if (sym_xHM == "P 42/n n m") {
    sym_xHM <- "P 42/n n m :1"
  } else if (sym_xHM == "P 42/n m c") {
    sym_xHM <- "P 42/n m c :1"
  } else if (sym_xHM == "P 42/n c m") {
    sym_xHM <- "P 42/n c m :1"
  } else if (sym_xHM == "I 41/a m d") {
    sym_xHM <- "I 41/a m d :1"
  } else if (sym_xHM == "I 41/a c d") {
    sym_xHM <- "I 41/a c d :1"
  } else if (sym_xHM == "R 3") {
    sym_xHM <- "R 3 :H"
  } else if (sym_xHM == "H 3") {
    sym_xHM <- "R 3 :H"
  } else if (sym_xHM == "R -3") {
    sym_xHM <- "R -3 :H"
  } else if (sym_xHM == "R 3 2") {
    sym_xHM <- "R 3 2 :H"
  } else if (sym_xHM == "R 3 m") {
    sym_xHM <- "R 3 m :H"
  } else if (sym_xHM == "R 3 c") {
    sym_xHM <- "R 3 c :H"
  } else if (sym_xHM == "R -3 m") {
    sym_xHM <- "R -3 m :H"
  } else if (sym_xHM == "R -3 c") {
    sym_xHM <- "R -3 c :H"
  } else if (sym_xHM == "P n -3") {
    sym_xHM <- "P n -3 :1"
  } else if (sym_xHM == "F d -3") {
    sym_xHM <- "F d -3 :1"
  } else if (sym_xHM == "P n -3 n") {
    sym_xHM <- "P n -3 n :1"
  } else if (sym_xHM == "F d -3 m") {
    sym_xHM <- "F d -3 m :1"
  } else if (sym_xHM == "F d -3 c") {
    sym_xHM <- "F d -3 c :1"
  } else {
    sym_xHM <- NULL
  }
  SG <- sym_xHM

  # Returned, corrected xHM symbol
  return(SG)
}


#' Number of space group settings
#'
#' Although a space group is uniquely defined, i.e. the symmetry operations defining it
#' is uniquely given, the choice of vectors that defines a unit cell for that symmetry
#' is not unique. The different choices are known as settings. Most space groups have only
#' one setting, but it is possible to find space groups with several settings. For example,
#' "C 1 2/c 1" has 18 settings. While the xHM symbol for setting 1 is "C 1 2/c 1", the
#' symbol for setting 2 is "A 1 2/n 1", etc.
#'
#' @param SG A character string or a number indicating the space group.
#' @return nsett The number of setting for the given space group.
#'
#' @examples
#' # P 1 21 1 (group number 4) has three settings
#' num_symm_settings(4)
#'
#' # Find the extended Hermann-Mauguin symbols
#' translate_SG(4,"number","xHM",1)$msg
#' translate_SG(4,"number","xHM",2)$msg
#' translate_SG(4,"number","xHM",3)$msg
#'
#' @export
num_symm_settings <- function(SG=NULL) {
  # Stop if no input given
  if (is.null(SG)) {
    stop("Space group input needed.")
  }

  # Whatever the input, turn it into a space group number
  flg <- FALSE
  tmp <- translate_SG(SG,"number","number")
  if (tmp$ans) {
    gn <- tmp$msg
    flg <- TRUE
  }
  if (!flg) {
    tmp <- translate_SG(SG,"ccp4","number")
    if (tmp$ans) {
      gn <- tmp$msg
      flg <- TRUE
    }
  }
  if (!flg) {
    tmp <- translate_SG(SG,"xHM","number")
    if (tmp$ans) {
      gn <- tmp$msg
      flg <- TRUE
    }
  }
  if (!flg) {
    tmp <- translate_SG(SG,"Hall","number")
    if (tmp$ans) {
      gn <- tmp$msg
      flg <- TRUE
    }
  }
  if (!flg) {
    tmp <- translate_SG(SG,"old","number")
    if (tmp$ans) {
      gn <- tmp$msg
      flg <- TRUE
    }
  }

  # A last attempt, if previous have failed
  if (!flg) {
    tmp <- findHM(SG)
    if (!is.null(tmp)) {
      gn <- translate_SG(tmp,"xHM","number")
      flg <- TRUE
    }
  }

  # Give up!
  if (!flg) {
    stop("Input does not correspond to any valid space group.")
  }

  # Now cycle until no more settings for that specific space group
  # number are found
  nsett <- 0
  lista <- list()
  lista$ans <- TRUE
  lista$msg <- NULL
  while (lista$ans)
  {
    nsett <- nsett+1
    lista <- translate_SG(value=gn,SG_in="number",SG_out="xHM",
                          set=nsett)
  }
  nsett <- nsett-1

  return(nsett)
}


#' Find specific space group setting
#'
#' Although a space group is uniquely defined, i.e. the symmetry
#' operations defining it are uniquely given, the choice of
#' vectors that defines a unit cell for that symmetry is not
#' unique. The different choices are known as settings. Most
#' space groups have only one setting, but it is possible to find
#' space groups with several settings. For example, "C 1 2/c 1"
#' has 18 settings. The xHM symbol for setting 1 is
#' "C 1 2/c 1", the symbol for setting 2 is "A 1 2/n 1", etc.
#'
#' @param SG A character string indicating the extended
#'           Hermann-Mauguin symbol for the space group.
#' @return set An integer equal to the specific setting
#'             corresponding to the given xHM symbol.
#'
#' @examples
#' # P2 (group number 4) has three settings
#' nsets <- find_symm_setting("P 1 2 1")
#' print(nsets)
#'
#' @export
find_symm_setting <- function(SG) {
  # Is input a character?
  ans <- is.character(SG)
  if (!ans) {
    msg <- paste("Input needs to be an extended",
                 "Hermann-Mauguin symbol.\n")
    cat(msg)

    return(NULL)
  }

  # Verify it's a valid symbol
  ans <- findHM(SG)
  if (!is.null(ans)) SG <- ans
  tmp <- translate_SG(SG,SG_in="xHM",SG_out="number")
  if (!tmp$ans) {
    msg <- "No space group corresponding to the given symbol.\n"
    cat(msg)

    return(NULL)
  }
  SGnum <- tmp$msg

  # Now check how many settings
  nsets <- num_symm_settings(SG)

  # Loop until match is found
  setN <- 0
  for (i in 1:nsets) {
    tmp <- translate_SG(SGnum,SG_in="number",SG_out="xHM",set=i)
    ans <- SG == tmp$msg
    if (ans) setN <- i
  }

  # Safety check (if no setting found)
  if (setN == 0) {
    msg <- "No setting corresponding to given xHM symbol.\n"
    cat(msg)

    return(NULL)
  }

  return(setN)
}

#' Space group compatible with given cell
#'
#' This function returns the space group, among those compatible
#' with the given unit cell, with the lowest symmetry group
#' number.
#'
#' A given unit cell is compatible with several symmetries.
#' For example, a cell with different sides and different
#' angles, also different from 90 degrees, is compatible with
#' the triclinic lattice system, which corresponds to space
#' groups P1 and P -1. A cell with different sides, alpha = 90,
#' gamma=90 and beta different from 90, is compatible with the
#' monoclinic lattice system, which corresponds to space groups
#' from number 3 ("P 1 2 1") to number 15 (C 1 2/c 1"). In the
#' first case this function returns "P 1", while in the second
#' case it returns "P 1 2 1".
#'
#' @param uc An object of class 'unit_cell'.
#' @return csym An object of class 'cryst_symm', corresponding
#'              to the selected, lowest symmetry.
#'
#' @examples
#' # Monoclinic cell
#' uc <- unit_cell(10,20,15,90,110,90)
#'
#' # The selected space group is "P 1 2 1"
#' csym <- lowest_uc_compatible_SG(uc)
#' print(csym)
#'
#' @export
lowest_uc_compatible_SG <- function(uc) {
  # Is input a 'unit_cell' object?
  ans <- check_unit_cell_validity(uc)
  if (!ans) {
    msg <- paste("Input is not a valid object",
                 "of class 'unit_cell'.\n")
    cat(msg)

    return(NULL)
  }

  # Starting group is "P 1"
  SG <- "P 1"

  # Cell parameters
  a <- uc$a
  b <- uc$b
  c <- uc$c
  aa <- uc$alpha[1]
  bb <- uc$beta[1]
  cc <- uc$gamma[1]

  # Determine lattice system
  dab <- abs(a-b)
  dac <- abs(a-c)
  dbc <- abs(b-c)
  daabb <- abs(aa-bb)
  daacc <- abs(aa-cc)
  dbbcc <- abs(bb-cc)
  da90 <- abs(aa-90)
  db90 <- abs(bb-90)
  dc90 <- abs(cc-90)
  dc120 <- abs(cc-120)

  ## MONOCLINIC ##
  # 'P 1 2 1'
  if (dab > 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 < 0.000001 & db90 > 0.000001 & dc90 < 0.000001) {
    SG <- "P 1 2 1"
  }

  # 'P 1 1 2'
  if (dab > 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 < 0.000001 & db90 < 0.000001 & dc90 > 0.000001) {
    SG <- "P 1 1 2"
  }

  # 'P 2 1 1'
  if (dab > 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 > 0.000001 & db90 < 0.000001 & dc90 < 0.000001) {
    SG <- "P 2 1 1"
  }

  ## ORTHOROMBIC ##
  # 'P 2 2 2'
  if (dab > 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 < 0.000001 & db90 < 0.000001 & dc90 < 0.000001) {
    SG <- "P 2 2 2"
  }

  ## TETRAGONAL ##
  # 'P 4'
  if (dab < 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 < 0.000001 & db90 < 0.000001 & dc90 < 0.000001) {
    SG <- "P 4"
  }

  ## HEXAGONAL (TRIGONAL) ##
  # 'P 3'
  if (dab < 0.000001 & dac > 0.000001 & dbc > 0.000001 &
      da90 < 0.000001 & db90 < 0.000001 & dc120 < 0.000001) {
    SG <- "P 3"
  }

  ## ROMBOHEDRAL ##
  # 'R3 :R'
  if (dab < 0.000001 & dac < 0.000001 & dbc < 0.000001 &
      da90 > 0.000001 & db90 > 0.000001 & dc90 > 0.000001 &
      daabb < 0.000001 & daacc < 0.000001 & dbbcc < 0.000001) {
    SG <- "R 3 :R"
  }

  ## CUBIC ##
  # 'P 2 3'
  if (dab < 0.000001 & dac < 0.000001 & dbc < 0.000001 &
      da90 < 0.000001 & db90 < 0.000001 & dc90 < 0.000001) {
    SG <- "P 2 3"
  }

  # Object of class 'cryst_symm'
  csym <- cryst_symm(SG)

  return(csym)
}


#' Cell parameter constrains from symmetry
#'
#' This function returns a set of constrains, as string character
#' expressions, imposed by the specific symmetry group on the
#' given unit cell.
#'
#' Space group symmetry imposes certain constraints on the
#' values that unit cell parameters can take. For example, the
#' symmetry represented by the monoclinic space group of extended
#' Hermann-Mauguin symbol "P 1 2 1" is compatible with a unit cell
#' in which alpha=gamma=90.
#'
#' There is just a handful of constrains for unit cells. Here they
#' are indicated with the following set of specific strings:
#' \itemize{
#'   \item \strong{'No constrains'} Like in a triclinic cell.
#'   \item \strong{'alpha=90'} The alpha angle is fixed at 90
#'                           degrees.
#'   \item \strong{'beta=90'} The beta angle is fixed at 90
#'                          degrees.
#'   \item \strong{'gamma=90'} The gamma angle is fixed at 90
#'                             degrees.
#'   \item \strong{'gamma=120'} The gamma angle is fixed at 120
#'                              degrees.
#'   \item \strong{'alpha=beta=gamma'} The three angle have the
#'                                     same value, different
#'                                     from 90 degrees.
#'   \item \strong{'a=b'} Cell side a is equal to cell side b.
#'   \item \strong{'a=b=c'} The three cell sides are equal.
#' }
#'
#' @param SG A character string indicating the extended
#'           Hermann-Mauguin symbol for the space group.
#' @return vcons A character vector. Each component is a string,
#'               like 'alpha=90' or 'a=b', that describes the
#'               type of constrain to be applied to a unit cell
#'               of a crystal structure with given space group
#'               symmetry (see above).
#' @examples
#' # P 1 1 2 (group number 3) corresponds to setting 2
#' SG <- translate_SG(3,set=2)
#'
#' # Constrains for this symmetry
#' stmp <- symm_to_cell_const(SG)
#' print(stmp)
#'
#' # R 3 (rombohedral setting)
#' stmp <- symm_to_cell_const("R 3 :R")
#' print(stmp)
#'
#' @export
symm_to_cell_const <- function(SG) {
  ans <- is.character(SG)
  if (!ans) {
    msg <- paste("Input needs to be an extended",
                 "Hermann-Mauguin symbol.\n")
    cat(msg)

    return(NULL)
  }
  ans <- findHM(SG)
  if (!is.null(ans)) SG <- ans
  tmp <- translate_SG(SG,SG_in="xHM",SG_out="number")
  if (!tmp$ans) {
    msg <- "No space group corresponding to the given symbol.\n"
    cat(msg)

    return(NULL)
  }

  # Space group number
  sym_number <- tmp$msg

  # Find correct setting
  setting <- find_symm_setting(SG)

  # Long list for all space groups and settings
  vcons <- "No constrains"

  #########################################################################################################################
  #########
  #########            TRICLINIC: 1 to 2
  #########
  #########################################################################################################################
  # Nothing to do for sym_number 1
  # Nothing to do for sym_number 2

  #########################################################################################################################
  #########
  #########            MONOCLINIC: 3 to 15
  #########
  #########################################################################################################################
  #
  # Generic constrains for system
  if (sym_number >= 3 & sym_number <= 15) {
    vcons <- c("alpha=90","gamma=90")
  }

  if (sym_number == 3 | sym_number == 4 | sym_number == 6 |
      sym_number == 10 | sym_number == 11) {
    if (setting == 2) vcons <- c("alpha=90","beta=90")
    if (setting == 3) vcons <- c("beta=90","gamma=90")

    return(vcons)
  }
  if (sym_number == 5 | sym_number == 12 | sym_number == 13 |
      sym_number == 14) {
    if (setting == 4 | setting == 5 | setting == 6)
      vcons <- c("alpha=90","beta=90")
    if (setting == 7 | setting == 8 | setting == 9)
      vcons <- c("beta=90","gamma=90")

    return(vcons)
  }
  if (sym_number == 9 | sym_number == 15) {
    if (setting == 7 | setting == 8 | setting == 9 |
        setting == 10 | setting == 11 | setting == 12)
      vcons <- c("alpha=90","beta=90")
    if (setting == 13 | setting == 14 | setting == 15 |
        setting == 16 | setting == 17 | setting == 18)
      vcons <- c("beta=90","gamma=90")

    return(vcons)
  }

  #########################################################################################################################
  #########
  #########            ORTHOROMBIC: 16 to 74
  #########
  #########################################################################################################################
  #
  # Generic constrains for system
  if (sym_number >= 16 & sym_number <= 74) {
    vcons <- c("alpha=90","beta=90","gamma=90")
  }

  #########################################################################################################################
  #########
  #########            TETRAGONAL: 75 to 142
  #########
  #########################################################################################################################
  # Generic constrains for system
  if (sym_number >= 75 & sym_number <= 142) {
    vcons <- c("a=b","alpha=90","beta=90","gamma=90")
  }

  #########################################################################################################################
  #########
  #########            TRIGONAL: 143 to 167
  #########
  #########################################################################################################################
  # Generic constrains for the system
  if (sym_number >= 143 & sym_number <= 167) {
    vcons <- c("a=b","alpha=90","beta=90","gamma=120")
  }
  if (sym_number == 146 | sym_number == 148 | sym_number == 155 |
      sym_number == 160 | sym_number == 161 | sym_number == 166 |
      sym_number == 167) {
    if (setting == 2) vcons <- c("a=b=c","alpha=beta=gamma")

    return(vcons)
  }

  #########################################################################################################################
  #########
  #########            HEXAGONAL: 168 to 194
  #########
  #########################################################################################################################
  # Generic constrains for the system
  if (sym_number >= 168 & sym_number <= 194) {
    vcons <- c("a=b","alpha=90","beta=90","gamma=120")
  }

  #########################################################################################################################
  #########
  #########            CUBIC: 195 to 230
  #########
  #########################################################################################################################
  # # Generic constrains for the system
  if (sym_number >= 195 & sym_number <= 230) {
    vcons <- c("a=b=c","alpha=90","beta=90","gamma=90")
  }

  # This is returned when setting is 1
  return(vcons)
}



#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
## Auxiliary functions, only used here
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

.translate_SG <- function(value,SG_in="number",SG_out="xHM",set=1) {
  # A few checks
  if (SG_in != "number" & SG_in != "ccp4" & SG_in != "Hall" & SG_in != "xHM" & SG_in != "old")
    stop("Wrong SG_in string. Valid strings are: number ccp4 Hall xHM old")
  if (SG_out != "number" & SG_out != "ccp4" & SG_out != "Hall" & SG_out != "xHM" & SG_out != "old")
    stop("Wrong SG_out string. Valid strings are: number ccp4 Hall xHM old")

  # If number is not in the range 1:230 stop
  if (SG_in == "number")
  {
    idx <- which(value == 1:230)
    if (length(idx) == 0)
    {
      msg <- "There is not a space group associated with the input number."
      return(list(msg=msg,ans=FALSE))
    }
  }

  # If input is "number", turn into a string
  if (SG_in == "number") value <- paste("number ",value)
  if (SG_in == "ccp4") value <- paste("symbol ccp4",value)

  # Complete string for non-number cases
  if (SG_in == "Hall") value <- paste("symbol Hall '",value,"'",sep="")
  if (SG_in == "xHM") value <- paste("symbol xHM  '",value,"'",sep="")
  if (SG_in == "old") value <- paste("symbol old  '",value,"'",sep="")

  # The whole content of "syminfo.lib" is already in object syminfo
  bsg <- grep(value,syminfo,fixed=TRUE)

  # Select correct one in the "number" case
  if (SG_in == "number")
  {
    bsg2 <- c()
    for (i in 1:length(bsg))
    {
      numero1 <- strsplit(syminfo[bsg[i]],"  ")[[1]][2]
      numero2 <- strsplit(value,"  ")[[1]][2]
      if (numero1 == numero2) bsg2 <- c(bsg2,bsg[i])
    }
    bsg <- bsg2
    rm(bsg2)
  }
  #print(paste("There are ",length(bsg)," settings for this space group"))

  # Select case based on set
  if (length(bsg) < set)
  {
    msg <- "Something wrong in your input:"
    msg <- paste(msg,"   1) the symbol or number input for this space group does not exist",sep="\n")
    msg <- paste(msg,
                 "   2) if your inpur was a number, perhaps for this space group there are not that many settings",
                 sep="\n")
    return(list(msg=msg,ans=FALSE))
  }
  bsg <- bsg[set]

  bini <- grep("begin_spacegroup",syminfo)
  prima <- bini[length(bini[bini < bsg])]
  dtmp <- bini[bini > bsg][1]
  if (!is.na(dtmp)) {
    dopo  <- dtmp-1   # Add 1 for those cases when
                      # prima = dopo
  } else {
    dopo <- length(syminfo)
  }

  if (SG_out == "number")
  {
    key <- "number"
  }
  if (SG_out != "number")
  {
    key <- paste("symbol",SG_out)
  }

  tmp <- syminfo[prima:dopo][grep(key,syminfo[prima:dopo])]
  if (key == "number") translated_value <- strsplit(tmp,"  ")[[1]][2]
  if (key != "number")
  {
    tmp2 <- strsplit(tmp," ")[[1]]
    if (tmp2[2] == "ccp4") translated_value <- tmp2[3]
    if (tmp2[2] == "Hall") translated_value <- strsplit(tmp,"'")[[1]][2]
    if (tmp2[2] == "xHM") translated_value <- strsplit(tmp,"'")[[1]][2]
    if (tmp2[2] == "old") translated_value <- strsplit(tmp,"'")[[1]][2]
    rm(tmp2)
  }
  rm(tmp)

  # If output requires number, turn character into numeric
  if (SG_out == "number" | SG_out == "ccp4") translated_value <- as.integer(translated_value)

  list_SG <- list(msg=translated_value,ans=TRUE)

  # "Funny" cases at the end of CCP4 syminfo file
  if (list_SG$ans & (substr(translated_value,1,2) == "NA")) {
    list_SG$ans <- FALSE
  }

  return(list_SG)
}

