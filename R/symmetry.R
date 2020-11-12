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
#' @param SG_in A string representing the space group format. Possible values are:
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
  dopo  <- bini[bini > bsg][1]-1   # Add 1 for those cases when prima = dopo
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
