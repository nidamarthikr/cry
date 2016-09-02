#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads and output an MTZ header
#'
#' @param filename A character string. The path to a valid mtz file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the mtz header.
#' @return A named list. Each name correspond to a valid field in the mtz
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/mtz/file.mtz"
#' ltmp <- .readMTZHeader(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export
readMTZHeader <- function(filename,messages=TRUE)
{
  # Reads and output MTZ header (version running independently of readMTZ)

  # Create a connection to binary file
  f <- file(filename, open="rb")

  # Reads initial record
  irdata <- .readIR(f)
  errF <- irdata[[1]]
  if (errF != 0) stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

  # Work out how many reflection records are contained in this MTZ file (= nref X ncols)
  numDataItems <- irdata[[2]] - 20 - 1

  # Load all reflection records (only use here is for getting to the right point of binary file)
  if (irdata[[3]] == 1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
  if (irdata[[3]] == -1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")

  # Load header information
  hData <- .readH(f,messages)

  # Close connection
  close(f)

  # Return header stuff as list
  return(hData[[1]])
}


#' Reads and output the full content of an MTZ file
#'
#' Reads mtz files and store both header information and reflection data
#' records in named lists. A third list is used, if the mtz file is an
#' unmerged file, for storing batch headers.
#'
#' @param filename A character string. The path to a valid mtz file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting data included in the mtz file.
#' @return A named list of length 3. The first element is called "reflections"
#'    and is a dataframe with as many columns as are included in the mtz file.
#'    The name of each column of the dataframe coincides with the name of the
#'    corresponding column in the mtz. The second element is called "header"
#'    and is a named list in which each name correspond to a valid field in
#'    the mtz header. The third element is called "batch_header" and is a
#'    list with as many elements as the number of batches (images)
#'    included in the mtz file. Each list element is, itself, a named list
#'    including all the useful variables stored in batch headers. If no batch
#'    headers are contained in the file (merged mtz), the batch_header
#'    element is NULL.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/mtz/file.mtz"
#' ltmp <- .readMTZ(filename)
#' print(names(ltmp))
#' print(class(ltmp$reflections))
#' str(ltmp$reflections)
#' print(class(ltmp$header))
#' print(class(ltmp$batch_header))
#'
#' refs <- ltmp$reflections
#' print(colnames(refs))
#' print(range(ltmp$H))
#' }
#' @export
readMTZ <- function(filename, messages = TRUE){

  ################################################################################
  ## a function for reading an MTZ file                                         ##
  ## N.B. assumes 4 bytes per data item, not sure if this is always true        ##
  ##                                                                            ##
  ## David Waterman and James Foadi. Diamond Light Source and Imperial College. ##
  ## Started off by David Waterman in June 2009.                                ##
  ## Extended by James Foadi in July 2009                                       ##
  ################################################################################

  # Create a connection to binary file
  f <- file(filename, open="rb")

  # Reads initial record
  irdata <- .readIR(f)
  errF <- irdata[[1]]
  if (errF != 0) stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

  # Work out how many reflection records are contained in this MTZ file (= nref X ncols)
  numDataItems <- irdata[[2]] - 20 - 1

  # Load all reflection records
  if (irdata[[3]] == 1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
  if (irdata[[3]] == -1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")

  # Load header information
  hData <- .readH(f,messages)
  hF <- hData[[2]]                                    # hF == 1 means no batch headers
  # Turn reflnData into a dataframe where columns are named
  tmpmatrix <- matrix(data=reflnData,nrow=hData[[1]]$NCOL[2],ncol=hData[[1]]$NCOL[1],byrow=TRUE,dimnames=list(NULL,hData[[1]]$COLUMN$labels))
  tmpdataframe <- as.data.frame(tmpmatrix)
  reflnData <- tmpdataframe
  rm(tmpmatrix,tmpdataframe)

  # If no batch headers are contained, program ends here
  if (hF == 1)
  {
    close(f)                  # Close connection
    if (hData[[1]]$NCOL[3] > 0) warning("MTZ appears to be a multi-record file, but it does not include batch headers.")

    # Output data are packed in a list
    data <- list(reflections=reflnData,header=hData[[1]],batch_header=NULL)

    # Create "mtz" class
    #class(data) <- "mtz"

    return(data)
  }

  # Batch headers loop
  bhData <- .readBH(f)

  # Close connection
  close(f)

  # Output data are packed in a list
  data <- list(reflections=reflnData,header=hData[[1]],batch_header=bhData)

  # Create "mtz" class
  #class(data) <- "mtz"

  return(data)
}
