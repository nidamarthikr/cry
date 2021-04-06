#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads and output an MTZ header
#'
#' An MTZ file is a binary file created to carry information on
#' x-ray diffraction experiments on crystals. It includes x-ray
#' diffraction data and information on the experiment and the
#' crystal.
#'
#' The function returns a named list whose components are
#' the \code{reflections}, the \code{header} and the
#' \code{batch_header}. The \code{header} is a named list whose
#' components are:
#' \describe{
#'   \item{TITLE}{A character string containing the title of the
#'                MTZ file.}
#'   \item{NCOL}{Number of columns in data frame
#'    \code{reflections}.}
#'   \item{CELL}{A numeric vector of length 6, containing the
#'    unit cell parameters.}
#'   \item{SORT}{An integer vector of length 5, containing the
#'    sort order of the first 5 columns of data.}
#'   \item{SYMINF}{Un-named list with 6 components: the number
#'    of symmetry operations (an integer), the number of
#'    primitive operations (an integer), the lattice type
#'    (a one-letter character), the space group number (an
#'    integer), the space group name (a 10-letter character
#'    string) and the point group name (a 6-letter character).}
#'   \item{RESO}{Minimum and maximum data resolution, stored as
#'    \eqn{{1/d^2}}.}
#'   \item{NDIF}{Number of datasets whose reflection data are
#'    present in the file.}
#'   \item{SYMM}{A character vector whose length depends on the
#'    type of symmetry. It describes the symmetry operations in
#'    human-readable format, International Tables style. Each
#'    string is 80 characters long.}
#'   \item{PROJECT}{A data frame whose rows provide an ID and a
#'    name (called "pname") for the projects for which the data
#'    contained have been produced.}
#'   \item{CRYSTAL}{A data frame whose rows provide an ID and a
#'    name (called "pname") for the crystals for which the data
#'    contained have been produced.}
#'   \item{DATASET}{A data frame whose rows provide an ID and a
#'    name (called "pname") for the datasets included in the
#'    reflections record.}
#'   \item{DCELL}{A data frame whose rows contain the
#'    \code{CRYSTAL} IDs and cell parameters of each crystal
#'    that has contributed to the data in the reflections record.}
#'   \item{DWAVEL}{A data frame whose rows contain the
#'    \code{DATASET} IDs and the wavelength with which the
#'     reflection data were collected.}
#'   \item{COLUMN}{A data frame describing the type of data
#'    included in the reflections record. The data frame includes
#'    the labels for each column, the dtype (see
#'    \code{\link{merged_reflections}}) for each column, min and
#'    max values and the \code{DATASET} ID.}
#'   \item{COLSRC}{A data frame with three columns. The first
#'    includes the labels of each reflections record column. The
#'    second includes a time stamp of when each data column was
#'    created. The third is the dataset ID as a string.}
#'   \item{COLGRP}{A character string vector where each component
#'    is an 80-letters string describing the name and type of
#'    data.}
#'   \item{HISTORY}{A character string vector of variable length.
#'    Each component is an 80-letter string summarising the steps
#'    that lead to the curren reflections record. HISTORY can
#'    contain a maximum of 30 lines.}
#'}
#'
#' @param filename A character string. The path to a valid mtz file.
#' @param message A logical variable. If TRUE the function prints
#'    a message highlighting what is included in the mtz header.
#'    Default value is \code{message=FALSE}.
#' @return A named list. Each name correspond to a valid field in the mtz
#'    header (see details).
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei_phases.mtz")
#' ltmp <- readMTZHeader(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#'
#' @export
readMTZHeader <- function(filename,message=FALSE)
{
  # Reads and output MTZ header (version running independently
  # of readMTZ)

  # Create a connection to binary file
  f <- file(filename, open="rb")

  # Reads initial record
  irdata <- .readIR(f)
  errF <- irdata[[1]]
  if (errF != 0)
  stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

  # Work out how many reflection records are contained in this MTZ file (= nref X ncols)
  numDataItems <- irdata[[2]] - 20 - 1

  # Load all reflection records (only use here is
  # for getting to the right point of binary file)
  if (irdata[[3]] == 1)
    reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
  if (irdata[[3]] == -1)
    reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")

  # Load header information
  hData <- .readH(f,message)

  # Close connection
  close(f)

  # Return header stuff as list
  return(hData[[1]])
}


#' Load an MTZ file
#'
#' Reads mtz files and store both header information and
#' reflection data records in named lists. A third list is
#' used, if the mtz file is an unmerged file, for storing batch
#' headers.
#'
#' @param filename A character string. The path to a valid mtz
#'                 file.
#' @param message A logical variable. If TRUE the
#'                function prints a message highlighting data
#'                included in the mtz file. Default value is
#'                \code{message=FALSE}.
#' @return A named list of length 3. The first element is
#'         called "reflections" and is a dataframe with as
#'         many columns as are included in the mtz file.
#'         The name of each column of the dataframe coincides
#'         with the name of the corresponding column in the mtz.
#'         The second element is called "header" and is a named
#'         list in which each name correspond to a valid field
#'         in the mtz header (see details in
#'         \code{\link{readMTZHeader}}).
#'         The third element is called
#'         "batch_header" and is a list with as many elements as
#'         the number of batches (images) included in the mtz
#'         file. Each list element is, itself, a named list
#'         including all the useful variables stored in batch
#'         headers. If no batch headers are contained in the file
#'         (merged mtz), the batch_header element is NULL.
#' @examples
#'
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei_phases.mtz")
#' ltmp <- readMTZ(filename)
#' print(names(ltmp))
#' print(class(ltmp$reflections))
#' str(ltmp$reflections)
#' print(class(ltmp$header))
#' print(class(ltmp$batch_header))
#'
#' refs <- ltmp$reflections
#' print(colnames(refs))
#' print(range(refs$H))
#'
#' @export
readMTZ <- function(filename,message=FALSE) {
  # Create a connection to binary file
  f <- file(filename, open="rb")

  # Reads initial record
  irdata <- .readIR(f)
  errF <- irdata[[1]]
  if (errF != 0) stop("readMTZ: MTZ binary file does not contain initial \"MTZ \" tag. It is either a badly formatted or a wrong MTZ file.")

  # Work out how many reflection records are contained
  # in this MTZ file (= nref X ncols)
  numDataItems <- irdata[[2]] - 20 - 1

  # Load all reflection records
  if (irdata[[3]] == 1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4)
  if (irdata[[3]] == -1) reflnData <- readBin(f, "numeric", n=numDataItems, size=4,endian="big")

  # Change "NaN" into "NA" which are properly handled in R
  idx <- which(is.nan(reflnData))
  if (length(idx) > 0) reflnData[idx] <- NA

  # Load header information
  hData <- .readH(f,message)
  hF <- hData[[2]]            # hF == 1 means no batch headers

  # Turn reflnData into a dataframe where columns are named
  tmpmatrix <- matrix(data=reflnData,nrow=hData[[1]]$NCOL[2],ncol=hData[[1]]$NCOL[1],
  byrow=TRUE,dimnames=list(NULL,hData[[1]]$COLUMN$labels))
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
  data <- list(reflections=reflnData,header=hData[[1]],
               batch_header=bhData)

  # Create "mtz" class
  #class(data) <- "mtz"

  return(data)
}


#' Write data to an MTZ file
#'
#' Write reflections and experimental information
#' to an MTZ file
#'
#' @param reflections A data frame containing all reflection
#'                    records in columns. This is usually derived
#'                    from modifications of a previously existing
#'                    data frame obtained using
#'                    \code{\link{readMTZ}}.
#' @param header A list whose components are other R objects. This
#'               is normally derived from the reading of another
#'               MTZ file using \code{\link{readMTZ}}. See further
#'               details at \code{\link{readMTZHeader}}.
#' @param batch_header A named list including information at data
#'                     collection time. This component is present
#'                     only for raw (unmerged) intensity data
#'                     produce after the diffraction images
#'                     integration. Merged MTZ reflection files
#'                     have \code{batch_header=NULL}.
#'                     Names and types depend on
#'                     the type of experiment (more information
#'                     on this can be found at
#'                     \href{http://www.ccp4.ac.uk}{CCP4}.)
#' @param filename A character string. The path to a valid mtz
#'                 file. If a file with the same name exists, it
#'                 will be deleted.
#' @param title A character string. The character string
#'              associated with the TITLE keyword in an MTZ file.
#'              This feature makes it easy to quickly identify the
#'              data file in \href{http://www.ccp4.ac.uk}{CCP4}
#'              programs. Default (NULL) is for the output file
#'              to have the same title as the input file.
#' @return A logical variable. If TRUE, the function has written
#'         a successful MTZ file.
#'
#' @examples
#' # Read the 1dei_phases data included in the package
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei_phases.mtz")
#' lMTZ <- readMTZ(filename)
#'
#' # Change dataset name
#' print(lMTZ$header$DATASET)
#' lMTZ$header$DATASET[2,2] <- "New CRY dataset"
#'
#' # Add one HISTORY line (string has to be 80-letters long)
#' addhist <- "From CRY 0.3.0 - run on Apr 2 20:12:00 2021"
#' n <- nchar(addhist)
#' nblanks <- 80-n
#' for (i in 1:nblanks) addhist <- paste0(addhist," ")
#' lMTZ$header$HISTORY <- c(lMTZ$header$HISTORY,addhist)
#'
#' # Write to a new MTZ file
#' wd <- tempdir()
#' fname <- file.path(wd,"new.mtz")
#' writeMTZ(lMTZ$reflections,lMTZ$header,fname)
#'
#' @export
writeMTZ <- function(reflections,header,filename,
                     title=NULL,
                     batch_header=NULL){
  # Create a connection to binary file
  f <- file(filename, open="wb")

  # Hex version of ascii "MTZ "
  rawNULL <- as.raw(0)
  rawBlank <- charToRaw(" ")
  rawMTZ <- c(as.raw(77),as.raw(84),as.raw(90),rawBlank)

  # Write MTZ file signature
  writeBin(rawMTZ,f)

  # Determine data length
  dl <- length(reflections[,1])
  cl <- length(reflections)
  data_length <- dl*cl

  # Location of header
  #headerLoc <- 80+data_length*4+1
  headerLoc <- data_length+20+1
  storage.mode(headerLoc) <- "integer"  # Turn headerLoc into
                                        # integer for correct
                                        # binary storage
  writeBin(headerLoc,f,size=4)

  # Machine stamp (to be improved)
  machineStamp <- "DA"
  writeChar(machineStamp,f,nchars=2)

  # To get to "DA" we have used 8+2 bytes. Reflections will be
  # written starting from byte 81. So we need to add 70 bytes.
  # Fill with NULLs for 69 bytes (2 to complete machine stamp,
  # plus 68 to get to a total of 80. We write one less as this
  # is automatically added when connection is closed)
  for (i in 1:69) writeBin(rawNULL,f)

  # Data (data frame -> matrix -> transpose(matrix) ->
  # vector -> write to binary file)
  #numDataItems <- (headerLoc-80-1)/4
  numDataItems <- data_length
  linedata <- as.vector(t(as.matrix(reflections)))
  writeBin(linedata,f,size=4)
  rm(linedata)

  # Build header 80-characters lines and write to binary file
  fullhdata <- ""
  hdata <- "VERS MTZ:V1.1                                                                   "
  fullhdata <- paste(fullhdata,hdata,sep="") # MTZ stamp
  if (is.null(title)) title <- header$TITLE
  hdata <- paste("TITLE",title)
  nblanks <- 80-nchar(hdata)
  for (i in 1:nblanks) hdata <- paste(hdata," ",sep="")
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # TITLE
  hdata <- sprintf("NCOL%9d%13d%9d                                             ",
                   header$NCOL[1],header$NCOL[2],header$NCOL[3])
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # NCOL
  hdata <- sprintf("CELL %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f               ",
                   header$CELL[1],header$CELL[2],header$CELL[3],header$CELL[4],
                   header$CELL[5],header$CELL[6])
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # CELL
  hdata <- sprintf("SORT %4d%4d%4d%4d%4d                                                       ",
                   header$SORT[1],header$SORT[2],header$SORT[3],
                   header$SORT[4],header$SORT[5])
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SORT
  hdata <- sprintf("SYMINF %3d%3d %1s%6d",
                   header$SYMINF[[1]],header$SYMINF[[2]],header$SYMINF[[3]],
                   header$SYMINF[[4]])
  nblanks <- 23-nchar(header$SYMINF[[5]])
  stmp <- ""
  for (i in 1:nblanks) stmp <- paste(stmp," ",sep="")
  stmp <- paste(stmp,header$SYMINF[[5]],sep="")
  hdata <- paste(hdata,stmp,sprintf("%6s                              ",
                                    header$SYMINF[[6]]),sep="")
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SYMINF
  for (i in 1:length(header$SYMM)) {
    hdata <- header$SYMM[i]
    fullhdata <- paste(fullhdata,hdata,sep="")                                                  # SYMM
  }
  hdata <- sprintf("RESO %8.6f   %8.6f                                                        ",
                   header$RESO[1],header$RESO[2])
  fullhdata <- paste(fullhdata,hdata,sep="")                                                  # RESO
  hdata <- "VALM NAN                                                                        "
  fullhdata <- paste(fullhdata,hdata,sep="")  # VALM

  # Complex part of header COLUMN+COLSRC+COL
  for (i in 1:length(header$COLUMN[,1])) {
    linea1 <- header$COLUMN$labels[i]
    linea2 <- header$COLUMN$types[i]
    linea5 <- as.numeric(header$COLUMN$id[i])
    # String depends on type of data
    if (header$COLUMN$types[i] == "H" |
        header$COLUMN$types[i] == "P" |
        header$COLUMN$types[i] == "B" |
        header$COLUMN$types[i] == "Y" |
        header$COLUMN$types[i] == "I") {
    linea3 <- as.integer(header$COLUMN$min[i])
    linea4 <- as.integer(header$COLUMN$max[i])
    hdata <- sprintf("COLUMN  %-30s%1s  %16d  %16d  %3d",
            linea1,linea2,linea3,linea4,linea5)
    } else {
      linea3 <- as.numeric(header$COLUMN$min[i])
      linea4 <- as.numeric(header$COLUMN$max[i])
      hdata <- sprintf("COLUMN  %-30s%1s  %16.3f  %16.3f  %3d",
            linea1,linea2,linea3,linea4,linea5)
    }
    #print(hdata)
    fullhdata <- paste0(fullhdata,hdata)    # COLUMN
    linea1 <- header$COLSRC$labels[i]
    linea2 <- header$COLSRC$created[i]
    linea3 <- as.numeric(header$COLSRC$id[i])
    hdata <- sprintf("COLSRC %-31s%27s %14d",
            linea1,linea2,linea3)
    #print(hdata)
    fullhdata <- paste0(fullhdata,hdata)    # COLSRC
  }
  hdata <- sprintf("NDIF %8d                                                                   ",
                   header$NDIF)
  fullhdata <- paste(fullhdata,hdata,sep="") # NDIF
  for (i in 1:header$NDIF) {
    hdata <- sprintf("PROJECT%8d %-64s",
     header$PROJECT$id[i],header$PROJECT$pname[i])
    fullhdata <- paste(fullhdata,hdata,sep="")
    hdata <- sprintf("CRYSTAL%8d %-64s",
     header$CRYSTAL$id[i],header$CRYSTAL$cname[i])
    fullhdata <- paste(fullhdata,hdata,sep="")
    hdata <- sprintf("DATASET%8d %-64s",
                     header$DATASET$id[i],header$DATASET$dname[i])
    fullhdata <- paste(fullhdata,hdata,sep="")
    hdata <- sprintf("DCELL%10d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f    ",
                     header$DCELL$id[i],header$DCELL$a[i],
                     header$DCELL$b[i],header$DCELL$c[i],
                     header$DCELL$alpha[i],header$DCELL$beta[i],
                     header$DCELL$gamma[i])
    fullhdata <- paste(fullhdata,hdata,sep="")
    hdata <- sprintf("DWAVEL%9d %10.5f                                                      ",
                     header$DWAVEL$id[i],header$DWAVEL$lambda[i])
    fullhdata <- paste(fullhdata,hdata,sep="")
  }
  plength <- length(header$BATCH)
  dnd <- 0
  while (plength > 0) {
    dnd <- dnd+1
    rem <- plength%%12
    istart <- (dnd-1)*12+1
    iend <- istart+11
    hdata <- sprintf("BATCH ")
    if (plength >= 12) {
      for (i in istart:iend) hdata <- paste(hdata,
                  sprintf("%6d",header$BATCH[i]),sep="")
      hdata <- paste(hdata,sprintf("  "),sep="")
      fullhdata <- paste(fullhdata,hdata,sep="")
    }
    if (plength < 12) {
      for (i in istart:length(header$BATCH))
        hdata <- paste(hdata,sprintf("%6d",header$BATCH[i]),
                       sep="")
      iend <- 12-rem
      for (i in 1:iend) hdata <- paste(hdata,
                                       sprintf("      "),sep="")
      hdata <- paste(hdata,sprintf("  "),sep="")
      fullhdata <- paste(fullhdata,hdata,sep="")
    }
    plength <- plength-12
  }
  hdata <- sprintf("END                                                                             ")
  fullhdata <- paste(fullhdata,hdata,sep="")
  #hdata <- sprintf("MTZHIST%4d                                                                     ",(length(data[[2]]$HISTORY)-1))
  #fullhdata <- paste(fullhdata,hdata,sep="")
  if (!is.null(header$HISTORY)) {
    for (i in 1:length(header$HISTORY)) {
      hdata <- header$HISTORY[i]
      fullhdata <- paste(fullhdata,hdata,sep="")                                                # HISTORY
    }
  }

  # End of MTZ or beginning of batch headers
  if (header$NCOL[3] == 0) {
  hdata <- sprintf("MTZENDOFHEADERS                                                                 ")
    fullhdata <- paste(fullhdata,hdata,sep="")
    writeChar(fullhdata,f,eos=NULL)

    # Close connection
    close(f)
  }

  # Conversion (just copying) to re-use old code
  data <- list(NULL,NULL,batch_header)

  # Carry on filling batch_header (if unmerged file)
  if (header$NCOL[3] != 0) {
    hdata <- sprintf("MTZBATS                                                                         ")
    fullhdata <- paste(fullhdata,hdata,sep="")
    writeChar(fullhdata,f,eos=NULL)

    # All BATCH HEADER stuff is written from here onward
    nbatches <- length(data[[2]]$BATCH)
    for (num in data[[2]]$BATCH[1]:data[[2]]$BATCH[nbatches])
    {
      hdata <- sprintf("BH %8d%8d%8d%8d                                             ",num,data[[3]][[num]]$NWORDS,data[[3]][[num]]$NINTGR,data[[3]][[num]]$NREALS)
      writeChar(hdata,f,eos=NULL)
      hdata <- paste("TITLE",data[[3]][[num]]$TITLE)                                                             # TITLE
      nblanks <- 80-nchar(hdata)
      for (i in 1:nblanks) hdata <- paste(hdata," ",sep="")
      writeChar(hdata,f,eos=NULL)
      ## INTEGERS
      writeBin(data[[3]][[num]]$NWORDS,f,size=4)                                                                 # NWORDS
      writeBin(data[[3]][[num]]$NINTGR,f,size=4)                                                                 # NINTGR
      writeBin(data[[3]][[num]]$NREALS,f,size=4)                                                                 # NREALS
      writeBin(data[[3]][[num]]$IORTYP,f,size=4)                                                                 # IORTYP
      writeBin(data[[3]][[num]]$LBCELL,f,size=4)                                                                 # LBCELL (6)
      writeBin(data[[3]][[num]]$MISFLG,f,size=4)                                                                 # MISFLG
      writeBin(data[[3]][[num]]$JUMPAX,f,size=4)                                                                 # JUMPAX
      writeBin(data[[3]][[num]]$NCRYST,f,size=4)                                                                 # NCRYST
      writeBin(data[[3]][[num]]$LCRFLG,f,size=4)                                                                 # LCRFLG
      writeBin(data[[3]][[num]]$LDTYPE,f,size=4)                                                                 # LDTYPE
      writeBin(data[[3]][[num]]$JSCAX,f,size=4)                                                                  # JSCAX
      writeBin(data[[3]][[num]]$NBSCAL,f,size=4)                                                                 # NBSCAL
      writeBin(data[[3]][[num]]$NGONAX,f,size=4)                                                                 # NGONAX
      writeBin(data[[3]][[num]]$LBMFLG,f,size=4)                                                                 # LBMFLG
      writeBin(data[[3]][[num]]$NDET,f,size=4)                                                                   # NDET
      writeBin(data[[3]][[num]]$LBSETID,f,size=4)                                                                # LBSETID
      writeBin(data[[3]][[num]]$INTPAD,f,size=4)                                                                 # INTPAD (8)
      ## REALS
      writeBin(data[[3]][[num]]$CELL,f,size=4)                                                                   # CELL (6)
      writeBin(as.vector(data[[3]][[num]]$UMAT),f,size=4)                                                        # UMAT (9)
      writeBin(as.vector(data[[3]][[num]]$PHIXYZ),f,size=4)                                                      # PHIXYZ (6)
      writeBin(data[[3]][[num]]$CRYDAT$ETAD,f,size=4)                                                            # CRYDAT$ETAD
      writeBin(data[[3]][[num]]$CRYDAT$ETADH,f,size=4)                                                           # CRYDAT$ETADH
      writeBin(data[[3]][[num]]$CRYDAT$ETADV,f,size=4)                                                           # CRYDAT$ETADH
      writeBin(as.vector(data[[3]][[num]]$CRYDAT$GENERIC),f,size=4)                                              # CRYDAT$GENERIC (9)
      writeBin(data[[3]][[num]]$DATUM,f,size=4)                                                                  # DATUM (3)
      writeBin(data[[3]][[num]]$PHISTT,f,size=4)                                                                 # PHISTT
      writeBin(data[[3]][[num]]$PHIEND,f,size=4)                                                                 # PHIEND
      writeBin(data[[3]][[num]]$SCANAX,f,size=4)                                                                 # SCANAX (3)
      writeBin(data[[3]][[num]]$TIME1,f,size=4)                                                                  # TIME1
      writeBin(data[[3]][[num]]$TIME2,f,size=4)                                                                  # TIME2
      writeBin(data[[3]][[num]]$BSCALE,f,size=4)                                                                 # BSCALE
      writeBin(data[[3]][[num]]$BBFAC,f,size=4)                                                                  # BBFAC
      writeBin(data[[3]][[num]]$SDBSCL,f,size=4)                                                                 # SDBSCL
      writeBin(data[[3]][[num]]$SDBFAC,f,size=4)                                                                 # SDBFAC
      writeBin(data[[3]][[num]]$PHIRANGE,f,size=4)                                                               # PHIRANGE
      writeBin(data[[3]][[num]]$BATPAD,f,size=4)                                                                 # BATPAD (11)
      writeBin(data[[3]][[num]]$E1,f,size=4)                                                                     # E1 (3)
      writeBin(data[[3]][[num]]$E2,f,size=4)                                                                     # E2 (3)
      writeBin(data[[3]][[num]]$E3,f,size=4)                                                                     # E3 (3)
      writeBin(data[[3]][[num]]$GONPAD,f,size=4)                                                                 # GONPAD (12)
      writeBin(data[[3]][[num]]$SOURCE,f,size=4)                                                                 # SOURCE (3)
      writeBin(data[[3]][[num]]$S0,f,size=4)                                                                     # S0 (3)
      writeBin(data[[3]][[num]]$BEMDAT$ALAMBD,f,size=4)                                                          # BEMDAT$ALAMBD
      writeBin(data[[3]][[num]]$BEMDAT$DELAMB,f,size=4)                                                          # BEMDAT$DELAMB
      writeBin(data[[3]][[num]]$BEMDAT$DELCOR,f,size=4)                                                          # BEMDAT$DELCOR
      writeBin(data[[3]][[num]]$BEMDAT$DIVHD,f,size=4)                                                           # BEMDAT$DIVHD
      writeBin(data[[3]][[num]]$BEMDAT$DIVVD,f,size=4)                                                           # BEMDAT$DIVVD
      writeBin(data[[3]][[num]]$BEMDAT$rest,f,size=4)                                                            # rest of BEMDAT (20)
      writeBin(data[[3]][[num]]$DX1,f,size=4)                                                                    # DX1
      writeBin(data[[3]][[num]]$THETA1,f,size=4)                                                                 # THETA1
      writeBin(as.vector(data[[3]][[num]]$DETLM1),f,size=4)                                                      # DETLM1 (4)
      writeBin(data[[3]][[num]]$DX2,f,size=4)                                                                    # DX2
      writeBin(data[[3]][[num]]$THETA2,f,size=4)                                                                 # THETA2
      writeBin(as.vector(data[[3]][[num]]$DETLM2),f,size=4)                                                      # DETLM2 (4)
      writeBin(data[[3]][[num]]$DETPAD,f,size=4)                                                                 # DETPAD (33)

      # Last line "BHCH"
      if (is.null(data[[3]][[num]]$GONLAB)) batchData <- "BHCH                                                                            "
      if (!is.null(data[[3]][[num]]$GONLAB))
      {
        lengthGONLAB <- length(data[[3]][[num]]$GONLAB)
        batchData <- "BHCH "
        for (i in 1:lengthGONLAB) batchData <- paste(batchData,sprintf("%-9s%-9s%-9s",
                                                                       data[[3]][[num]]$GONLAB[1],data[[3]][[num]]$GONLAB[2],data[[3]][[num]]$GONLAB[3]),sep="")
        nblanks <- 80-nchar(hdata)
        for (i in 1:nblanks) batchData <- paste(batchData," ",sep="")
      }
      writeChar(batchData,f,eos=NULL)
    }
    # End of MTZ file
    batchData <- "MTZENDOFHEADERS                                                                 "
    writeChar(batchData,f,eos=NULL)

    # Close connection
    close(f)
  }
}





#' Mean and standard deviation in resolution shells.
#'
#' Calculates averages and standard deviations of the input vector quantity for all
#' reflections, corresponding to shells of resolution.
#'
#' @param nbin A positive integer. The number of resolution shells.
#' @param resos A vector of real quantities. These are the resolutions (in angstroms)
#'  corresponding to the data vector, II. If the data vector is missing, the averages
#'  will be computed just for resos.
#' @param II A vector of real quantities. This is the key quantity whose averages
#'  and standard deviations are calculated. If \code{II} is set to \code{NULL}, resolutions
#'  averages and standard deviations will be the calculated quantities.
#' @param m Minimum (highest) resolution (in angstroms). Data with resolution smaller than m
#'  will be ignored when calculating the averages.
#' @param M Maximum (lowest) resolution (in angstroms). Data with resolution larger than M
#'  will be ignored when calculating the averages.
#' @return A named list of length 4. \code{counts} is a vector of integers, the number of
#'  reflections in each resolution shell. \code{mids} is the representative inverse resolution
#'  for each resolution shell; the value is decided by the function \code{hist}. \code{ave}
#'  is the average value in each resolution shell and \code{sd} is the corresponding standard
#'  deviation.
#' @details Binning is done with inverse resolutions in order to have lower resolutions
#'  correspond to small numbers and high resolutions to large numbers. The output
#'  \code{mids}, \code{ave} and \code{sd} correspond to inverse resolutions.
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei_phases.mtz")
#' lmtz <- readMTZ(filename)
#' hkl <- lmtz$reflections[,1:3]
#' II <- lmtz$reflections[,4]
#' cpars <- lmtz$header$CELL
#' resos <- hkl_to_reso(hkl[,1],hkl[,2],hkl[,3],
#'                      cpars[1],cpars[2],cpars[3],
#'                      cpars[4],cpars[5],cpars[6])
#' ltmp <- avei_vs_res(20,resos,II)
#' plot(ltmp$mids,ltmp$ave,type="b",pch=16)
#' @export
avei_vs_res <- function(nbin,resos,II=NULL,m=max(resos),M=min(resos))
{
  # All calculations done with reciprocal of resos
  ss <- (1/resos)^2
  m <- (1/m)^2
  M <- (1/M)^2

  # Resize vectors if needed
  idx <- which(ss >= m & ss <= M)

  # Return NULL if no elements are included inside min and max range
  if (length(idx) == 0)
  {
    return(list(counts=NULL,mids=NULL,ave=NULL,sd=NULL))
  }

  # If there are elements carry on
  ss <- ss[idx]
  if (length(II) > 0) II <- II[idx]

  # Break points
  breaks <- seq(m,M,length=(nbin+1))

  # Histogram
  hst <- hist(ss,breaks=breaks,plot=FALSE)

  # Partition intensity in resolution bins
  if (length(II) > 0)
  {
    aveI <- c()
    sdI  <- c()
    for (i in 1:nbin)
    {
      idx <- which(ss >= breaks[i] & ss < breaks[i+1])
      if (i == nbin) idx <- which(ss >= breaks[i] & ss <= breaks[i+1])
      if (length(idx) > 0)
      {
        aveI <- c(aveI,mean(II[idx],na.rm=TRUE))
        sdI  <- c( sdI,  sd(II[idx],na.rm=TRUE))
      }
      if (length(idx) == 0)
      {
        aveI <- c(aveI,NA)
        sdI <- c(sdI,NA)
      }
    }
  }
  if (length(II) == 0)
  {
    aveI <- NULL
    sdI <- NULL
  }

  return(list(counts=hst$counts,mids=hst$mids,ave=aveI,sd=sdI))
}



#########
### Auxiliary functions (not exported)
#########


# Cuts string of words into individual words without blanks
.tagRead <- function(stringa)
{
  # Given a string whose words are separated by blanks, returns
  # a list whose first element is all that matter. This is a vector
  # containing the individual words composing the string.
  blanks <- "[[:blank:]]+"
  ltmp <- strsplit(stringa,blanks)
  if (length(ltmp) > 0) {
    return(ltmp[[1]])
  } else {
    return(ltmp)
  }
}


# Reads the identification record of an MTZ binary file 1) "MTZ " 2) Machine
# stamp. It returns a named list of length 3. The first element of the list is called
# errF and can have values 0 and 1: 0 means that the identification record is
# ok, otherwise errF will be 1. The second element is called headerLoc and is
# the position in the mtz binary file of the start of the header section. The
# third field id called edn and is a flag to indicate the machine's endian.
# Normally edn = 1; for, rare, big endian edn = -1. This function also
# positions the connection flow for file reading at the start of the header
# block in the mtz binary file.
.readIR <- function(f)
{
  # Reads identification record of an MTZ binary file.
  # 1) "MTZ "
  # 2) Machine stamp

  # Hex version of ascii "MTZ "
  rawBlank <- charToRaw(" ")
  rawMTZ <- c(as.raw(77),as.raw(84),as.raw(90),rawBlank)

  # errF is zero if identification record is OK; otherwise is 1
  errF <- 0

  # Check initial word is "MTZ "
  mtzStamp <- readBin(f,"raw",n=4)
  for (i in 1:4)
  {
    if (mtzStamp[i] != rawMTZ[i])
    {
      errF <- 1
      return(list(errF,NULL))
    }
  }

  # Location of start for data records. If subsequently-read
  # machineStamp is not "DA", then use big endian
  edn <- 1
  headerLoc <- readBin(f,"integer",n=1,size=4)

  # Machine stamp
  machineStamp <- readChar(f,4) # number formats of the
                                # architecture f was written on

  # If machineStamp != "DA" re-read headerLoc with big endian
  # and change edn into -1
  if (machineStamp != "DA")
  {
    seek(f,4)
    headerLoc <- readBin(f,"integer",n=1,size=4,endian="big")
    edn <- -1
  }

  # Prepare for data records reading: pPosition connection
  # at beginning of reflection data (21st byte=4*20+1)
  seek(f,80)

  return(list(errF=errF,headerLoc=headerLoc,edn=edn))
}


# Given a connection f corresponding to the binary MTZ file immediately after
# end of reflection records, this function keeps reading the file to extract
# information included in the header.
# This function returns a named list of length 2. The first element of the list
# is called l_data and is another list with as many elements as the various fields
# in the mtz header. The second element is called hF and is an integer. 0
# means that the reading has gone ok.
.readH <- function(f,message)
{
  # Given a connection f corresponding to the binary MTZ file immediately
  # after end of reflection records, this function keeps reading the file
  # to extract information included in the header.
  col_symm <- c()
  col_labels <- c()
  col_types <- c()
  col_min <- c()
  col_max <- c()
  col_id <- c()
  colsrc_labels <- c()
  colsrc_created <- c()
  colsrc_id <- c()
  colgrp <- c()
  d_project <- data.frame(id=1,pname=NA)
  d_crystal <- data.frame(id=1,cname=NA)
  d_dataset <- data.frame(id=1,dname=NA)
  d_dcell <- data.frame(id=1,a=NA,b=NA,c=NA,
                        alpha=NA,beta=NA,gamma=NA)
  d_dwavel <- data.frame(id=1,lambda=NA)
  col_batch <- c()
  l_data <- list()
  endTag=""
  while (endTag != "END")
  {
    hdata <- readChar(f,80)
    #print(hdata)
    fields <- .tagRead(hdata)
    endTag <- fields[1]

    # Distribute values according to keyword
    if (endTag == "TITLE") {
      stmp <- ""
      for (i in 2:length(fields)) {
        stmp <- paste(stmp,fields[i])
      }
      stmp <- trimws(stmp,"both")
      l_data$TITLE <- stmp
    }
    if (endTag == "NCOL") l_data$NCOL <- c(as.integer(fields[2]),           # Number of columns
                                           as.integer(fields[3]),           # Number of reflection records
                                           as.integer(fields[4]))           # Number of batches
    if (endTag == "CELL") l_data$CELL <- c(as.numeric(fields[2]),           # Cell parameter a
                                           as.numeric(fields[3]),           # Cell parameter b
                                           as.numeric(fields[4]),           # Cell parameter c
                                           as.numeric(fields[5]),           # Cell parameter alpha
                                           as.numeric(fields[6]),           # Cell parameter beta
                                           as.numeric(fields[7]))           # Cell parameter gamma
    if (endTag == "SORT") l_data$SORT <- c(as.integer(fields[2]),           # Sort order 1
                                           as.integer(fields[3]),           # Sort order 2
                                           as.integer(fields[4]),           # Sort order 3
                                           as.integer(fields[5]),           # Sort order 4
                                           as.integer(fields[6]))           # Sort order 5
    if (endTag == "SYMINF")
    {
      # Need to treat this case differently from the others
      cplace <- c()
      for (i in 1:nchar(hdata))
      {
        if (substr(hdata,i,i) == "'") cplace <- c(cplace,i)
      }
      if (length(cplace) != 0 & length(cplace) != 2 &
          length(cplace) != 4)
           stop("Wrongly formatted SYMINF line in header")
      if (length(cplace) == 0)
      {
        sgname <- fields[length(fields)-1]
        pgname <- fields[length(fields)]
      }
      if (length(cplace) == 2)
      {
        sgname <- substr(hdata,cplace[1],cplace[2])
        pgname <- fields[length(fields)]
      }
      if (length(cplace) == 4)
      {
        sgname <- substr(hdata,cplace[1],cplace[2])
        pgname <- substr(hdata,cplace[3],cplace[4])
      }
      #sgname <- paste("'",strsplit(hdata,"'",fixed=TRUE)[[1]][2],"'",sep="")
      #pgname <- paste("'",strsplit(hdata,"'",fixed=TRUE)[[1]][4],"'",sep="")
      l_data$SYMINF <- list(as.integer(fields[2]),    # Number of symmetry operations
                            as.integer(fields[3]),    # Number of primitive symmetry operations
                            fields[4] ,    # Lattice type
                            as.integer(fields[5]),    # Space group number
                            sgname ,    # Space group name
                            #fields[7])     # Point group name
                            pgname)     # Point group name
    }
    #if (endTag == "SYMINF") l_data$SYMINF <- list(as.integer(fields[2]),    # Number of symmetry operations
    #                                              as.integer(fields[3]),    # Number of primitive symmetry operations
    #                                              fields[4],                # Lattice type
    #                                              as.integer(fields[5]),    # Space group number
    #                                              fields[6],                # Space group name
    #                                              fields[7])                # Point group name
    if (endTag == "SYMM") col_symm <- c(col_symm,hdata)
    if (endTag == "RESO") l_data$RESO <- c(as.numeric(fields[2]),           # Minimum resolution (stored as 1/d-squared)
                                           as.numeric(fields[3]))           # Maximum resolution (stored as 1/d-squared)
    if (endTag == "NDIF") l_data$NDIF <- c(as.numeric(fields[2]))           # Number of datasets included
    if (endTag == "PROJECT")
    {
      id <- as.integer(fields[2])                                             # ID of dataset
      nome <- paste(fields[3:length(fields)],collapse=" ")                    # Project name
      d_project <- rbind(d_project,data.frame(id=id,pname=nome))
    }
    if (endTag == "CRYSTAL")
    {
      id <- as.integer(fields[2])                                             # ID of dataset
      nome <- paste(fields[3:length(fields)],collapse=" ")                    # Crystal name
      d_crystal <- rbind(d_crystal,data.frame(id=id,cname=nome))
    }
    if (endTag == "DATASET")
    {
      id <- as.integer(fields[2])                                             # ID of dataset
      nome <- paste(fields[3:length(fields)],collapse=" ")                    # Dataset name
      d_dataset <- rbind(d_dataset,data.frame(id=id,dname=nome))
    }
    if (endTag == "DCELL")
    {
      id <- as.integer(fields[2])                                             # ID of dataset
      d_dcell <- rbind(d_dcell,data.frame(id=id,a=as.numeric(fields[3]),
                                          b=as.numeric(fields[4]),
                                          c=as.numeric(fields[5]),
                                          alpha=as.numeric(fields[6]),
                                          beta=as.numeric(fields[7]),
                                          gamma=as.numeric(fields[8]))) # Cell parameters (one set for each dataset)
    }
    if (endTag == "DWAVEL")
    {
      id <- as.integer(fields[2])                                             # ID of dataset
      lambda <- as.numeric(fields[3])
      if (is.na(lambda)) lambda <- -1
      d_dwavel <- rbind(d_dwavel,
                        #                  data.frame(id=id,lambda=as.numeric(fields[3])))       # Wavelength (one for each dataset)
                        data.frame(id=id,lambda=lambda))       # Wavelength (one for each dataset)
    }
    if (endTag == "COLUMN")
    {
      col_labels <- c(col_labels,fields[2])
      col_types <- c(col_types,fields[3])
      col_min <- c(col_min,fields[4])
      col_max <- c(col_max,fields[5])
      col_id <- c(col_id,fields[6])
      if (message) print(hdata)
    }
    if (endTag == "COLSRC")
    {
      colsrc_labels <- c(colsrc_labels,fields[2])
      colsrc_created <- c(colsrc_created,fields[3])
      colsrc_id <- c(colsrc_id,fields[4])
    }
    if (endTag == "COLGRP") colgrp <- c(colgrp,hdata)
    if (endTag == "BATCH")
    {
      col_batch <- c(col_batch,as.integer(fields[2:length(fields)]))
    }
  }

  l_data$SYMM <- col_symm
  if (length(d_project[,1]) != 1) d_project <- na.omit(d_project)
  rownames(d_project) <- 1:length(d_project[,1])
  l_data$PROJECT <- d_project
  if (length(d_crystal[,1]) != 1) d_crystal <- na.omit(d_crystal)
  rownames(d_crystal) <- 1:length(d_crystal[,1])
  l_data$CRYSTAL <- d_crystal
  if (length(d_dataset[,1]) != 1) d_dataset <- na.omit(d_dataset)
  rownames(d_dataset) <- 1:length(d_dataset[,1])
  l_data$DATASET <- d_dataset
  if (length(d_dcell[,1]) != 1) d_dcell <- na.omit(d_dcell)
  rownames(d_dcell) <- 1:length(d_dcell[,1])
  l_data$DCELL <- d_dcell
  if (length(d_dwavel[,1]) != 1) d_dwavel <- na.omit(d_dwavel)
  for (i in 1:length(d_dwavel[,1]))
  {
    if (!is.na(d_dwavel[i,2])) if (d_dwavel[i,2] == -1) d_dwavel[i,2] <- NA
  }
  rownames(d_dwavel) <- 1:length(d_dwavel[,1])
  l_data$DWAVEL <- d_dwavel

  # The .I_null_check() avoids values being turned into factor
  # levels
  l_data$COLUMN <- data.frame(labels=.I_null_check(col_labels),
                              types=.I_null_check(col_types),
                              min=.I_null_check(col_min),
                              max=.I_null_check(col_max),
                              id=.I_null_check(col_id))
  l_data$COLSRC <- data.frame(labels=.I_null_check(colsrc_labels),
                    created=.I_null_check(colsrc_created),
                              id=.I_null_check(colsrc_id))
  l_data$COLGRP <- colgrp

  l_data$BATCH <- col_batch

  # History
  col_history <- c()
  bheaderTag=" "
  while (bheaderTag != "MTZBATS" &
         bheaderTag != "MTZENDOFHEADERS")
  {
    hdata <- readChar(f,80)
    col_history <- c(col_history,hdata)
    bheaderTag <- .tagRead(hdata)[1]
  }
  if (substr(col_history[1],1,6) != "MTZEND") {
   l_data$HISTORY <- col_history[1:(length(col_history)-1)]
  } else {
    l_data$HISTORY <- NULL
  }

  if (bheaderTag == "MTZENDOFHEADERS")
  {
    hF <- 1
  }
  else
  {
    hF <- 0
  }

  data <- list(l_data=l_data,hF=hF)
  return(data)
}


# Given a connection f corresponding to the binary MTZ file immediately
# after "MTZBATS", this function keeps reading the file to extract the
# batch header content.
# It returns a list whose elements include information on each batch header.
# This is the same as the one documented in the MTZLIB library of CCP4
# There are still 19 integers which I can't recognize and, therefore, do not
# output. From the documentation there seem to be included also a character
# X 3, which I can't find in the binary. It returns A list with as many elements
# as the number of batches (images) included in the mtz file. Each list element is,
# itself a named list with all the useful variables stored in batch headers.
.readBH <- function(f)
{
  # Given a connection f corresponding to the binary MTZ file immediately
  # after "MTZBATS", this function keeps reading the file to extract the
  # batch header content.
  # It returns a list whose elements include information on each batch header.
  # This is the same as the one documented in the MTZLIB library of CCP4
  # There are still 19 integers which I can't recognize and, therefore, do not output.
  # From the documentation there seem to be included also a character X 3, which I can't
  # find in the binary.

  # Hex version of ascii NULL and "B"
  rawNULL <- as.raw(0)
  rawB <- as.raw(66)

  # Flag to stop reading batch headers
  bflag <- 0

  # List containing all batch headers
  l_data <- list()

  # Main cycle through all batch headers
  while (bflag != 1)
  {
    # Each batch header is temporarily stored in a list on its own.
    # This is initialized after each cycle
    single_batch <- list()

    batchData <- readChar(f,80)
    batchData <- .tagRead(batchData)              # Initial string starting by "BH" and containing batch number, number of records, number of integers and number of reals
    idx <- as.integer(batchData[2])
    if (batchData[1] == "BH")
    {
      batchData <- readChar(f,80)
      batchData <- .tagRead(batchData)              # Title of batch header
      if (length(batchData) == 1)
      {
        single_batch$TITLE <- ""
      }
      if (length(batchData) > 1)
      {
        parola <- ""
        for (i in 2:length(batchData)) parola <- paste(parola,batchData[i])
        single_batch$TITLE <- parola
      }

      # Numbers for batch header, 29 integers and 156 reals
      batchData <- readBin(f,"integer",n=29,size=4)
      single_batch$NWORDS <- batchData[1]
      single_batch$NINTGR <- batchData[2]
      single_batch$NREALS <- batchData[3]
      single_batch$IORTYP <- batchData[4]
      single_batch$LBCELL <- batchData[5:10]
      single_batch$MISFLG <- batchData[11]
      single_batch$JUMPAX <- batchData[12]
      single_batch$NCRYST <- batchData[13]
      single_batch$LCRFLG <- batchData[14]
      single_batch$LDTYPE <- batchData[15]
      single_batch$JSCAX <- batchData[16]
      single_batch$NBSCAL <- batchData[17]
      single_batch$NGONAX <- batchData[18]
      single_batch$LBMFLG <- batchData[19]
      single_batch$NDET <- batchData[20]
      single_batch$LBSETID <- batchData[21]
      single_batch$INTPAD <- batchData[22:29]
      batchData <- readBin(f,"numeric",n=156,size=4)
      single_batch$CELL <- batchData[1:6]
      single_batch$UMAT <- matrix(data=batchData[7:15],nrow=3,ncol=3)
      single_batch$PHIXYZ <- matrix(data=batchData[16:21],nrow=3,ncol=2)
      single_batch$CRYDAT$ETAD <- batchData[22]
      single_batch$CRYDAT$ETADH <- batchData[23]
      single_batch$CRYDAT$ETADV <- batchData[24]
      single_batch$CRYDAT$GENERIC <- matrix(data=batchData[25:33],nrow=3,ncol=3)
      single_batch$DATUM <- batchData[34:36]
      single_batch$PHISTT <- batchData[37]
      single_batch$PHIEND <- batchData[38]
      single_batch$SCANAX <- batchData[39:41]
      single_batch$TIME1 <- batchData[42]
      single_batch$TIME2 <- batchData[43]
      single_batch$BSCALE <- batchData[44]
      single_batch$BBFAC <-  batchData[45]
      single_batch$SDBSCL <- batchData[46]
      single_batch$SDBFAC <- batchData[47]
      single_batch$PHIRANGE <- batchData[48]
      single_batch$BATPAD <- batchData[49:59]
      single_batch$E1 <- batchData[60:62]
      single_batch$E2 <- batchData[63:65]
      single_batch$E3 <- batchData[66:68]
      single_batch$GONPAD <- batchData[69:80]
      single_batch$SOURCE <- batchData[81:83]
      single_batch$S0 <- batchData[84:86]
      single_batch$BEMDAT$ALAMBD <- batchData[87]
      single_batch$BEMDAT$DELAMB <- batchData[88]
      single_batch$BEMDAT$DELCOR <- batchData[89]
      single_batch$BEMDAT$DIVHD  <- batchData[90]
      single_batch$BEMDAT$DIVVD  <- batchData[91]
      single_batch$BEMDAT$rest   <- batchData[92:111]
      # First detector
      single_batch$DX1 <- batchData[112]
      single_batch$THETA1 <- batchData[113]
      single_batch$DETLM1 <- matrix(data=batchData[114:117],nrow=2,ncol=2)
      # Second detector
      single_batch$DX2 <- batchData[118]
      single_batch$THETA2 <- batchData[119]
      single_batch$DETLM2 <- matrix(data=batchData[120:123],nrow=2,ncol=2)

      single_batch$DETPAD <- batchData[124:156]

      # Last line
      batchData <- readChar(f,80)
      fields <- .tagRead(batchData)
      if (length(fields) > 1) single_batch$GONLAB <- fields[2:length(fields)]
      if (length(fields) == 1) single_batch$GONLAB <- NULL

      # Store this batch header
      l_data[[idx]] <- single_batch
    }
    else
    {
      bflag <- 1
    }
  }

  return(l_data)
}


# Patch function taking care of changes with which R makes use of the I() variant
# to read characters as factors.
.I_null_check <- function(x) {
  if (is.null(x)) {
    return(list())
  } else {
    return(I(x))
  }
}


#' Generate Miller indices
#'
#' Function to create a data frame with complete set of Miller
#' indices, up to a given resolution (in angstroms).
#'
#' Miller indices are named H, K, L in the data frame. Only
#' values of (H,K,L) corresponding to a resolution d(h,k,l) >=
#' reso (in angstroms), are included. The full list does not
#' include systematic absences corresponding to the specific
#' symmetry of the crystal.
#'
#' @param uc An object of class "unit_cell".
#' @param SG A character string or a number indicating the
#'           extended Hermann-Mauguin symbol for the space group.
#' @param reso A real number. The highest data resolution, in
#'             angstroms.
#' @return hkl A data frame with columns H, K, L corresponding
#'             to the three Miller indices, and a columns S
#'             corresponding to their inverse resolutions (in
#'             angstroms).
#' @examples
#' # C 2 monoclinic space group
#' SG <- "C 1 2 1"
#'
#' # Create an arbitrary cell compatible with C 2
#' uc <- unit_cell(10,15,10,90,110,90)
#'
#' # Generate Miller indices to 5 angstroms resolution
#' reso <- 5
#' hkl <- generate_miller(uc,SG,reso)
#'
#' # Display first 10 indices
#' hkl[1:10,]
#'
#' @export
generate_miller <- function(uc,SG,reso) {
  # check input
  ans <- check_unit_cell_validity(uc)
  if (!ans) {
    msg <- "'uc' is not a valid 'unit_cell' object.\n"
    cat(msg)

    return(NULL)
  }
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

  # Cell parameters
  a <- uc$a
  b <- uc$b
  c <- uc$c
  aa <- uc$alpha
  bb <- uc$beta
  gg <- uc$gamma

  # Inverse, squared resolution
  ss <- (1/reso)^2

  # Coefficients of squared resolution function
  ctmp <- squared_resolution_coeffs(a,b,c,aa,bb,gg)
  ca <- ctmp[1]
  cb <- ctmp[2]
  cc <- ctmp[3]
  cd <- ctmp[4]
  ce <- ctmp[5]
  cf <- ctmp[6]

  # Find max h, k, l to generate all (h,k,l) within minimal box

  # Coefficient determinant (common to both h, and k and l)
  D <- ca*(cb*cc-cf^2)-cd*(cc*cd-ce*cf)+ce*(cd*cf-cb*ce)

  # Max h
  A <- (cb*cc-cf^2)/D
  B <- (ce*cf-cc*cd)/D
  C <- (cd*cf-cb*ce)/D
  t <- sqrt(ss/(ca*A^2+cb*B^2+cc*C^2+2*cd*A*B+2*ce*A*C+2*cf*B*C))
  max_h <- ceiling(A*t)

  # Max k
  A <- (ce*cf-cc*cd)/D
  B <- (ca*cc-ce^2)/D
  C <- (cd*ce-ca*cf)/D
  t <- sqrt(ss/(ca*A^2+cb*B^2+cc*C^2+2*cd*A*B+2*ce*A*C+2*cf*B*C))
  max_k <- ceiling(B*t)

  # Max l
  A <- (cd*cf-cb*ce)/D
  B <- (cd*ce-ca*cf)/D
  C <- (ca*cb-cd^2)/D
  t <- sqrt(ss/(ca*A^2+cb*B^2+cc*C^2+2*cd*A*B+2*ce*A*C+2*cf*B*C))
  max_l <- ceiling(C*t)

  # Dataframe with all h, k, l
  hkl <- expand.grid(H=-max_h:max_h,K=-max_k:max_k,
                     L=-max_l:max_l)

  # Add inverse of resolution (s) as last column of data frame
  s <- 1/hkl_to_reso(hkl$H,hkl$K,hkl$L,a,b,c,aa,bb,gg)
  hkl <- cbind(hkl,data.frame(S=s))

  # Strictly reflections with resolution lower than reso
  hkl <- hkl[hkl$S <= 1/reso,]

  # Eliminate systematic absences from data frame
  hkl <- deplete_systematic_absences(hkl,SG)

  # Re-number rows
  row.names(hkl) <- 1:length(hkl[,1])

  return(hkl)
}

# Given cell parameters return coefficients for the expression
# of squared resolution:
#
#  s^2 = a*h^2+b*k^2+c*l^2+2*d*h*k+2*e*h*l+2*f*k*l
#
squared_resolution_coeffs <- function(a,b,c,aa,bb,cc) {
  # Copy original names into new ones (to adapt old code)
  cell_a <- a
  cell_b <- b
  cell_c <- c
  cell_alpha <- aa
  cell_beta <- bb
  cell_gamma <- cc

  ctmp <- unname(lattice_stuff(cell_a,cell_b,cell_c,
                               cell_alpha,cell_beta,cell_gamma))
  den <- 1-ctmp[4]^2-ctmp[5]^2-
         ctmp[6]^2+2*ctmp[4]*ctmp[5]*ctmp[6]
  a <- ctmp[1]^2/(cell_a^2*den)
  b <- ctmp[2]^2/(cell_b^2*den)
  c <- ctmp[3]^2/(cell_c^2*den)
  d <- (ctmp[4]*ctmp[5]-ctmp[6])/(cell_a*cell_b*den)
  e <- (ctmp[4]*ctmp[6]-ctmp[5])/(cell_a*cell_c*den)
  f <- (ctmp[5]*ctmp[6]-ctmp[4])/(cell_b*cell_c*den)

  return(c(a,b,c,d,e,f))
}


#' Deplete systematic absences
#'
#' Remove systematically-absent reflections from a data frame
#' in which Miller indices are in the first three columns.
#' The systematically-absent reflections are determined by the
#' specific space group.
#'
#' Crystallography symmetry forces constraints on data from
#' x-ray diffraction. One of these constraints consists in the
#' full cancellation of reflections with certain Miller indices.
#' It is said that the reflection with that specific Miller index
#' is systematically absent. For example, in data corresponding
#' to a crystal with space group C 2, general reflections like
#' (h,k,l) must obey h+k=2n (even number). Thus, the Miller
#' indices (2,3,1) are a systematic absence because 2+3=5 (odd).
#'
#' @param hkl A data frame with first three columns H, K, L
#'        corresponding to the three Miller indices. This is
#'        normally the 'record' data frame in an object of
#'        class "merged_reflections".
#' @param SG A character. The extended Hermann-Mauguin symbol
#'           of the crystallographic space group.
#' @return hkl The same data frame acquired from input, depleted
#'             of all systematic absences.
#'
#' @examples
#' # C 2 monoclinic space group
#' SG <-"C 1 2 1"
#'
#' # Create an arbitrary cell compatible with C 2
#' uc <- unit_cell(10,15,10,90,110,90)
#'
#' # Crete the related reciprocal cell
#' ruc <- create_rec_unit_cell(uc)
#'
#' # Create a full data frame of Miller indices
#' hkl <- expand.grid(H=-4:4,K=-4:4,L=-4:4)
#'
#' # Get rid of systematic absences
#' new_hkl <- deplete_systematic_absences(hkl,SG)
#'
#' # Compare first 10 items of original and depleted arrays
#' hkl[1:10,]
#' new_hkl[1:10,]
#'
#' @export
deplete_systematic_absences <- function(hkl,SG) {
  # Verify inputs
  ans <- dim(hkl)
  if(is.null(ans)) {
    msg <- paste("Input 'hkl' needs to be a nXm array,",
                 "with m >= 3.\n")
    cat(msg)

    return(NULL)
  }
  ans <- ans[2] >= 3
  if (!ans) {
    msg <- paste("Input 'hkl' needs to be a nXm array,",
                 "with m >= 3.\n")
    cat(msg)

    return(NULL)
  }
  ans <- is.character(SG)
  if (!ans) {
    msg <- paste("Input needs to be an extended",
                 "Hermann-Mauguin symbol.\n")
    cat(msg)

    return(NULL)
  }
  ans <- findHM(SG)
  if (!is.null(ans)) SG <- ans

  # Make a m X 3 matrix of original Miller indices
  hkl2 <- as.matrix(hkl[,1:3])

  # Delete incorrect Miller indices if they are
  # systematically-absent
  idx <- sysabs(hkl2,SG)
  hkl <- hkl[idx,]

  return(hkl)
}


#' Locate systematic absences
#'
#' Given an mX3 matrix of Miller indices, this function returns
#' those indices corresponding to valid reflections, i.e. to
#' reflections which are not systematic absences.
#'
#' Crystallography symmetry forces constraints on data from
#' x-ray diffraction. One of these constraints consists in the
#' full cancellation of reflections with certain Miller indices.
#' It is said that the reflection with that specific Miller index
#' is systematically absent. For example, in data corresponding
#' to a crystal with space group C 2, general reflections like
#' (h,k,l) must obey h+k=2n (even number). Thus, the Miller
#' indices (2,3,1) are a systematic absence because 2+3=5 (odd).
#'
#' @param hkl An mX3 matrix or a data frame whose rows are the
#'            three integers corresponding to the Miller
#'            indices.
#' @param SG A character. The extended Hermann-Mauguin symbol
#'           of the crystallographic space group.
#' @return idx A vector of integers corresponding to the
#'             position, in the array \code{mhkl}, in which the
#'             Miller indices ARE NOT systematically absent.
#'             The position of systematically-absent reflections
#'             can be found using !idx.
#'
#' @examples
#' # C 2 monoclinic space group (special setting)
#' csym <- cryst_symm(15,set=5)
#' print(csym$SG)
#'
#' # Create a full data frame of Miller indices
#' hkl <- expand.grid(H=-4:4,K=-4:4,L=-4:4)
#'
#' # Index corresponding to valid reflections
#' # (not systematic absences)
#' idx <- sysabs(hkl,csym$SG)
#'
#' # Indices for all reflections
#' fulldx <- 1:length(hkl[,1])
#'
#' # Index corresponding to systematic absences
#' jdx <- fulldx[-idx]
#'
#' # A couple of systematic absences
#' hkl[jdx[1:2],]
#'
#' @export
sysabs  <- function(hkl,SG) {
  # Verify inputs
  ans <- dim(hkl)
  if(is.null(ans)) {
    msg <- "Input 'hkl' needs to be a nX3 array.\n"
    cat(msg)

    return(NULL)
  }
  ans <- ans[2] >= 3
  if (!ans) {
    msg <- paste("Input 'hkl' needs to be a nXm array",
                 " with m >= 3.\n")
    cat(msg)

    return(NULL)
  }
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

  # First add a fourth column with 0 to matrix
  nrefs <- nrow(hkl)
  hkl <- cbind(hkl,matrix(rep(0,times=nrefs),nrow=nrefs))
  colnames(hkl)[4] <- "FLAG"

  # Long list for all space groups and settings

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
  # Nothing to do for sym_number 3
  if (sym_number == 4)
  {
    if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 2) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 5)
  {
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 10) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 11) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 6
  if (sym_number == 7)
  {
    if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 2,3] <- NA
    if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2,3] <- NA
  }
  if (sym_number == 8)
  {
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 10) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 9)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 7)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 8)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 9)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 10)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 11)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 12)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 13)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 14)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 15)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 16)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 17)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 18)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  # Nothing to do for sym_number 10, 11
  if (sym_number == 12)
  {
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 5) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 6) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 7) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 8) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 9) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 13)
  {
    if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 2) hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 4) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 5) hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 7) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 8) hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 9) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 14)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 7)
    {
      hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 8)
    {
      hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 9)
    {
      hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 15)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 7)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 8)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 9)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 10)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 11)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 12)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 13)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 14)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 15)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 16)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 17)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 18)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  #########################################################################################################################
  #########
  #########            ORTHOROMBIC: 16 to 74
  #########
  #########################################################################################################################
  # Nothing to do for sym_number 16
  if (sym_number == 17)
  {
    if (setting == 1) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 2) hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 18
  if (sym_number == 19)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 20)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 21)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- 3
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- 3
    if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- 3
  }
  if (sym_number == 22)
  {
    hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 23) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 24)
  {
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 25
  if (sym_number == 26)
  {
    if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 2) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 27)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 28)
  {
    if (setting == 1) hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 2 | setting == 3) hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 4) hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 5) hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 6) hkl[hkl[,3] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 29)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,2] == 0 & hkl[,1] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,3] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 30)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 31)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 32)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 33)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 34)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 35)
  {
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 36)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 37)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 38)
  {
    if (setting == 1) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 4) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 5) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 6) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 39)
  {
    if (setting == 1)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 40)
  {
    if (setting == 1)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 41)
  {
    if (setting == 1)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 42)
  {
    hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 43)
  {
    if (setting == 1)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%4 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%4 != 0,4] <- NA
    }
  }
  if (sym_number == 44) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 45)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 46)
  {
    if (setting == 1)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  # Nothing to do for sym_number 47
  if (sym_number == 48)
  {
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 49)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 50)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 51)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 52)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 53)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 54)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 55)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 56)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 57)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 58)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 59)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 60)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 61)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 62)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 63)
  {
    if (setting == 1) hkl[hkl[,2] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 2) hkl[hkl[,1] == 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 3) hkl[hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 4) hkl[hkl[,2] == 0 & hkl[,1]%%2 != 0,4] <- NA
    if (setting == 5) hkl[hkl[,1] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 6) hkl[hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    if (setting == 1 | setting == 2) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 3 | setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 5 | setting == 6) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 64)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    }
    if (setting == 1 | setting == 2) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 3 | setting == 4) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting == 5 | setting == 6) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 65)
  {
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,1])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,1])%%2 != 0,4] <- NA
  }
  if (sym_number == 66)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 1) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting == 2) hkl[(hkl[,2]+hkl[,1])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,1])%%2 != 0,4] <- NA
  }
  if (sym_number == 67)
  {
    if (setting == 1 | setting == 2 | setting == 6)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 4 | setting == 5)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
  }
  if (sym_number == 68)
  {
    if (setting == 1 | setting == 2 | setting == 3 | setting == 4)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 5 | setting == 6 | setting == 7)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 8)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 9 | setting == 10 | setting == 11)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 12)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting >= 1 & setting <= 4) hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    if (setting >= 5 & setting <= 8) hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    if (setting >= 9 & setting <= 12) hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 69)
  {
    hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 70)
  {
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (3*(hkl[,1]+hkl[,2]))%%4 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (3*(hkl[,2]+hkl[,3]))%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (3*(hkl[,1]+hkl[,3]))%%4 != 0,4] <- NA
    hkl[(hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 71)
  {
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 72)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 73)
  {
    if (setting == 1)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 74)
  {
    if (setting == 1 | setting == 2)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 3)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    if (setting == 4)
    {
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
      hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    }
    if (setting == 5)
    {
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    }
    if (setting == 6)
    {
      hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
      hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    }
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  #########################################################################################################################
  #########
  #########            TETRAGONAL: 75 to 142
  #########
  #########################################################################################################################
  # Nothing to do for sym_number 75
  if (sym_number == 76) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  if (sym_number == 77) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (sym_number == 78) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  if (sym_number == 79) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 80)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 81
  if (sym_number == 82) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  # Nothing to do for sym_number 83
  if (sym_number == 84) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (sym_number == 85) hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  if (sym_number == 86)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 87) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 88)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 89
  if (sym_number == 90)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 91) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (sym_number == 92)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 93) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (sym_number == 94)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 95) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
  if (sym_number == 96)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 97) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 98)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 99
  if (sym_number == 100)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 101)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 102)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 103)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 104)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 105)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 106)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 107) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 108)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 109)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 110)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 111
  if (sym_number == 112) hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  if (sym_number == 113)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  if (sym_number == 114)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 115
  if (sym_number == 116)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 117)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 118)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 119) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 120)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 121) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 122)
  {
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 123
  if (sym_number == 124)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] ==  hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 125)
  {
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & +hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & +hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 126)
  {
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 127)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
  }
  if (sym_number == 128)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 129)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
  }
  if (sym_number == 130)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 131)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 132)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 133)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 134)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 135)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 136)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & (hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 137)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 138)
  {
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,2]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & (hkl[,1]+hkl[,2])%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 139) hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  if (sym_number == 140)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 141)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & (hkl[,1]+hkl[,3])%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  if (sym_number == 142)
  {
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%4 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] == 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] != 0 & hkl[,3] == 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,1]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & (2*hkl[,1]+hkl[,3])%%4 != 0,4] <- NA
    hkl[(hkl[,1]+hkl[,2]+hkl[,3])%%2 != 0,4] <- NA
  }
  #########################################################################################################################
  #########
  #########            TRIGONAL: 143 to 167
  #########
  #########################################################################################################################
  # Nothing to do for sym_number 143
  if (sym_number == 144) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 145) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 146)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    # Nothing to do for setting 2
  }
  # Nothing to do for sym_number 147
  if (sym_number == 148)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    # Nothing to do for setting 2
  }
  # Nothing to do for sym_number 149
  # Nothing to do for sym_number 150
  if (sym_number == 151) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 152) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 153) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 154) hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%3 != 0,4] <- NA
  if (sym_number == 155)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    # Nothing to do for setting 2
  }
  # Nothing to do for sym_number 156
  # Nothing to do for sym_number 157
  if (sym_number == 158)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 159)
  {
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 160)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    # Nothing to do for setting 2
  }
  if (sym_number == 161)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    # Nothing else for setting 2
  }
  # Nothing to do for sym_number 162
  if (sym_number == 163)
  {
    hkl[hkl[,1] == hkl[,2] & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  # Nothing to do for sym_number 164
  if (sym_number == 165)
  {
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  if (sym_number == 166)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
  }
  if (sym_number == 167)
  {
    if (setting == 1)
    {
      hkl[(-hkl[,1]+hkl[,2]+hkl[,3])%%3 != 0,4] <- NA
    }
    hkl[hkl[,1] == 0 & hkl[,2] != 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] != 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
    hkl[hkl[,1] == 0 & hkl[,2] == 0 & hkl[,3] != 0 & hkl[,3]%%2 != 0,4] <- NA
  }
  #########################################################################################################################
  #########
  #########            HEXAGONAL: 168 to 194
  #########
  #########################################################################################################################
  #########################################################################################################################
  #########
  #########            CUBIC: 195 to 230
  #########
  #########################################################################################################################

  # Delete systematic absences
  idx <- which(complete.cases(hkl))

  return(idx)
}

#' Change COLSRC date and time stamp
#'
#' Function to update the \code{created} column of the data
#' frame \code{COLSRC} with current date and time.
#'
#' The COLSRC data frame of an MTZ header has a column called
#' \code{created} which displays the date and time at which the
#' MTZ file data columns were created. When writing out a
#' modified list obtained from reading an MTZ file, one might
#' want to change the \code{created} column with the current
#' date and time. Other specific types of change can be operated
#' by handling the \code{COLSRC} data frame in an *ad hoc* manner.
#'
#' @param hdr A data frame. The \code{COLSRC} data frame included
#'            in the \code{header} component of the named list
#'            obtained with \code{\link{readMTZ}} or
#'            \code{\link{readMTZHeader}}.
#' @return The \code{hdr} input data frame with the \code{created}
#'         column of the \code{COLSRC} data frame changed to
#'         display the current date and time.
#'
#' @examples
#' # Read a sample MTZ file
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei_phases.mtz")
#' lMTZ <- readMTZ(filename)
#'
#' # Original COLSRC
#' print(lMTZ$header$COLSRC)
#'
#' # Update date and time stamp
#' lMTZ$header <- change_COLSRC(lMTZ$header)
#'
#' # New COLSRC
#' print(lMTZ$header$COLSRC)
#'
#' @export
change_COLSRC <- function(hdr) {
  # Get date and time from system
  tt <- strsplit(as.character(Sys.time())," ")[[1]]
  ss <- strsplit(tt[1],"-")[[1]]
  g <- paste0(ss[3],"/",ss[2],"/",ss[1])
  stmp <- paste0("CREATED_",g,"_",tt[2])

  # Change 'created' field of COLSRC
  for (i in 1:length(hdr$COLSRC[,1])) {
    hdr$COLSRC[i,2] <- stmp
  }

  return(hdr)
}
