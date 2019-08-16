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
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"insulin_merged.mtz")
#' ltmp <- readMTZHeader(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
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
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"insulin_merged.mtz")
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
  data <- list(reflections=reflnData,header=hData[[1]],batch_header=bhData)

  # Create "mtz" class
  #class(data) <- "mtz"

  return(data)
}


#' Mean and standard deviation in resolution shells.
#'
#' Calculates averages and standard deviations of the input vector quantity for all
#' reflections, corresponding to shells of resolution.
#'
#' @param nbin A positive integer. The number of resolution shells.
#' @param II A vector of real quantities. This is the key quantity whose averages
#'  and standard deviations are calculated. If \code{II} is set to \code{NULL}, resolutions
#'  averages and standard deviations will be the calculated quantities.
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
#' filename <- file.path(datadir,"insulin_merged.mtz")
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

  return(ltmp[[1]])
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

  # Location of start for data records. If subsequently-read machineStamp
  # is not "DA", then use big endian
  edn <- 1
  headerLoc <- readBin(f,"integer",n=1,size=4)

  # Machine stamp
  machineStamp <- readChar(f,4) # number formats of the architecture f was
  # written on

  # If machineStamp != "DA" re-read headerLoc with big endian and change edn into -1
  if (machineStamp != "DA")
  {
    seek(f,4)
    headerLoc <- readBin(f,"integer",n=1,size=4,endian="big")
    edn <- -1
  }

  # Prepare for data records reading: pPosition connection at beginning
  # of reflection data (21st byte=4*20+1)
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
.readH <- function(f,messages)
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
  d_dcell <- data.frame(id=1,a=NA,b=NA,c=NA,alpha=NA,beta=NA,gamma=NA)
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
      if (length(cplace) != 0 & length(cplace) != 2 & length(cplace) != 4) stop("Wrongly formatted SYMINF line in header")
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
      if (messages) print(hdata)
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

  # The .I_null_check() avoids values being turned into factor levels
  l_data$COLUMN <- data.frame(labels=.I_null_check(col_labels),types=.I_null_check(col_types),
                              min=.I_null_check(col_min),max=.I_null_check(col_max),
                              id=.I_null_check(col_id))
  l_data$COLSRC <- data.frame(labels=.I_null_check(colsrc_labels),
                              created=.I_null_check(colsrc_created),
                              id=.I_null_check(colsrc_id))
  l_data$COLGRP <- colgrp

  l_data$BATCH <- col_batch

  # History
  col_history <- c()
  bheaderTag=" "
  while (bheaderTag != "MTZBATS" & bheaderTag != "MTZENDOFHEADERS")
  {
    hdata <- readChar(f,80)
    col_history <- c(col_history,hdata)
    bheaderTag <- .tagRead(hdata)[1]
  }
  l_data$HISTORY <- col_history[1:(length(col_history)-1)]

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
