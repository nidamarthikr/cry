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

  # Change "NaN" into "NA" which are properly handled in R
  idx <- which(is.nan(reflnData))
  if (length(idx) > 0) reflnData[idx] <- NA

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
#' @param csym An object of class "cryst_symm".
#' @param reso A real number. The highest data resolution, in
#'             angstroms.
#' @param set An integer number corresponding to the specific
#'            setting for the given space group. Default is 1.
#' @return hkl A data frame with columns H, K, L corresponding
#'             to the three Miller indices.
#' @examples
#' Create a C 2 (monoclinic) space group
#' csym <- cryst_symm("C 2")
#'
#' # Create an arbitrary cell compatible with C 2
#' uc <- unit_cell(10,15,10,90,110,90)
#'
#' # Generate Miller indices to 5 angstroms resolution
#' hkl <- generate_miller(uc,csym,reso)
#'
#' # Display first 10 indices
#' hkl[1:10,]
#'
#' @export
generate_miller <- function(uc,csym,reso,set=1) {
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

  # Strictly reflections with resolution lower than reso
  hkl <- hkl[1/hkl_to_reso(hkl$H,hkl$K,hkl$L,
                               a,b,c,aa,bb,gg) <= 1/reso,]

  # Add inverse of resolution (s) as last column of data frame
  s <- 1/hkl_to_reso(hkl$H,hkl$K,hkl$L,a,b,c,aa,bb,gg)
  hkl <- cbind(hkl,data.frame(s=s))

  # Eliminate systematic absences from data frame
  hkl <- deplete_systematic_absences(hkl,csym,set)

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
#' @param csym An object of class "cryst_symm".
#' @param set An integer number corresponding to the specific
#'            setting for the given space group. Default is 1.
#' @return hkl The same data frame acquired from input, depleted
#'             of all systematic absences.
#'
#' @examples
#' # Create a C 2 (monoclinic) space group
#' csym <- cryst_symm("C 2")
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
#' new_hkl <- deplete_systematic_absences(hkl,sym)
#'
#' # Compare first 10 items of original and depleted arrays
#' hkl[1:10,]
#' new_hkl[1:10,]
#'
#' @export
deplete_systematic_absences <- function(hkl,csym,set=1) {
  # Make a m X 3 matrix of original Miller indices
  hkl2 <- as.matrix(hkl[,1:3])

  # Delete incorrect Miller indices if they are
  # systematically-absent
  idx <- sysabs(hkl2,csym,set)
  hkl <- hkl[idx,]

  return(hkl)
}


#' Delete systematic absences (Miller indices)
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
#' @param csym An object of class "cryst_symm".
#' @param set An integer number corresponding to the specific
#'            setting for the given space group. Default is 1.
#' @return idx A vector of integers corresponding to the
#'             position, in the array \code{mhkl}, in which the
#'             Miller indices ARE NOT systematically absent.
#'             The position of systematically-absent reflections
#'             can be found using !idx.
#'
#' @examples
#' # Create a C 2 (monoclinic) space group
#' csym <- cryst_symm("C 2")
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
#' # Index corresponding to valid reflections
#' # (not systematic absences)
#' idx <- sysabs(hkl,csym)
#'
#' # Indices for all reflections
#' fulldx <- 1:length(hkl[,1])
#'
#' # Index corresponding to systematic absences
#' jdx <- fulldx[-idx,]
#'
#' # A couple of systematic absences
#' hkl[jdx[1:2],]
#'
#' @export
sysabs  <- function(hkl,csym,set=1) {
  # Extract symmetry number and setting
  value <- csym$SG
  sym_number <- translate_SG(value,SG_in="xHM",
                             SG_out="number",set=set)
  if (!sym_number$ans) stop(sym_number$msg)
  sym_number <- sym_number$msg
  setting <- set

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
    if (setting == 2) hkl[(hkl[,2]+hkl[,l])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,l])%%2 != 0,4] <- NA
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
    if (setting == 2) hkl[(hkl[,2]+hkl[,l])%%2 != 0,4] <- NA
    if (setting == 3) hkl[(hkl[,1]+hkl[,l])%%2 != 0,4] <- NA
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
