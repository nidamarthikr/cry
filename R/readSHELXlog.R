#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads and SHELXD log files
#'
#' @param filename A character string. The path to a valid log file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the log header.
#' @return A named list. Each name correspond to a valid field in the log
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/shelxd/shelxd.log"
#' ltmp <- .readSHELXDlog(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export

readSHELXlog <- function(filename)
{
  # data<-readLines(filename)
  header <- scan(filename, nlines = 3, what = character(), quiet = TRUE)
  logfile <- header <- scan(filename, nlines = 3, what = character(), quiet = TRUE)
  # if condition lo search for the right SHELX log file.
  #if(header[3] == grep("SHELXC", header, value = TRUE)) #&& length(header[3] != 0))
  x <- header[3] == grep("SHELXC", header, value = TRUE)
  if (length(x) != 0 && header[3] == grep("SHELXC", header, value = TRUE))
  {
    data <- readLines(filename)
    ## Extract all the row containg the information to be plotted.

    ## Resolution
    res <- grep("Resl.", data, value = TRUE)
    res_1 <- gsub("Resl.|Inf.|-", "", res)
    res_2 <- gsub("[[:blank:]]+", " ", res_1)
    res_split <- strsplit(res_2, split = '[[:blank:]]+')
    res_df <- as.data.frame(res_split, col.names = "Res", stringsAsFactors = FALSE)
    # res_df = res_df[!res_df<=1]

    ## N(data)
    N_data <- grep("N\\(data\\) ", data, value = TRUE)
    N_data_1 <- gsub("N\\(data\\)", "", N_data)
    N_data_split <- strsplit(N_data_1, split = '[[:blank:]]+')
    N_data_df <- as.data.frame(N_data_split, col.names = "N_data", stringsAsFactors = FALSE)
    # N_data_df = N_data_df[!N_data_df<=1]

    ## Chi-sq
    Chi_sq <- grep("Chi-sq ", data, value = TRUE)
    Chi_sq_1 <- gsub("Chi-sq", "", Chi_sq)
    Chi_sq_split <- strsplit(Chi_sq_1, split = '[[:blank:]]+')
    Chi_sq_df <- as.data.frame(Chi_sq_split, col.names = "Chi_sq", stringsAsFactors = FALSE)

    ## <I/sig>
    I_sig <- grep("<I/sig> ", data, value = TRUE)
    I_sig_1 <- gsub("<I/sig>", "", I_sig)
    I_sig_split <- strsplit(I_sig_1, split = '[[:blank:]]+')
    I_sig_df <- as.data.frame(I_sig_split, col.names = "I_sig", stringsAsFactors = FALSE)

    ## %Complete
    Complete <- grep("%Complete ", data, value = TRUE)
    Complete_1 <- gsub("%Complete", "", Complete)
    Complete_split <- strsplit(Complete_1, split = '[[:blank:]]+')
    Complete_df <- as.data.frame(Complete_split, col.names = "Complete", stringsAsFactors = FALSE)

    ## <d\"/sig>
    d_sig <- grep("<d\"/sig> ", data, value = TRUE)
    d_sig_1 <- gsub("<d\"/sig>", "", d_sig[1])
    d_sig_split <- strsplit(d_sig_1, split = '[[:blank:]]+')
    d_sig_df <- as.data.frame(d_sig_split, col.names = "d_sig", stringsAsFactors = FALSE)

    ## CC(1/2)
    cc <- grep("CC\\(1/2\\) ", data, value = TRUE)
    cc_1 <- gsub("CC\\(1/2\\)", "", cc)
    cc_split <- strsplit(cc_1, split = '[[:blank:]]+')
    cc_df <- as.data.frame(cc_split, col.names = "CC1_2", stringsAsFactors = FALSE)

    shelxc_df <- data.frame(res_df, N_data_df, Chi_sq_df, I_sig_df, Complete_df, d_sig_df, cc_df)
    shelxc_df[,c(1,2,3,4,5,6,7)] <- sapply(shelxc_df[,c(1,2,3,4,5,6,7)], as.numeric)
    extract_data <- shelxc_df[!apply(shelxc_df == "", 1, all),]
  }

  x <- header[3] == grep("SHELXD", header, value = TRUE)
  if (length(x) != 0 && header[3] == grep("SHELXD", header, value = TRUE))
  {
    data <- readLines(filename)
    try_data<-grep("Try",data, value = TRUE)
    # Remove punctuation, blanks and CPU word.
    tmp1 <- gsub(",|CPU|/", "", try_data)
    tmp1_2 <- gsub("[[:blank:]]+", " ", tmp1)
    test <- strsplit(tmp1_2, split = '[[:blank:]]+')
    tt <- as.data.frame(test, stringsAsFactors = FALSE)
    tt_df <- as.data.frame(t(tt), row.names = "", stringsAsFactors = FALSE)
    extract_data <-cbind(tt_df$V7,tt_df$V8)
    colnames(extract_data) <- c("CCall","CCweak")
    extract_data <- as.data.frame(extract_data, stringsAsFactors = FALSE)
    # extract_data <-cbind(ccall_ccw[,2],ccall_ccw[,6])
    # colnames(extract_data) <- c("CCall","CCweak")
  }
  x <- header[3] == grep("SHELXE", header, value = TRUE)
  if (length(x) != 0 && header[3] == grep("SHELXE", header, value = TRUE))
  {

    ### Extract cycles of autotracing ###
    data <- readLines(filename)
    extract_data <- list()
    data_reg <- grep("Contrast", data, value = TRUE)
    tmp1 <- gsub("<|>|=|,|for dens.mod.", "", data_reg)
    tmp1_2 <- gsub("[[:blank:]]", " ", tmp1)
    test <- strsplit(tmp1, split = '[[:blank:]]+')
    tt <- as.data.frame(test, stringsAsFactors = FALSE)
    tt_df <- as.data.frame(t(tt), row.names = "", stringsAsFactors = FALSE)
    # Remove unused columns
    tt_df <- tt_df[c("V3", "V5", "V7", "V9")]
    names(tt_df) <- c("wt", "Contrast", "Connect","cycle")
    extr_data <- as.data.frame(tt_df, stringsAsFactors = FALSE)
    extract_data$CYCLE <- extr_data

    ### Extract estimated mean FOM and mapCC as a function of resolution ###
    ## Resolution
    res <- grep("d    inf", data, value = TRUE)
    res_1 <- gsub("d|inf|-", "", res)
    res_2 <- gsub("[[:blank:]]+", " ", res_1)
    res_split <- strsplit(res_2, split = '[[:blank:]]+')
    res_df <- as.data.frame(res_split, col.names = "Res", stringsAsFactors = FALSE)

    ## FOM
    fom <- grep("<FOM>", data, value = TRUE)
    fom1 <- gsub("<FOM>", "", fom)
    fom_split <- strsplit(fom1, split = '[[:blank:]]+')
    fom_df <- as.data.frame(fom_split, col.names = "FOM", stringsAsFactors = FALSE)

    ## MapCC
    mapCC <- grep("<mapCC>", data, value = TRUE)
    mapCC1 <- gsub("<mapCC>", "", mapCC)
    mapCC_split <- strsplit(mapCC1, split = '[[:blank:]]+')
    mapCC_df <- as.data.frame(mapCC_split, col.names = "mapCC", stringsAsFactors = FALSE)

    ## N
    N <- grep("N    ", data, value = TRUE)
    N1 <- gsub("N", "", N)
    N_split <- strsplit(N1, split = '[[:blank:]]+')
    N_df <- as.data.frame(N_split, col.names = "N", stringsAsFactors = FALSE)


    ## Create a data frame with these variables
    shelxe_df <- data.frame(res_df, fom_df, mapCC_df, N_df)
    extract_data2 <- shelxe_df[!apply(shelxe_df == "", 1, all),]

    extract_data$FOM_mapCC <- extract_data2

    ### Extract Density (in map sigma units) at input heavy atom sites ###
    #dt <- readLines(filename)
    # Find the line number associated with the string "Site  "
    nsite <- which(grepl("Site ", data))
    # number os rows containing the first table
    nrow_1 <- (nsite[2] - (nsite[1])-1) -1
    #data_dt <- read.table(filename, skip = nsite[1]-1, nrows = nros_s, header = TRUE, blank.lines.skip = TRUE,
    #                   fill = TRUE)

    # number of line containg the data for the Sites
    #nlineSite <- which(data_dt == "Site")

    Site1 <- read.table(filename, skip = nsite[1]-1, nrows = nrow_1,
                        header = TRUE, blank.lines.skip = TRUE)


    extract_data$Site1 <- Site1
    # Extract the second table containg the sites
    #data2 <- read.table(filename, skip = nsite[2], blank.lines.skip = TRUE,
    #                   fill = TRUE)

    #nlineBest <- which(data2 == "Best")
    nbest <- which(grepl("Best", data))
    nrow_2 <- (nbest - nsite[2] -1)

    Site2 <- read.table(data, skip = nsite[2] -1 , nrows = nrow_2 - 1, header = TRUE,
                        blank.lines.skip = TRUE)

    extract_data$Site2 <- Site2
    ## Create a list which will contains all data frame extrated from shelxe.
    #extract_data <- list(extr_data, extract_data2, Site1, Site2)

  }

  return(extract_data)

  # Convert a character data frame ti a numeric one.

  ## Specifi case for shelxe, becasue there are four datasets
  if (header[3] == grep("SHELXE", header, value = TRUE)) {
    logfile2 <- list()
    CYCLE <- logfile$CYCLE
    logfile2$CYCLE <- as.data.frame(sapply(CYCLE,as.numeric))
    FOM_mapCC <- logfile$FOM_mapCC
    logfile2$FOM_mapCC <- as.data.frame(sapply(FOM_mapCC,as.numeric))
    Site1 <- logfile$Site1
    logfile2$Site1 <- as.data.frame(sapply(Site1,as.numeric))
    Site2 <- logfile$Site2
    logfile2$Site2 <- as.data.frame(sapply(Site2,as.numeric))
  }

  else logfile2 <- as.data.frame(sapply(logfile, as.numeric))
  return(logfile2)
}

