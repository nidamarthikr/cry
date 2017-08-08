#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads SHELX log file file
#'
#' @param filename A character string. The path to a valid shelx log file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return A named list. Each name correspond to a valid field in the shelx log
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/shelxlog/file.log"
#' ltmp <- .readXDS(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export


#filename <- "/Users/ritagiordano/Desktop/Rita/crydir/SHELX/fae_kappa_shelxc.log"
# filename <- "/Users/ritagiordano/Desktop/Rita/crydir/SHELX/fae_kappa_shelxd.log"

filename <- "/Users/ritagiordano/Desktop/Rita/crydir/SHELX/fae_kappa_shelxe_i.log"

## Build the function use .tagRead from reflection_dot.R
# Create a function to extract CCall/CCweak from SHELXE log file
.extract_RE <- function(filename, message = FALSE)
{
  header <- scan(filename, nlines = 3, what = character(), quiet = TRUE)
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
}

#### Maybe this function is not very useful ...
## .exportCSV
#' Export shelx table to csv file.
#'
#' @param filenameDF A data frame. The name of the sheld data sets to plot CCall vs Ccweak
#' @param path Path where the file will be saved.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return plot
#' @examples
#' \dontrun{
#' filename <- shelxd_dataframe
#' ltmp <- .plotSHELXD(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export

.exportCSV <- function(filenameDF, path, message = TRUE) {
  write.csv(filenameDF, file = path, sep = " ", row.names = FALSE)
}


## .plotSHELXC
#' Plot SHELXC log file
#'
#' @param filenameDF A character string. The name of the shelc data frame.
#' @param var The variable to be plotted against the resolution
#' @param ylabel The y label.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return plot
#' @examples
#' \dontrun{
#' filename <- shelxc_dataframe
#' ltmp <- .plotSHELXC(filenameDF)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export
.plotSHELXC <- function(filenameDF, var, message = TRUE) {

  gg <- ggplot(filenameDF, aes(1/(Res)^2, var))
  if(var == filenameDF$d_sig && length(var) >= 1) {
    gp <- gg + geom_point() + geom_line() + theme_bw() +
      xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Delta*F/sig(Delta*F)))
  }
  if(var == filenameDF$Chi_sq && length(var) >= 1){
    gp <- gg + geom_point() + geom_line() + theme_bw() +
      xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Chi^2))
  }
  if(var == filenameDF$I_sig && length(var) >= 1){
    gp <- gg + geom_point() + geom_line() + theme_bw() +
      xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(I/sig(I)))
  }
  if(var == filenameDF$Complete && length(var) >= 1) {
    gp <- gg + geom_point() + geom_line() + theme_bw() +
      xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(Completeness))
  }
  if(var == filenameDF$CC1_2 && length(var) >= 1){
    gp <- gg + geom_point() + geom_line() + theme_bw() +
      xlab(expression(h^2 * (ring(A)^-2))) + ylab(expression(CC_(1/2)))
  }
  # gg + geom_point() + geom_line() + theme_bw() +
  #   xlab(expression(h^2 * (A^-2))) + ylab(expression(Delta*F/sig(Delta*F)))

  return(gp)

}


## .plotSHELXD
#' Plot SHELXD log file
#'
#' @param filenameDF A list containg the SHELXE data frame. The name of the shelxe data sets
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return plot
#' @examples
#' \dontrun{
#' filename <- shelxd_dataframe
#' ltmp <- .plotSHELXD(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export
.plotSHELXD <- function(filenameDF, message = TRUE) {
  library("ggplot2")
  ggplot(filenameDF,aes(CCall,CCweak)) + geom_point() + theme_bw()
}


## .plotSHELXE
#' Plot SHELXE log file
#'
#' @param filenameDF A character string. The name of the sheld data sets to plot CCall vs Ccweak
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return plot
#' @examples
#' \dontrun{
#' filename <- shelxd_dataframe
#' plot <- .plotSHELXD(filename)
#' }
#' @export
.plotSHELXE <- function(filenameDF_i, filenameDF_o, message = TRUE) {

  #if(filenameDF_i == filenameDF_i$CYCLE){
    # create a new column wich will be the total number of cycles
    lcycle_i <- length(shelxE_i$CYCLE[,1])
    lcycle_o <- length(shelxE_o$CYCLE[,1])
    ncycle_i <- seq(1,lcycle_i)
    ncycle_o <- seq(1,lcycle_o)
    CYCLE_i <- cbind(shelxeDF$CYCLE, ncycle_i)
    # Change the name of the new variable
    names(CYCLE_i)[names(CYCLE_i) == 'ncycle_i'] <- 'ncycle'
    CYCLE_o <- cbind(shelxE_o$CYCLE, ncycle_o)
    names(CYCLE_o)[names(CYCLE_o) == 'ncycle_o'] <- 'ncycle'

    # Since the DF can have different length, join the inverted and original DF
    # to create a new data frame. All the column must have the same name.
    df <- rbind(CYCLE_o, CYCLE_i)
    df$dataset <- c(rep("Original", nrow(CYCLE_o)), rep("Inverted", nrow(CYCLE_i)))


    # Plot Contrast vs. Cycle
    gg <- ggplot(df, aes(ncycle, Contrast, color=dataset, group = dataset)) +
          geom_line() + geom_point() + theme_bw() + xlab("Cycle")
    return(gg)
  #}
}

