#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads XDS file file
#'
#' @param filename A character string. The path to a valid shelx log file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the shelx log header.
#' @return A named list. Each name correspond to a valid field in the shelx log
#'    header.
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/xds/file.HKL"
#' ltmp <- .readXDS(filename)
#' print(names(ltmp))
#' print(ltmp$CELL)
#' print(ltmp$SYMM)
#' }
#' @export

#filename = "/Users/ritagiordano/Desktop/Rita/crydir/autoprocessing_FAE_w1_run3_1/grenades_fastproc/XDS_ASCII.HKL"
# Extract variable from XDS or XSCALE header
.extract_RE_xds <- function(filename, message = TRUE)
{

  f_xds<-file(filename,open="r")
  # Determine the number of row present in the file.
  nrow_xds <- length(f_xds)
  rl_xds<-readLines(f_xds)
  xds_data <- list()
  # if(variable == "RESOLUTION_RANGE")
  # {
  ## Resolution range
  tmp1 <- grep("INCLUDE_RESOLUTION_RANGE=", rl_xds, value = TRUE)
  tmp2 <- gsub("!INCLUDE_RESOLUTION_RANGE=", "", tmp1)
  tmp3 <- strsplit(tmp2, split = '[[:blank:]]+')
  xds_data$MAX_RESOLUTION <- tmp3[[1]][3]
  xds_data$MIN_RESOLUTION <- tmp3[[1]][2]
  #xds_data$RESOLUTION_RANGE <- c(as.integer(tmp3[[1]][3]),
  #                              as.integer(tmp3[[1]][2]))
  # max_res <- tmp3[[1]][3]
  # min_res <- tmp3[[1]][2]
  #xds_data <- list(max_resolution = max_res, min_resolution = min_res)
  #var <- paste("max resolution =",tmp3[[1]][3], ";", "min resolution =", tmp3[[1]][2])

  ## Oscillation range
  tmp_OR <- grep("!OSCILLATION_RANGE=", rl_xds, value = TRUE)
  tmp2_OR <- gsub("!OSCILLATION_RANGE=", "", tmp_OR)
  tmp3_OR <- strsplit(tmp2_OR, split = '[[:blank:]]+')
  xds_data$OSCILLATION_RANGE <- tmp3_OR[[1]][2]

  ## SPACE_GROUP_NUMBER
  tmp_SG <- grep("!SPACE_GROUP_NUMBER=", rl_xds, value = TRUE)
  tmp2_SG <- gsub("!SPACE_GROUP_NUMBER=", "", tmp_SG)
  tmp3_SG <- strsplit(tmp2_SG, split = '[[:blank:]]+')
  xds_data$SPACE_GROUP_NUMBER <- tmp3_SG[[1]][2]

  ## !UNIT_CELL_CONSTANTS
  tmp_UC <- grep("!UNIT_CELL_CONSTANTS=", rl_xds, value = TRUE)
  tmp2_UC <- gsub("!UNIT_CELL_CONSTANTS=", "", tmp_UC)
  tmp3_UC <- strsplit(tmp2_UC, split = '[[:blank:]]+')
  xds_data$UNIT_CELL_CONSTANTS <- c(tmp3_UC[[1]][2], tmp3_UC[[1]][3], tmp3_UC[[1]][4],
                                    tmp3_UC[[1]][5], tmp3_UC[[1]][6], tmp3_UC[[1]][7])

  ## !UNIT_CELL_A-AXIS
  tmp_UCA <- grep("!UNIT_CELL_A-AXIS=", rl_xds, value = TRUE)
  tmp2_UCA <- gsub("!UNIT_CELL_A-AXIS=", "", tmp_UCA)
  tmp3_UCA <- strsplit(tmp2_UCA, split = '[[:blank:]]+')
  xds_data$UNIT_CELL_A_AXIS <- c(tmp3_UCA[[1]][2], tmp3_UCA[[1]][3], tmp3_UCA[[1]][4])

  ## !UNIT_CELL_B-AXIS
  tmp_UCB <- grep("!UNIT_CELL_B-AXIS=", rl_xds, value = TRUE)
  tmp2_UCB <- gsub("!UNIT_CELL_B-AXIS=", "", tmp_UCB)
  tmp3_UCB <- strsplit(tmp2_UCB, split = '[[:blank:]]+')
  xds_data$UNIT_CELL_B_AXIS <- c(tmp3_UCB[[1]][2], tmp3_UCB[[1]][3], tmp3_UCB[[1]][4])

  ## !UNIT_CELL_C-AXIS
  tmp_UCC <- grep("!UNIT_CELL_C-AXIS=", rl_xds, value = TRUE)
  tmp2_UCC <- gsub("!UNIT_CELL_C-AXIS=", "", tmp_UCC)
  tmp3_UCC <- strsplit(tmp2_UCC, split = '[[:blank:]]+')
  xds_data$UNIT_CELL_C_AXIS <- c(tmp3_UCC[[1]][2], tmp3_UCC[[1]][3], tmp3_UCC[[1]][4])

  ## !X-RAY_WAVELENGTH=
  tmp_UXW <- grep("!X-RAY_WAVELENGTH=", rl_xds, value = TRUE)
  tmp2_UXW <- gsub("!X-RAY_WAVELENGTH=", "", tmp_UXW)
  tmp3_UXW <- strsplit(tmp2_UXW, split = '[[:blank:]]+')
  xds_data$X_RAY_WAVELENGTS <- tmp3_UXW[[1]][2]

  ## !DETECTOR_DISTANCE
  tmp_UDD <- grep("!DETECTOR_DISTANCE=", rl_xds, value = TRUE)
  tmp2_UDD <- gsub("!DETECTOR_DISTANCE=", "", tmp_UDD)
  tmp3_UDD <- strsplit(tmp2_UDD, split = '[[:blank:]]+')
  xds_data$DETECTOR_DISTANCE <- tmp3_UDD[[1]][2]

  #data <- list(xds_data = xds_data)
  return(xds_data)
}



#filename = "/Users/ritagiordano/Desktop/Rita/crydir/autoprocessing_FAE_w1_run3_1/XDSAPP/results/xa_w1_run3_anom_CORRECT.LP"
.extract_correct <- function(filename, message = TRUE, alien = TRUE, path)
{
  f_correct<-file(filename,open="r")
  # Determine the number of row present in the file.
  nrow_correct <- length(f_correct)
  rl_correct<-readLines(f_correct)
  reg = "\\d+"
  correct_digit <- grep(reg, rl_correct, value = TRUE)
  correct_digitDF <- data.frame(correct_digit, stringsAsFactors = FALSE)

  correct_data <- list()
  if (alien == TRUE)
  {
    alien <- grep("alien", rl_correct, value = TRUE)
    keeps <- c(".", "+")
    tmp <-gsub(paste0(".*?($|'|", paste(paste0("\\",keeps), collapse = "|"),
                "|[^[:punct:]]).*?"), "\\1", alien)

    #tmp <- gsub("[[:blank:]]+", " ", alien)
    test <- strsplit(tmp, split = '[[:blank:]]+')
    tt <- as.data.frame(test, stringsAsFactors = FALSE)
    tt_df <- as.data.frame(t(tt), row.names = "", stringsAsFactors = FALSE)
    tt_df$V1 <- NULL
    colnames(tt_df)<-c("h","k","l","RES", "Z", "Intensity","Sigma", "alien")
    alien_df <- data.frame(lapply(tt_df[1:7], as.numeric), "alien" = tt_df$alien,
                           stringsAsFactors = FALSE)
    # Write REMOVE.HKL
    write.table(alien_df, file = "/Users/ritagiordano/Desktop/crydir/REMOVE.HKL",
                sep ="\t", append = FALSE, row.names = FALSE, col.names = FALSE)
  }


  }

#################################################################################################
# rxstring = "RESOLUTION\\s+RANGE\\s+I/Sigma\\s+Chi\\^2\\s+R-FACTOR\\s+R-FACTOR\\s+NUMBER\\s+ACCEPTED\\s+REJECTED\\s*\\n\\s*observed\\s
#   \\s*\\n\\s*((?:\\n\\s*-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*\\s+\\d+\\s+\\d+\\s+\\d+\\s*)+)"
# regex = "RESOLUTION RANGE  I\\/Sigma  Chi\\^2  R\\-FACTOR  R\\-FACTOR  NUMBER ACCEPTED REJECTED\n"
# reg = "\\d+\\s\\d+\\s+\\d+\\s+\\d+"
#
# tex <- readLines("/Users/ritagiordano/Desktop/Rita/crydir/test.txt")
#
# tblLines<-readLines("/Users/ritagiordano/Desktop/Rita/crydir/test.txt")
# tblDF<-data.frame(tblLines)
# tblDF$isHeader<-substr(tblDF$tblLines,1,1) =="R"
# tblDF$tblNumber<-cumsum(tblDF$isHeader)
# tblDF$rowToRemove<-substr(tblDF$tblLines,1,1) %in% c("R","o")
#
# myTbls<-list()
# for(i in unique(tblDF$tblNumber)){
#   myTbls[[i]]<-read.table(
#     text=paste(tblDF$tblLines[tblDF$tblNumber==i&!tblDF$rowToRemove],collapse="\n")
#   )
# }
#
