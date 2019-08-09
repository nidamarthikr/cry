#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Reads XDS CORRECT file, check if there are alien
#' reflections and create a REMOVE.HKL.
#'
#' @param filename A character string. The path to a valid xds file.
#' @param messages A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the xds header.
#' @param alien A logical variable, which the default value is TRUE. This variable
#' @return A named list. Each name correspond to a valid field in the xds
#'    header.
#' @description ISIGMA will cut the data following the wished I over sigma
#' @examples
#' \dontrun{
#' filename <- "path/to/a/valid/HKL/CORRECT.LP"
#' correct <- readCORRECT(filename, alien = TRUE)
#' }
#' @export


readCORRECT <- function(filename, message = TRUE, alien = TRUE, path)
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
