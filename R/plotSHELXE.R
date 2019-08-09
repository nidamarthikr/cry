#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Plot SHELXD log files
#'
#' @param filenameDF_i A data frame containing the inverted hand.
#' @param filenameDF_o A data frame containing the original hand.
#' @return A graphics containg the solution founded by SHELXD.
#'
#' @examples
#' \dontrun{
#' plot <- plotSHELXE(filenameDF_i, filenameDF_o, message = TRUE)
#' plot
#' }
#' @export

plotSHELXE <-function(filenameDF_i, filenameDF_o, message = TRUE)
{
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
}
