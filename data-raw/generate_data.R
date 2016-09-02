#
# This file is part of the cry package
#

# The file "sysdata.rda" contained under the "R/" directory, includes vital
# information for the whole cry package. This information is hidden to the
# user.
#
# This file can be recreated at any time starting from the 4 following
# files:
#
# atomsf.lib
# syminfo.lib
# elements_list.dat
# amino_list.dat
#
# Simply source this file from the cry project directory:
#
#     source("data-raw/generate_data.R")
#
# and use function generate_sysdata:
#
#     generate_sysdata()
#
# File "sysdata.rda" will be added under the "R" directory.
#
# Information from atomsf.lib is temporarily unloaded. It will be included
# into "sysdata.rda" later.

generate_sysdata <- function()
{
  require(devtools)

  # Symmetry information
  symfile <- file.path("data-raw","syminfo.lib")
  syminfo <- scan(file=symfile,what="character",sep="\n",quiet=TRUE)

  # Chemical elements
  elefile <- file.path("data-raw","elements_list.dat")
  .ATOMS_data.frame <- read.table(file=elefile,header=TRUE,as.is=1:2)

  # Alphabet
  .ALPHABET <- c(LETTERS," ")

  # Amino acids
  aminofile <- file.path("data-raw","amino_list.dat")
  .AMINO_ACIDS_data.frame <- read.table(file=aminofile)

  # Save all data
  use_data(syminfo,.ATOMS_data.frame,.ALPHABET,.AMINO_ACIDS_data.frame,
           internal=TRUE)
}
