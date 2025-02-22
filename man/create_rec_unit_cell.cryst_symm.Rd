% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_S3.R
\name{create_rec_unit_cell.cryst_symm}
\alias{create_rec_unit_cell.cryst_symm}
\title{Reciprocal unit cell from a 'cryst_symm' object}
\usage{
\method{create_rec_unit_cell}{cryst_symm}(ar, ...)
}
\arguments{
\item{ar}{An object of class 'cryst_symm'.}

\item{...}{Additional arguments passed to the
create_rec_unit_cell.}
}
\value{
An object of class "rec_unit_cell". It is a named list
        of length 6 whose last three slots are of class
        'angle'. The cell parameters are calculated from those
        of the corresponding unit cell. The default unit cell
        parameters are a=10, b=20, c=15, alpha=70, beta=80,
        gamma=100. When constrains due to symmetry are
        required, b and c might be equalled to a, alpha, beta
        and gamma might be set to 90, gamma might be set to
        120 and the three angles might be set equal to each
        other.
}
\description{
Method to create an object of class "rec_unit_cell" starting
from an object of class 'cryst_symm'.
}
\details{
The symmetry of a space group imposes constrains on the
parameters of unit cells. For example, the cubic group P 2 3
means that all cell sides have to be equal and all angles
have to be equal to 90 degrees. This function suggests the
appropriate reciprocal cell compatible with the given space
group.
}
\examples{
# Symmetry "C 1 2/c 1"
csym <- cryst_symm("C 1 2/c 1")

# Reciprocal unit_cell
ruc <- create_rec_unit_cell(csym)
print(ruc)

}
