% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert}
\alias{convert}
\title{Convert chemical structure notation}
\usage{
convert(
  input,
  input_type = c("SMILES", "INCHI"),
  output_type = c("INCHI", "INCHIKEY", "SMILES")
)
}
\arguments{
\item{input}{a valid SMILE or InChI}

\item{input_type}{either \code{"SMILES"} or \code{"INCHI"}, denoting the input type}

\item{output_type}{either \code{"SMILES"}, \code{"INCHI"} or \code{"INCHIKEY"}, denoting the output type}
}
\value{
The converted chemical structure.
}
\description{
Convert between SMILES and InChI and to InchiKey chemical structure notations.
}
\details{
This functionality is not currently supported on Windows.
}
\examples{
if (Sys.info()["sysname"] != 'Windows'){
  convert("C[C@H](C(=O)O)N",'SMILES','INCHI')
}
}
