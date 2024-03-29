% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptors.R
\name{chemicalDescriptors}
\alias{chemicalDescriptors}
\title{Calculate chemical descriptors}
\usage{
chemicalDescriptors(SMILES)
}
\arguments{
\item{SMILES}{a character vector of valid SMILES}
}
\value{
A tibble containing the chemical descriptor information for the input SMILES.
}
\description{
Calculate chemical descriptors from SMILES.
}
\examples{
chemicalDescriptors(amino_acids$SMILES)
}
