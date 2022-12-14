% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metabolite_database.R
\name{nEntries}
\alias{nEntries}
\alias{nEntries,MetaboliteDatabase-method}
\alias{filterMR}
\alias{filterMR,MetaboliteDatabase-method}
\alias{filterER}
\alias{filterER,MetaboliteDatabase-method}
\alias{filterIP}
\alias{filterIP,MetaboliteDatabase-method}
\alias{filterEntries}
\alias{filterEntries,MetaboliteDatabase-method}
\alias{filterMF}
\alias{filterMF,MetaboliteDatabase-method}
\title{Metabolite database utilities}
\usage{
nEntries(db)

\S4method{nEntries}{MetaboliteDatabase}(db)

filterMR(db, lower, upper)

\S4method{filterMR}{MetaboliteDatabase}(db, lower, upper)

filterER(db, rule)

\S4method{filterER}{MetaboliteDatabase}(db, rule)

filterIP(db, rule)

\S4method{filterIP}{MetaboliteDatabase}(db, rule)

filterEntries(db, IDs)

\S4method{filterEntries}{MetaboliteDatabase}(db, IDs)

filterMF(db, mf)

\S4method{filterMF}{MetaboliteDatabase}(db, mf)
}
\arguments{
\item{db}{S4 object of class \code{MetaboliteDatabase}}

\item{lower}{lower mass boundary}

\item{upper}{upper mass boundary}

\item{rule}{a filtering expression}

\item{IDs}{a numeric vector of entry IDs}

\item{mf}{a molecular formula to filter}
}
\value{
An S4 object of class \code{MetaboliteDatabase}.
}
\description{
Utilities for working with metabolite databases.
}
\examples{
metabolite_database <- metaboliteDB(amino_acids)

## Return the number of database entries
nEntries(metabolite_database)

## Filter database entries
filterEntries(metabolite_database,c(1:5))

## Filter database using a mass range
filterMR(metabolite_database,100,120)

## Filter the database by an element frequency rule
filterER(metabolite_database,C > 2)

## Filter the database by an ionisation product rule
filterIP(metabolite_database,HBA2>0 & Total_Charge==0)

## Filter a database by a molecular formula
filterMF(metabolite_database,"C3H7NO2")
}
