% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/test.R
\name{deseqTest}
\alias{deseqTest}
\title{DESeq Test}
\usage{
deseqTest(count.data, group1, group2, filter)
}
\arguments{
\item{count.data}{Normalised CountDataSet object to test}

\item{group1}{First group to test, the reference or control}

\item{group2}{Second group to test, the treatment}

\item{filter}{Boolean, whether to apply HTSFilter}
}
\value{
dataframe containing test results
}
\description{
Testa a normalised CountDataSet object using DESeq
}

