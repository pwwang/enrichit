% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichit.R
\name{CheckUserList}
\alias{CheckUserList}
\title{Check user list against gene sets}
\usage{
CheckUserList(userlist, genesets, use_matched_only, show_all_unmatched = FALSE)
}
\arguments{
\item{userlist}{Vector of user-provided genes}

\item{genesets}{List of gene sets}

\item{use_matched_only}{Logical, whether to return only matched genes
Otherwise the genes in the user list will be all used.
This will affect the number of genes in the user list when computing the
p-value.}

\item{show_all_unmatched}{Logical, whether to show all unmatched genes. Or instead
show only the first 10 unmatched genes.}
}
\value{
A vector of user-provided genes that are also in the gene sets
}
\description{
Check user list against gene sets
}
\keyword{internal}
