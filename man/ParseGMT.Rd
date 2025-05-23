% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmt.R
\name{ParseGMT}
\alias{ParseGMT}
\title{Read GMT file}
\usage{
ParseGMT(gmtfile, swap_name_desc_if_needed = TRUE)
}
\arguments{
\item{gmtfile}{Path to the GMT file}

\item{swap_name_desc_if_needed}{Logical, whether to swap name and description fields.
They will be swapped only if:
\itemize{
\item \code{swap_name_desc_if_needed} is \code{TRUE}; and
\item The descriptions are not empty; and
\item The descriptions are shorter than the names; and
\item The descriptions are not ID-like (i.e., hsa00001, or 123456).
}}
}
\value{
A list of gene sets, the names of the list are the gene set names,
and each element is a vector of gene names.
}
\description{
Read GMT file
}
\examples{
\donttest{
if (FALSE) {
# Example GMT file content
gmtfile <- system.file("extdata", "KEGG_2021_Human.gmt", package = "enrichit")
# Read the GMT file
ParseGMT(gmtfile)
}
}
}
