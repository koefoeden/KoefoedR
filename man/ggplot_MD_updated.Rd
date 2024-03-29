% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCOP_mRNA_seq.R
\name{ggplot_MD_updated}
\alias{ggplot_MD_updated}
\title{Plot MD-figures for all samples in DGElist.}
\usage{
ggplot_MD_updated(object, samples = "all", ncol = 3)
}
\arguments{
\item{object}{The DGElist object with transcript/gene counts and sample information
encoded in the y$sample object}

\item{samples}{The samples to plot as a character vector}

\item{ncol}{Number of columns in the combined plot as integer vector}
}
\value{
A facetted ggplot
}
\description{
Plot MD-figures for all samples in DGElist.
}
