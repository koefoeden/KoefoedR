% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCOP_mRNA_seq.R
\name{get_MD_data_for_all_samples}
\alias{get_MD_data_for_all_samples}
\title{Get MD data for all samples in long tibble}
\usage{
get_MD_data_for_all_samples(object)
}
\arguments{
\item{object}{The DGElist object with transcript/gene counts and sample information
encoded in the y$sample object}
}
\value{
A long tibble with three columns: sample label, mean and diff values
}
\description{
Get MD data for all samples in long tibble
}
