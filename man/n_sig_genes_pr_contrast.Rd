% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCOP_mRNA_seq.R
\name{n_sig_genes_pr_contrast}
\alias{n_sig_genes_pr_contrast}
\title{Get number of significant genes/GO-terms per contrast}
\usage{
n_sig_genes_pr_contrast(results_DF_list)
}
\arguments{
\item{results_DF_list}{An list of either DE genes or GO-terms results,
with an element for each contrast.}
}
\value{
A facetted ggplot
}
\description{
Get number of significant genes/GO-terms per contrast
}
