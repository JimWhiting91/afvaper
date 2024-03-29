% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_null_cutoff.R
\name{find_null_cutoff}
\alias{find_null_cutoff}
\title{Find Null Cutoffs Based on Randomised Null Permuted Allele Frequency Change Vectors}
\usage{
find_null_cutoff(null_res = NULL, cutoffs = c(0.99))
}
\arguments{
\item{null_res}{A list of NULL allele frequency change vector matrices as output by calc_AF_vectors() with null_perms set to a numeric value.}

\item{cutoffs}{Numeric vector of alpha thresholds for null distribution cutoffs, for e.g. 0.95 = 5\%, 0.99 = 1\% e.t.c.}
}
\value{
A matrix of null cutoffs with one row per eigenvalue and one column per significance threshold
}
\description{
This function takes the output the eigen analysis and returns a matrix of summed values where each eigenvalue is the sum of itself and the eigenvalues that preceed it. This summed value is what is plotted and tested for significance. This function can be run over a list of eigen_analyse_vectors with lapply().
}
