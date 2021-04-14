#' Find Null Cutoffs Based on Randomised Null Permuted Allele Frequency Change Vectors
#'
#' This function takes the output the eigen analysis and returns a matrix of summed values where each eigenvalue is the sum of itself and the eigenvalues that preceed it. This summed value is what is plotted and tested for significance. This function can be run over a list of eigen_analyse_vectors with lapply().
#'
#' @param null_res A list of NULL allele frequency change vector matrices as output by calc_AF_vectors() with null_perms set to a numeric value.
#' @param cutoffs Numeric vector of alpha thresholds for null distribution cutoffs, for e.g. 0.95 = 5%, 0.99 = 1% e.t.c.
#'
#' @importFrom stats quantile
#'
#' @return A matrix of null cutoffs with one row per eigenvalue and one column per significance threshold
#' @export
find_null_cutoff <- function(null_res = NULL,
                             cutoffs = c(0.99)){
  # Do eigen analysis
  eigen_res <- lapply(null_res,eigen_analyse_vectors)

  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)

  # Sum them
  null_sum <- sum_eigenvals(eigen_res)

  # Make an outmatrix...
  out_mat <- matrix(ncol=length(cutoffs),nrow=eigen_count)
  for(i in 1:ncol(out_mat)){
    for(j in 1:nrow(out_mat)){
      out_mat[j,i] <- quantile(null_sum[,j],probs = cutoffs[i])
    }
  }

  # Tidy and return
  colnames(out_mat) <- paste0(cutoffs*100,"%")
  rownames(out_mat) <- paste0("Eigenvector ",1:nrow(out_mat))
  return(out_mat)
}

