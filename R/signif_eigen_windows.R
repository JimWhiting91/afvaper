#' Identify Windows with Parallel Evolution Evidence Based On Significant Eigenvalues
#'
#' Uses the results from eigen_analyse_vectors() and cutoffs derived from find_null_cutoffs() to return a list of significant windows associated with each eigenvector. Windows are associated with the highest order eigenvector for which their summed eigenvalues exceed the significance threshold.
#' @param eigen_res A list of results from eigen_analyse_vectors(). This can be easily created by combining eigen_analyse_vectors with lapply. Each element of the input list should correspond to the eigen decomposition results for a given genome window.
#' @param cutoffs Numeric vector of alpha thresholds for null distribution cutoffs. These must be equal in length to the number of eigenvalues, such that the first element is the threshold for eigenvalue 1, the second for eigevalue 2 e.t.c.
#'
#' @importFrom stats p.adjust
#'
#' @return A list of character vectors, where each element of the list represents the significant windows associated with each eigenvalue.
#' @export
signif_eigen_windows <- function(eigen_res = NULL,
                                 cutoffs = NULL){

  if(is.numeric(cutoffs) == FALSE){
    stop("Error: cutoffs should be a numeric vector with a cutoff per eigenvector")
  }
  # if(!(is.null(pval_alpha)) & is.null(null_vectors)){
  #   stop("Error: If using p-values, please provide null_vectors for calculation")
  # }

  # # Calculate pvalues if we are using them
  # if(!(is.null(pval_alpha))){
  #   tmp_pvals <- eigen_pvals(eigen_res,null_vectors)
  #   for(i in 1:ncol(tmp_pvals)){
  #     tmp_pvals[,i] <- p.adjust(tmp_pvals[,i],method=p.adjust.method)
  #   }
  # }

  # Fetch sums
  eigen_sums <- sum_eigenvals(eigen_res)

  # Run for each eigenvector
  eigenvec_out <- lapply(1:ncol(eigen_sums),function(vec){

    # Fetch the significant
    # if(is.null(pval_alpha)){
      signif_windows <- rownames(eigen_sums[eigen_sums[,vec] > cutoffs[vec],,drop=FALSE])
    # } else {
    #  signif_windows <- rownames(tmp_pvals[tmp_pvals[,vec] <= pval_alpha,])
    # }
    # Return
    if(length(signif_windows) > 0){
      return(signif_windows)
    } else {
      return(NA)
    }
  })

  # Remove windows that are already significant based on previous eigenvalues...
  for(i in 2:length(eigenvec_out)){
    eigenvec_out[[i]] <- eigenvec_out[[i]][!(eigenvec_out[[i]] %in% unlist(lapply(1:(i-1),function(x){eigenvec_out[[x]]})))]
  }

  # Set names
  names(eigenvec_out) <- paste0("Eigenvector ",1:length(eigenvec_out))
  return(eigenvec_out)
}
