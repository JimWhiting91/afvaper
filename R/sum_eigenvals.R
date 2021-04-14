#' Sum Eigenvalues Across eigen_analyse_vector Results
#'
#' This function takes the output the eigen analysis and returns a matrix of summed values where each eigenvalue is the sum of itself and the eigenvalues that preceed it. This summed value is what is plotted and tested for significance. This function can be run over a list of eigen_analyse_vectors with lapply().
#'
#' @param eigen_res A single matrix eigen analysis results as output by eigen_analyse_vectors()
#'
#' @return A matrix of summed eigenvalues
#' @export
sum_eigenvals <- function(eigen_res){

  # Check the input
  if(names(eigen_res[[1]]) != c("eigenvals","eigenvecs","A_matrix")){
    stop("Error, input must be results output by eigen_analyse_vectors()")
  }

  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)

  # Set up an out matrix
  out_mat <- matrix(nrow=length(eigen_res),ncol=eigen_count)
  rownames(out_mat) <- names(eigen_res)

  # Fill the mat
  for(i in 1:eigen_count){
    out_mat[,i] <- sapply(eigen_res,function(x){return(x$eigenvals[i])})
  }

  # Sum the mat
  for(i in 2:eigen_count){
    out_mat[,i] <- rowSums(out_mat[,c(i-1,i)])
  }

  # Clean
  colnames(out_mat) <- paste0("Eigenvalue_",1:eigen_count)

  # Return
  return(out_mat)
}
