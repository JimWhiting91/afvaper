#' Calculate associated p-values for AF-vapeR eigenvalues
#'
#' Takes a list containing results from eigen_analyse_vectors() and calculates p-values based on null vectors.
#' @param eigen_res A list of results from eigen_analyse_vectors(). This can be easily created by combining eigen_analyse_vectors with lapply. Each element of the input list should correspond to the eigen decomposition results for a given genome window.
#' @param null_vectors A list of NULL allele frequency change vector matrices as output by calc_AF_vectors() with null_perms set to a numeric value. This is only necessary if wanting to calculate and plot p.values.
#'
#' @import qvalue
#' 
#' @return Matrix of p-values for each window (rows) for each eigenvector (columns)
#' @export
eigen_pvals <- function(eigen_res,null_vectors){

  # Get eigen res for null vectors
  null_res <- lapply(null_vectors,eigen_analyse_vectors)

  # Organise nulls to matrix
  null_mat <- matrix(nrow=length(null_res),ncol=length(null_res[[1]]$eigenvals))
  for(i in 1:ncol(null_mat)){
    null_mat[,i] <- sapply(null_res,function(x){x[[1]][[i]]})
  }
  
  # Sum over
  for(i in 2:ncol(null_mat)){
    null_mat[,i] <- rowSums(null_mat[,(i-1):i])
  }

  # Organise obs to matrix
  obs_mat <- matrix(nrow=length(eigen_res),ncol=length(eigen_res[[1]]$eigenvals))
  for(i in 1:ncol(obs_mat)){
    obs_mat[,i] <- sapply(eigen_res,function(x){x[[1]][[i]]})
  }
  
  # Sum over
  for(i in 2:ncol(obs_mat)){
    obs_mat[,i] <- rowSums(obs_mat[,(i-1):i])
  }
  
  # Get p-vals
  p_vals <- matrix(ncol=ncol(obs_mat),nrow=nrow(obs_mat))
  for(i in 1:ncol(p_vals)){
    p_vals[,i] <- qvalue:::empPvals(stat = obs_mat[,i],stat0 = null_mat[,i])
  }

  colnames(p_vals) <- paste0("Eigenvalue_",1:ncol(p_vals))
  rownames(p_vals) <- names(eigen_res)
  return(p_vals)
}
