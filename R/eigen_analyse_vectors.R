#' Eigen Analysis of Allele Frequency Change Vector Matrix
#'
#' Performs eigen decomposition of a single matrix of normalised allele frequency change vectors. Can be run across the output of calc_AF_vectors() by combining with lapply()
#'
#' @param vector_input A matrix of allele frequency change vectors
#'
#' @return A list object with three elements detailing the results: "eigenvals" contains the eigenvalue distribution, "eigenvecs" contains the eigenvectors, "A_matrix" contains the projected SNP scores for each eigenvector
#' @export
eigen_analyse_vectors <- function(vector_input){
  # Calculate matrix C
  C <- vector_input%*%t(vector_input)  
  # eigen decomposition of C
  eig <- eigen(C)
  # eigenvectors (Q) of C 
  vecs <- eig$vectors
  vals <- eig$values
  a1 <- t(vecs)
  A <- t(vector_input) %*% solve(a1)

  # Tidy up
  rownames(vecs) <- rownames(vector_input)
  colnames(vecs) <- paste0("Eigenvector_",1:ncol(vecs))
  names(vals) <-  paste0("Eigenvector_",1:length(vals))

  # Return
  out <- list(vals,vecs,A)
  names(out) <- c("eigenvals","eigenvecs","A_matrix")
  return(out)
}
