#' Merge AF-vapeR Eigen Decomposition Results From Multiple Chromosomes
#'
#' Takes a list of multiple
#' @param eigen_res_list A list of results from eigen_analyse_vectors(). This can be easily created by combining eigen_analyse_vectors with lapply. Each element of the input list should correspond to the eigen decomposition results for a given genome window.
#'
#' @return Single list containing eigen analysis results for all chromosomes/scaffolds
#' @export
merge_eigen_res <- function(eigen_res_list){
  return(unlist(eigen_res_list,recursive = FALSE))
}
