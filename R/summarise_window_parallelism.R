#' Summarise Parallel Evolution in Windows with Significant Eigenvalues
#'
#' Produces a summary table describing how replicate population pairs load on relevant eigenvectors. This information allows you to make sense of what patterns are driving large eigenvalues. For each significant window, the table summarises how many replicates load (above a user-defined threshold) onto the relevant eigenvector, and whether they load in the same direction (parallel) or opposite direction (anti-parallel). For significant windows on eigenvectors 2+, the table includes information on all higher order eigenvectors (e.g. for eigenvector 2, there is a row for both eigenvector 1 and 2). This is because the analysis is focussed on the sum of eigenvalues, therefore all higher order eigenvalues are relevant to interpretation.
#' @param window_id A character string corresponding to a focal chromosomal window. This window must be in the output of names(eigen_res). Significant windows can be found using the function signif_eigen_windows().
#' @param eigen_res A list of results from eigen_analyse_vectors(). This can be easily created by combining eigen_analyse_vectors with lapply. Each element of the input list should correspond to the eigen decomposition results for a given genome window.
#' @param loading_cutoff Numeric value corresponding to the loading threshold that replicate population pairs need to exceed to be considered as associated with a given eigenvector. This threshold applies to positive or negative loadings as absolute loadings are considered when comparing to the threshold. The threshold can be considered as the correlation coefficient between the eigenvector scores for a given replicate population pair, and its observed allele frequency change values. As such, higher loading thresholds correspond to a stronger agreement, and closer association, between a replicate population pair and the relevent eigenvector.
#' @param eigenvector Integer value for the focal eigenvector to be considered, e.g. 1 = eigenvector 1
#'
#' @return A data.frame table summarising parallel evolution
#' @export
summarise_window_parallelism <- function(window_id,
                                         eigen_res,
                                         loading_cutoff = 0.3,
                                         eigenvector = 1){
  # Make the final_out table
  final_out <- data.frame(window_id=NA,
                          eigenvector=NA,
                          parallel_lineages=NA,
                          parallel_pops=NA,
                          antiparallel_pops=NA)

  # Loop over all eigenvectors up to focal
  for(eigN in 1:eigenvector){

    # Fetch loadings for eigenvectors
    loadings <- eigen_res[[window_id]]$eigenvecs[,eigN]

    # Output
    out <- data.frame(window_id=window_id,
                      eigenvector=paste0("Eig",eigN),
                      parallel_lineages = length(loadings[abs(loadings) > loading_cutoff]))

    # Catch cases where all parallel but negative
    if(length(loadings[loadings < 0]) == length(loadings)){
      loadings <- -1 * loadings
    }

    # Now describe parallel/antiparallel
    out$parallel_pops <- paste(names(loadings[loadings >= loading_cutoff]),collapse = ",")
    out$antiparallel_pops <- paste(names(loadings[loadings <= (-1*loading_cutoff)]),collapse = ",")

    # Catch empty parallel_pops columns
    if(length(loadings[loadings >= loading_cutoff]) == 0){
      out$parallel_pops <- out$antiparallel_pops
      out$antiparallel_pops <- ""
    }

    # Rbind together
    final_out <- rbind(final_out,out)
  }

  # Tidy and return
  final_out <- na.omit(final_out)
  rownames(final_out) <- 1:nrow(final_out)

  return(na.omit(final_out))
}
