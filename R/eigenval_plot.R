#' Plot Windowed Eigenvalues Along a Chromosome
#'
#' Takes as input a list of results from eigen_analyse_vectors() and plots all results along a chromosome. The chromosome positions are derived from the names() variable of the eigen_res input, which should be in the standard genome interval format of chrX:X-Y
#' @param eigen_res A list of results from eigen_analyse_vectors(). This can be easily created by combining eigen_analyse_vectors with lapply. Each element of the input list should correspond to the eigen decomposition results for a given genome window.
#' @param cutoffs Numeric vector of alpha thresholds for null distribution cutoffs. These must be equal in length to the number of eigenvalues, such that the first element is the threshold for eigenvalue 1, the second for eigevalue 2 e.t.c.
#' @param null_vectors A list of NULL allele frequency change vector matrices as output by calc_AF_vectors() with null_perms set to a numeric value. This is only necessary if wanting to calculate and plot p.values.
#' @param plot.pvalues Boolean flag as to whether or not to plot p.values instead of summed eigenvalues.
#'
#' @import ggplot2
#' @import data.table
#'
#' @return A list of ggplot figures, where each element corresponds to figure for each eigenvector.
#' @export
eigenval_plot <- function(eigen_res=NULL,
                          cutoffs=NULL,
                          null_vectors=NULL,
                          plot.pvalues=FALSE){

  # Check inputs
  if(is.null(eigen_res)){
    stop("Error, no results supplied from eigen_analyse_vectors")
  }
  if(names(eigen_res[[1]]) != c("eigenvals","eigenvecs","A_matrix")){
    stop("Error, input must be results output by eigen_analyse_vectors()")
  }
  if(plot.pvalues & is.null(null_vectors)){
    stop("Error, if plotting p.values you must provide a null vectors input")
  }
  if(names(eigen_res[[1]]) != c("eigenvals","eigenvecs","A_matrix")){
    stop("Error, input must be results output by eigen_analyse_vectors()")
  }

  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)

  # Find window start and ends...
  pos_string <- sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][2])})
  chr_string <- as.character(sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][1])}))
  start <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][1])}))
  plot_dd <- data.frame(chr = chr_string,
                        pos = start)

  # Add in eigenvals
  eigenval_sums <- matrix(ncol=length(eigen_res[[1]]$eigenvals),nrow=nrow(plot_dd))
  for(i in 1:ncol(eigenval_sums)){
    eigenval_sums[,i] <- as.numeric(sapply(eigen_res,function(j){return(j$eigenvals[i])}))
  }
  colnames(eigenval_sums) <- paste0("Eigenvector ",1:ncol(eigenval_sums))

  # Merge and sum
  plot_dd <- cbind(plot_dd,eigenval_sums)
  for(i in 2:eigen_count){
    plot_dd[,paste0("Eigenvector ",i)] <- rowSums(plot_dd[,c(paste0("Eigenvector ",i-1),
                                                             paste0("Eigenvector ",i))])
  }

  # Get pvals if so
  if(plot.pvalues){

    # Get eigen res for null vectors
    null_res <- lapply(null_vectors,eigen_analyse_vectors)

    # Organise nulls to matrix
    null_mat <- matrix(nrow=length(null_res),ncol=length(null_res[[1]]$eigenvals))
    for(i in 1:ncol(null_mat)){
      null_mat[,i] <- sapply(null_res,function(x){x[[1]][[i]]})
    }
    colnames(null_mat) <- paste0("Eigenvector ",1:ncol(null_mat))

    # Sum over
    for(i in 2:ncol(null_mat)){
      null_mat[,i] <- rowSums(null_mat[,(i-1):i])
    }

    for(eig in paste0("Eigenvector ",1:ncol(null_mat))){
      obs_eig <- plot_dd[,eig]
      null_vec <- null_mat[,eig]
      pvals <- sapply(obs_eig,function(x){return(length(null_vec[null_vec > x])+1)}) / (length(null_vec)+1)
      plot_dd[,eig] <- -log10(pvals)
    }

  }

  # Visualise each eigenval along the chr
  plots <- lapply(colnames(eigenval_sums),function(val){

    p1 <- ggplot(plot_dd,aes(x=pos,y=plot_dd[,val]))+
      geom_step()+
      theme_bw()+
      ggtitle(val)+
      theme(axis.text=element_text(size=16),
            axis.title = element_text(size=18),
            title = element_text(size=20))+
      scale_x_continuous(breaks=seq(0,max(plot_dd$pos),2000000),
                         labels=seq(0,max(plot_dd$pos),2000000)/1000000)

    if(plot.pvalues){
      p1 <- p1 + labs(y=expression(-log[10](p)),x="Chr Position (Mb)")
    } else {
      p1 <- p1 + labs(y=gsub("vector","value",val),x="Chr Position (Mb)")
    }

    if(!(is.null(cutoffs)) & !plot.pvalues){
      p1 <- p1 + geom_hline(yintercept=cutoffs[val],colour="red2")
    }

    return(p1)
  })
  names(plots) <- colnames(plot_dd)[3:ncol(plot_dd)]
  return(plots)
}