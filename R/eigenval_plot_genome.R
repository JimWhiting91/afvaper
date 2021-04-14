#' Plot Windowed Eigenvalues For A List of Chromosomes
#'
#' Takes as input a list of lists containing results from eigen_analyse_vectors() and plots all results. The chromosome positions are derived from the names() variable of the eigen_res input, which should be in the standard genome interval format of chrX:X-Y.
#' @param eigen_res_list A list of lists containing results from eigen_analyse_vectors(). The length() of the eigen_res_list input should be equal to the number of chromosomes, such that each element of the input is itself a list derived from eigen_analyse_vectors()
#' @param cutoffs Numeric vector of alpha thresholds for null distribution cutoffs. These must be equal in length to the number of eigenvalues, such that the first element is the threshold for eigenvalue 1, the second for eigevalue 2 e.t.c.
#' @param null_vectors A list of NULL allele frequency change vector matrices as output by calc_AF_vectors() with null_perms set to a numeric value. This is only necessary if wanting to calculate and plot p.values.
#' @param plot.pvalues Boolean flag as to whether or not to plot p.values instead of summed eigenvalues.
#'
#' @import ggplot2
#' @import data.table
#'
#' @return A list of ggplot figures, where each element corresponds to figure for each eigenvector.
#' @export
eigenval_plot_genome <- function(eigen_res_list=NULL,
                                 cutoffs=NULL,
                                 null_vectors=NULL,
                                 plot.pvalues=F){

  # Fetch chr to plot
  plot_chr <- unlist(lapply(eigen_res_list,function(x){return(names(x))}))
  plot_chr <- strsplit(plot_chr,":")
  plot_chr <- unique(sapply(plot_chr,function(x){return(x[1])}))

  # Fetch all the pos...
  total_plot_dd <- data.frame(rbindlist(lapply(plot_chr,function(chr){

    # Fetch eigen_res
    eigen_res <- eigen_res_list[[chr]]

    # Make plot_dd
    # Count eigenvectors
    eigen_count <- length(eigen_res[[1]]$eigenvals)

    # Find window start and ends...
    pos_string <- sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][2])})
    chr_string <- as.character(sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][1])}))
    start <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][1])}))
    end <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][2])}))
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
    return(plot_dd)
  })))

  # Tidy up
  colnames(total_plot_dd)[3:ncol(total_plot_dd)] <- paste0("Eigenvector ",1:(ncol(total_plot_dd)-2))

  # Factorise according to input vector
  total_plot_dd$chr_F <- factor(total_plot_dd$chr,levels=chrs)

  # Get pvalues if that is what we are plotting...
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
      obs_eig <- total_plot_dd[,eig]
      null_vec <- null_mat[,eig]
      pvals <- sapply(obs_eig,function(x){return(length(null_vec[null_vec > x])+1)}) / (length(null_vec)+1)
      total_plot_dd[,eig] <- -log10(pvals)
    }

  }


  # Visualise each eigenval along the chr
  plots <- lapply(grep("Eigen",colnames(total_plot_dd),value=T),function(val){

    p1 <- ggplot(total_plot_dd,aes(x=pos,y=total_plot_dd[,val]))+
      geom_step()+
      theme_minimal()+
      ggtitle(val)+
      theme(axis.text=element_text(size=16),
            axis.title = element_text(size=18),
            title = element_text(size=20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing.x=unit(0.1, "lines"),
            strip.text = element_text(angle=90,hjust=1,size=14))+
      facet_wrap(~chr_F,nrow = 1,strip.position = "bottom",scales="free_x")

    if(plot.pvalues){
      p1 <- p1 + labs(y=expression(-log[10](p)),x="Genome Position")
    } else {
      p1 <- p1 + labs(y=gsub("vector","value",val),x="Genome Position")
    }

    if(!(is.null(cutoffs)) & !(plot.pvalues)){
      p1 <- p1 + geom_hline(yintercept=cutoffs[val],colour="red2")
    }
    return(p1)
  })
  names(plots) <- grep("Eigen",colnames(total_plot_dd),value=T)
  return(plots)
}

