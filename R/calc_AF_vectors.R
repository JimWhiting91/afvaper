#' Calculate Allele Frequency Change Vectors from a VCF
#'
#' Takes as input a population map (popmap) and a VCF (vcfR object) and returns an m x n matrix, where m = the number of SNPs and n = the number of populations
#'
#' @param vcf A vcf stored as a vcfR object
#' @param popmap A two column data.frame, column 1 lists individuals and column 2 population assignment
#' @param window_size Integer value describing how many SNPs to include per window
#' @param vectors A list object that descibes each vector. Each list element should be a character vector with two strings, each corresponding to a population in the popmap
#' @param normalise Boolean to determine whether allele frequencies should be normalised. Normalisation is essential for vector analysis, but setting normalise=FALSE permits comparison with the raw allele frequencies
#' @param end_cutoff Numeric value describing how many SNPs the last window can have. By default, this equals the window_size, such that the last window must have the same number of SNPs as all windows. This is highly recommended.
#' @param null_perms Integer value describing how many null permutations to run. When set to NULL, function returns observed allele frequency vectors. When set to any integer value, the function returns randomised permutated vectors that are used for estimation of the null distribution.
#' @param n_cores number of cores to run with mclapply. If n_cores=1, runs with lapply instead.
#'
#' @import vcfR
#' @import parallel
#' @import data.table
#' @importFrom stats na.omit
#'
#' @return A list object where each element corresponds to an m (number of SNPs per window) x n (number of populations) matrix of allele frequency change vectors
#' @export
calc_AF_vectors <- function(vcf=NULL,
                            window_size=200,
                            popmap=NULL,
                            vectors=NULL,
                            n_cores=1,
                            normalise=TRUE,
                            end_cutoff=window_size,
                            null_perms=NULL){

  # Set up parallel if needs be
  this.lapply <- if (n_cores>1) { function (...) parallel::mclapply(...,mc.cores=n_cores) } else { lapply }

  # Check the inputs
  if(is.null(popmap)){
    stop("Error, no popmap supplied...")
  }
  if(is.null(vectors)){
    stop("Error, no vector list supplied...")
  }
  if(is.null(vcf)){
    stop("Error, no VCF (vcfR object) supplied...")
  }

  # Check vector pops are in popmap
  if(any(!(unlist(popmap) %in% unique(popmap[,2])))){
    stop("Error, there are populations in the vectors input that are not in the popmap...")
  }

  # Check popmap individuals against the VCF
  if(any(!(colnames(vcf@gt)[2:ncol(vcf@gt)] %in% popmap[,1]))){
    stop("Error, there are individuals in the popmap that are not in the VCF...")
  }

  # Check for only one chr and make SNP windows
  if(length(unique(vcf@fix[,1])) > 1){
    stop("Error, VCF contains more than one chromosome. Please run over individual chromosomes/scaffolds.")
  }

  # Filter away any invariants
  vcf_inds <- colnames(vcf@gt)[2:ncol(vcf@gt)]
  if(!(any(vcf_inds %in% popmap[,1]))){
    vcf_inv <- vcf[,c(1,which(vcf_inds %in% popmap[,1]))]
    vcf_inv <- vcf_inv[is.polymorphic(vcf_inv,na.omit = T),]
    message(paste0(nrow(vcf@fix) - nrow(vcf_inv@fix)," invariants removed"))
    vcf <- vcf_inv
  }

  # Make windows
  winds <- seq(1,nrow(vcf@fix),window_size)
  winds2 <- winds+window_size-1
  winds2[length(winds2)] <- nrow(vcf@fix)

  # Trim any windows that are less than a user-defined cutoff SNPs...
  if(max(winds2)-max(winds) < end_cutoff){
    winds <- winds[1:(length(winds)-1)]
    winds2 <- winds2[1:(length(winds2)-1)]
  }

  # How many windows to estimate
  if(is.null(null_perms)){
    message(paste0("Calculating AF vectors for ",length(winds)," windows with ",window_size," SNPs each"))
  } else {
    message(paste0("Calculating NULL AF vectors for ",null_perms," windows with ",window_size," SNPs each"))
  }

  # Filter for the popmap
  # Do a check here on popmap stability
  if(any(!(popmap[,1] %in% vcf_inds))){
    if(length(popmap[!(popmap[,1] %in% vcf_inds),1]) < 6){
      message(paste0("Warning, popmap inds:",paste(popmap[!(popmap[,1] %in% vcf_inds),1],collapse = ",")," are not in VCF"))
    } else {
      message(paste0("Warning, ",length(popmap[!(popmap[,1] %in% vcf_inds),1])," popmap inds are not in VCF"))
    }
    message("Continuing analysis, but these individuals will be ignored")
  } else {
    message("Popmap check passed, all popmap inds are in VCF")
  }

  # Do a check here on VCF inds stability
  if(any(!(vcf_inds %in% popmap[,1]))){
    if(length(vcf_inds[!(vcf_inds %in% popmap[,1])]) < 6){
      message(paste0("Warning, VCF inds:",paste(vcf_inds[!(vcf_inds %in% popmap[,1])],collapse = ",")," are not in popmap"))
    } else {
      message(paste0("Warning, ",length(vcf_inds[!(vcf_inds %in% popmap[,1])])," VCF inds are not in popmap"))
    }
    message("Continuing analysis, but these individuals will be ignored")
  } else {
    message("VCF check passed, all VCF inds are in popmap")
  }

  # Are we normalising
  if(!(normalise)){
    message("Warning: Allele frequencies will not be normalised, do not use for eigen analysis")
  }

  # Only calculate here if we aren't running null permutations
  if(is.null(null_perms)){

    # Get the allele frequencies
    AF_mat <- vcfR2AF(vcf,popmap,n_cores)

    # Build the full output matrix
    full_vectors <- matrix(nrow=length(vectors),ncol=nrow(AF_mat))
    for(m in 1:nrow(full_vectors)){
      full_vectors[m,] <- AF_mat[,vectors[[m]][2]] - AF_mat[,vectors[[m]][1]]
    }
    colnames(full_vectors) <- paste0(unique(vcf@fix[,1]),"_",vcf@fix[,2])
    rownames(full_vectors) <- names(vectors)


    # Run in parallel over all windows
    window_list <- this.lapply(1:length(winds),function(x){

      # Subset for the window
      sub_vectors <- full_vectors[,winds[x]:winds2[x]]

      # Normalise for correlation matrix, otherwise keep as raw...
      if(normalise){
        for(m in 1:nrow(sub_vectors)){
          sub_vectors[m,] <- sub_vectors[m,]/sqrt(sum(na.omit(sub_vectors[m,])^2))
        }
      }

      # Any remaining NA values, reduce to 0
      sub_vectors[is.na(sub_vectors)] <- 0
      return(sub_vectors)
    },mc.cores=n_cores)

    # Fetch start and end pos as names
    names(window_list) <- sapply(1:length(winds),function(x){
      return(paste0(vcf@fix[,1][1],":",vcf@fix[,2][winds[x]],"-",vcf@fix[,2][winds2[x]]))
    })

  } else {

    # Make null windows
    snp_count <- nrow(vcf@fix)
    null_sites <- (1:snp_count)[1:(snp_count-window_size)]
    null_winds <- sort(sample(x = null_sites,size = null_perms,replace = F))
    null_winds2 <- null_winds + window_size - 1

    window_list <- this.lapply(1:null_perms,function(x){

      # Make a randomised popmap
      null_popmap <- popmap
      null_popmap[,1] <- sample(null_popmap[,1],replace = F)

      # Subset the vcf
      vcf_sub <- vcf[null_winds[x]:null_winds2[x],]

      # Get AFs
      null_AFs <- vcfR2AF(vcf_sub,null_popmap,1)

      # Build the full output matrix
      null_vectors <- matrix(nrow=length(vectors),ncol=nrow(null_AFs))
      for(m in 1:nrow(null_vectors)){
        null_vectors[m,] <- null_AFs[,vectors[[m]][2]] - null_AFs[,vectors[[m]][1]]
      }
      colnames(null_vectors) <- paste0(unique(vcf_sub@fix[,1]),"_",vcf_sub@fix[,2])
      rownames(null_vectors) <- names(vectors)

      # Normalise for correlation matrix, otherwise keep as raw...
      if(normalise){
        for(m in 1:nrow(null_vectors)){
          null_vectors[m,] <- null_vectors[m,]/sqrt(sum(na.omit(null_vectors[m,])^2))
        }
      }

      # Any remaining NA values, reduce to 0
      null_vectors[is.na(null_vectors)] <- 0
      return(null_vectors)
    },mc.cores=n_cores)

    # Fetch start and end pos as names
    names(window_list) <- sapply(1:length(null_winds),function(x){
      return(paste0(vcf@fix[,1][1],":",vcf@fix[,2][null_winds[x]],"-",vcf@fix[,2][null_winds2[x]]))
    })

  }

  # Return
  return(window_list)
}
