#' Calculate Allele Frequencies from a VCF
#'
#' Takes as input a population map (popmap) and a VCF (vcfR object) and returns an m x n matrix, where m = the number of SNPs and n = the number of populations
#'
#' @param vcf a vcf stored as a vcfR object
#' @param popmap a two column data.frame, column 1 lists individuals and column 2 population assignment
#' @param n_cores number of cores to run with mclapply. If n_cores=1, runs with lapply instead.
#' @param allele character vector of either "REF" or "ALT" for which allele frequency to return
#'
#' @import vcfR
#' @import parallel
#'
#' @return A matrix of population-specific allele frequencies
#' @export
vcfR2AF <- function(vcf,
                    popmap,
                    n_cores=1,
                    allele="REF"){

  # Set up parallel if needs be
  this.lapply <- if (n_cores>1) { function (...) parallel::mclapply(...,mc.cores=n_cores) } else { lapply }

  # First calculate all population allele frequencies
  # Extract gt manually from VCF
  gt <- extract.gt(vcf)

  # Convert to numbers using the appropriate separator
  geno_format <- strsplit(gt[1,1],"")[[1]][2]
  gt[gt == paste0("0",geno_format,"0")] <- 0
  gt[gt == paste0("0",geno_format,"1")] <- 0.5
  gt[gt == paste0("1",geno_format,"1")] <- 1
  suppressWarnings(class(gt) <- "numeric")

  # Filter the genotypes based on popmap
  gt <- gt[,colnames(gt) %in% popmap[,1]]

  # For all populations, calculate the allele frequency of the ALT allele (it doesn't matter which we use)
  pops <- unique(unlist(popmap[,2]))

  pop_AF <- this.lapply(pops,function(pop){

    # Get only these pops...
    gt_tmp <- gt[,which(colnames(gt) %in% popmap[popmap[,2]==pop,1])]

    # Get the AFs
    allele_counts <- rowSums(gt_tmp,na.rm = T)

    # Get counts, and factor in missingness...
    total_counts <- ncol(gt_tmp) - apply(gt_tmp, 1, function(x) sum(is.na(x)))

    return(allele_counts/total_counts)
  })

  # Make a new matrix
  AF_mat <- matrix(ncol=length(pops),nrow=nrow(gt))
  for(i in 1:ncol(AF_mat)){
    AF_mat[,i] <- pop_AF[[i]]
  }
  colnames(AF_mat) <- pops

  # Get REF if needed
  if(allele=="REF"){
    AF_mat <- 1 - AF_mat
  }
  return(AF_mat)
}
