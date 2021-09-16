# Tests for summarise_window_parallelism
vcf_in <- vcfR::read.vcfR(system.file("full_parallel.vcf.gz",package="afvaper"))
popmap_in <- read.table(system.file("full_parallel.popmap",package="afvaper"))
# Truncate for speed
vcf_in2 <- vcf_in[1:10001,]
vectors <- lapply(2:5,function(x) return(c("pop1",paste0("pop",x))))
names(vectors) <- paste0("pop",2:5)
test_vectors <- calc_AF_vectors(vcf = vcf_in2,
                                popmap = popmap_in,
                                vectors = vectors,
                                n_cores = 1)
null_vectors <- calc_AF_vectors(vcf = vcf_in2,
                                popmap = popmap_in,
                                vectors = vectors,
                                n_cores = 1,
                                null_perms=100)

eigen_res <- lapply(test_vectors,eigen_analyse_vectors)

null_cutoffs <- find_null_cutoff(null_vectors,cutoffs = c(0.95,0.99,0.999))

# Find our signif windows
signif_windows <- signif_eigen_windows(eigen_res,cutoffs = null_cutoffs[,1])

# Summarise
eig1_summary <- summarise_window_parallelism(signif_windows[[1]],
                                             eigen_res = eigen_res,
                                             eigenvector = 1)
eig2_summary <- summarise_window_parallelism(signif_windows[[2]],
                                             eigen_res = eigen_res,
                                             eigenvector = 2)
eig4_summary <- summarise_window_parallelism(signif_windows[[4]],
                                             eigen_res = eigen_res,
                                             eigenvector = 4)

test_that("outputs are as expected...", {
  expect_equal(class(eig1_summary),c("data.frame"))
  expect_equal(unique(eig2_summary$eigenvector),paste0("Eig",1:2))
  expect_equal(unique(eig4_summary$eigenvector),paste0("Eig",1:4))
  expect_equal(nrow(eig1_summary),length(signif_windows[[1]]))
  expect_equal(colnames(eig1_summary),c("window_id","eigenvector","eigenvalue","parallel_lineages","parallel_pops","antiparallel_pops"))
  expect_equal(colnames(eig2_summary),c("window_id","eigenvector","eigenvalue","eigenvalue_sum","parallel_lineages","parallel_pops","antiparallel_pops"))
})
