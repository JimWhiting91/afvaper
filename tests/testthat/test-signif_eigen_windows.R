# Tests for signif_eigen_windows
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

test_that("signif windows output is as expected", {
  expect_equal(class(signif_windows), c("list"))
  expect_equal(class(signif_windows[[1]]), c("character"))
  expect_equal(length(signif_windows), length(vectors))
})
