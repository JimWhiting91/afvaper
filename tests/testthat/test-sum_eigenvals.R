# Tests for sum_eigenvals
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
summed_eigenvals <- sum_eigenvals(eigen_res)

test_that("output looks right...", {
  expect_equal(class(summed_eigenvals),c("matrix","array"))
  expect_equal(nrow(summed_eigenvals),length(test_vectors))
  expect_equal(ncol(summed_eigenvals),length(vectors))
  expect_false(any(as.numeric(summed_eigenvals) > length(vectors)+0.0001))
})