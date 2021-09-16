# Tests for eigen_pvals
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
emp_pvals <- eigen_pvals(eigen_res,null_vectors)

test_that("empirical pvals are expected", {
  expect_equal(nrow(emp_pvals), length(test_vectors))
  expect_equal(ncol(emp_pvals), length(vectors))
  expect_false(any(as.numeric(emp_pvals) > 1))
  expect_false(any(as.numeric(emp_pvals) < 0))
})
