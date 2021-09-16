# Tests for find_null_cutoff
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

test_that("we can find null cutoffs and get expected output...", {
  expect_equal(class(null_cutoffs), c("matrix","array"))
  expect_equal(nrow(null_cutoffs),length(vectors))
  expect_equal(ncol(null_cutoffs),3)
  expect_false(any(as.numeric(null_cutoffs) > length(vectors)+0.0001))
  expect_equal(colnames(null_cutoffs),paste0(c(95,99,99.9),"%"))
  expect_equal(rownames(null_cutoffs),paste0("Eigenvector ",1:length(vectors)))
})
