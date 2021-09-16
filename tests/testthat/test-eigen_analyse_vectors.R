# Tests for eigen_analyse_vectors
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

eigen_res <- lapply(test_vectors,eigen_analyse_vectors)

test_that("output format is correct", {
  expect_equal(length(eigen_res),length(test_vectors))
  expect_false(!(all(sapply(eigen_res,length)==3)))
  expect_false(!(all(names(eigen_res[[1]]) %in% c("eigenvals","eigenvecs","A_matrix"))))
  expect_false(!(all(sapply(eigen_res,function(x) return(nrow(x[[3]])))==200)))
  expect_false(!(all(sapply(eigen_res,function(x) return(length(x[[1]])))==length(vectors))))
  expect_false(!(all(sapply(eigen_res,function(x) return(nrow(x[[2]])))==length(vectors))))
  expect_false(!(all(sapply(eigen_res,function(x) return(ncol(x[[2]])))==length(vectors))))
  expect_false(!(all(round(sapply(eigen_res,function(x) return(sum(x[[1]]))),4)==length(vectors))))
})
