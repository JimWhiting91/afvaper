# Tests for vcfR2AF
vcf_in <- vcfR::read.vcfR(system.file("full_parallel.vcf.gz",package="afvaper"))
vcf_in2 <- vcf_in[1:10001,]
popmap_in <- read.table(system.file("full_parallel.popmap",package="afvaper"))

# Fetch our AF matrix
AF_mat <- vcfR2AF(vcf_in2,popmap_in)

test_that("AFs look right...", {
  expect_equal(class(AF_mat), c("matrix","array"))
  expect_false(any(as.numeric(AF_mat) < 0))
  expect_false(any(as.numeric(AF_mat) > 1.0001))
  expect_equal(colnames(AF_mat),paste0("pop",1:5))
  expect_equal(nrow(AF_mat),nrow(vcf_in2@gt))
})
