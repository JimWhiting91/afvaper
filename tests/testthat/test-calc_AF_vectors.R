# Tests for calc_AF_vectors()
vcf_in <- vcfR::read.vcfR(system.file("full_parallel.vcf.gz",package="afvaper"))
popmap_in <- read.table(system.file("full_parallel.popmap",package="afvaper"))

test_that("calc_AF_vectors() calculates AF vectors from VCF", {
  vectors <- lapply(2:5,function(x) return(c("pop1",paste0("pop",x))))
  names(vectors) <- paste0("pop",2:5)
  test_vectors <- calc_AF_vectors(vcf = vcf_in,
                                  popmap = popmap_in,
                                  vectors = vectors,
                                  n_cores = 1)
  
  # Are vectors normalised?
  expect_equal(sum(test_vectors[[1]][1,]^2)^0.5, 1)
  # Do we have all vectors?
  expect_equal(nrow(test_vectors[[1]]),length(vectors))
  expect_equal(length(test_vectors),floor(nrow(vcf_in@fix)/200))
})

test_that("calc_AF_vectors() can calculate without normalising", {

  # Truncate for speed
  vcf_in2 <- vcf_in[1:601,]
  vectors <- lapply(2:5,function(x) return(c("pop1",paste0("pop",x))))
  names(vectors) <- paste0("pop",2:5)
  test_vectors <- calc_AF_vectors(vcf = vcf_in2,
                                  popmap = popmap_in,
                                  vectors = vectors,
                                  normalise = F,
                                  n_cores = 1)
  
  # Are vectors normalised?
  expect_false(isTRUE(all.equal(sum(test_vectors[[1]][1,]^2)^0.5, 1)))
  
  # Make sure that there are no values bigger or smaller than 1/-1
  all_vals <- as.numeric(unlist(test_vectors))
  expect_false(any(abs(all_vals)>1))
})

# Test null perms
test_that("calc_AF_vectors() can perform null perms as expected...", {
  
  # Truncate for speed
  vcf_in2 <- vcf_in[1:600,]
  popmap_in <- read.table(system.file("full_parallel.popmap",package="afvaper"))
  vectors <- lapply(2:5,function(x) return(c("pop1",paste0("pop",x))))
  names(vectors) <- paste0("pop",2:5)
  test_vectors <- calc_AF_vectors(vcf = vcf_in2,
                                  popmap = popmap_in,
                                  vectors = vectors,
                                  normalise = T,
                                  null_perms = 4,
                                  n_cores = 1)
  
  # Fetch start_pos
  pos <- names(test_vectors)
  pos <- sapply(strsplit(pos,":"),'[[',2)
  start_pos <- as.integer(sapply(strsplit(pos,"-"),'[[',1))
  
  # Are positions shuffled...
  expect_false(all(start_pos %in% c(1,201,401)))
})
