
#from https://github.com/gdemin/maditr/blob/master/tests/tinytest.R
if(requireNamespace("tinytest", quietly=TRUE)){
  library(tinytest)
  #options(covr = FALSE)
  data.table::setDTthreads(0)
  tinytest::test_package("fastdid", remove_side_effects=FALSE)
  data.table::setDTthreads(NULL)
}