
test.latentreg.unidicho = function () {
  dir = system.file(package="LatentREGpp")
  file = "1000x50-1.csv"
  folder = "/dataset/1D/dicho/"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 1, verbose = FALSE)
  expect_identical(est$iterations, 38)
}

test_that(desc = "Unidicho test", code = test.latentreg.unidicho())