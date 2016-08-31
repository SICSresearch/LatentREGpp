
test.latentreg.unidicho = function () {
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/1D/dicho/"
  file = "1000x50-1.csv"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 1)
  expect_identical(est$iterations, 38)
}

test.latentreg.multidicho = function () {
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/dicho/"
  file = "1000x55-1.csv"
  size.cluster = c(20, 20, 15)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 3, clusters = size.cluster)
  expect_identical(est$iterations, 46)
}

test_that(desc = "Unidimensional dichotomous test", code = test.latentreg.unidicho())
test_that(desc = "Multidimensional dichotomous test", code = test.latentreg.multidicho())
