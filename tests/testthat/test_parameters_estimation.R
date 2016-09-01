
test.latentreg.unidicho = function () {
  print("1D 1000x50 dicho")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/1D/dicho/"
  file = "1000x50-1.csv"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 1)
  expect_identical(est$iterations, 38)
}

test.latentreg.unipoly = function () {
  print("1D 1000x50 poly")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/1D/poly/"
  file = "1000x50-1.csv"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 1)
  expect_identical(est$iterations, 41)
}

test.latentreg.multidicho = function () {
  print("3D 1000x55 dicho")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/dicho/"
  file = "1000x55-1.csv"
  size.cluster = c(20, 20, 15)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 3, clusters = size.cluster)
  expect_identical(est$iterations, 46)
}

test.latentreg.multipoly = function () {
  print("3D 1000x60 poly")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/poly/"
  file = "1000x60-1.csv"
  size.cluster = c(20, 20, 20)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = latentreg(data = data, dim = 3, clusters = size.cluster)
  expect_identical(est$iterations, 71)
}

test_that(desc = "latentreg: Unidimensional dichotomous test", code = test.latentreg.unidicho())
test_that(desc = "latentreg: Unidimensional polytomous test", code = test.latentreg.unipoly())
test_that(desc = "latentreg: Multidimensional dichotomous test", code = test.latentreg.multidicho())
test_that(desc = "latentreg: Multidimensional polytomous test", code = test.latentreg.multipoly())
