context("Item parameters estimation")

test.itemfit.unidicho = function () {
  print("1D 1000x50 dicho")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/1D/dicho/"
  file = "1000x50-1.csv"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = itemfit(data = data, dim = 1, model = "1PL", EMepsilon = 0.001)
  expect_identical(est$iterations, 12)
}

test.itemfit.unipoly = function () {
  print("1D 1000x50 poly")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/1D/poly/"
  file = "1000x50-1.csv"
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = itemfit(data = data, dim = 1, model = "1PL")
  expect_identical(est$iterations, 35)
}

test.itemfit.multidicho = function () {
  print("3D 1000x55 dicho")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/dicho/"
  file = "1000x55-1.csv"
  size.cluster = c(20, 20, 15)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = itemfit(data = data, dim = 3, clusters = size.cluster, EMepsilon = 0.001, model = "1PL")
  expect_identical(est$iterations, 12)
}

test.itemfit.multipoly = function () {
  print("3D 1000x60 poly")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/poly/"
  file = "1000x60-1.csv"
  size.cluster = c(20, 20, 20)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  est = itemfit(data = data, dim = 3, clusters = size.cluster, EMepsilon = 0.01, model = "1PL")
  expect_identical(est$iterations, 5)
}

test_that(desc = "itemfit: Unidimensional dichotomous test", code = test.itemfit.unidicho())
test_that(desc = "itemfit: Unidimensional polytomous test", code = test.itemfit.unipoly())
test_that(desc = "itemfit: Multidimensional dichotomous test", code = test.itemfit.multidicho())
test_that(desc = "itemfit: Multidimensional polytomous test", code = test.itemfit.multipoly())
