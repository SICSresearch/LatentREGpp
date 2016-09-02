context("Input data")

test.invalidmatrix = function ( ) {
  print("Invalid matrixes")
  mt = matrix(NA, nrow = 100, ncol = 10)
  expect_error(latentreg(data = mt, dim = 1))
  mt = matrix(0, nrow = 100, ncol = 10)
  mt[1, 1] = 1
  mt[1, 2] = 2
  expect_error(latentreg(data = mt, dim = 1))
  mt = matrix(-1, nrow = 100, ncol = 10)
  expect_error(latentreg(data = mt, dim = 1))
}

test.invalidclusters = function ( ) {
  print("Invalid clusters")
  dir = system.file(package="LatentREGpp")
  folder = "/dataset/3D/dicho/"
  file = "1000x55-1.csv"
  size.cluster = c(20, 20, 15, 10)
  data_dir = paste(c(dir, folder, file), collapse = "")
  data = read.table(file = data_dir, sep = ";")
  expect_error(latentreg(data = data, dim = 3, clusters = size.cluster))
}

test_that(desc = "latentreg: Invalid matrix", code = test.invalidmatrix())
test_that(desc = "latentreg: Invalid clusters", code = test.invalidclusters())
