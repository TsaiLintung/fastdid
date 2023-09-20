library(pkgKitten)
library(tinytest)
library(roxygen2)

pkgKitten::kitten("hihi")
tinytest::setup_tinytest("hihi")

test_all("hihi")

build_install_test("hihi")