test_that("rkriging weight option works as expected", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                     model = gstat::vgm("Sph"), weights = FALSE)
  output2 <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                      model = gstat::vgm("Sph"), weights = TRUE)

  expect_equal(output, output2[[1]])
})

test_that("rkriging weight option works as expected (sklm)", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                     model = gstat::vgm("Sph"), sklm = TRUE,
                     weights = FALSE)
  output2 <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                      model = gstat::vgm("Sph"), sklm = TRUE,
                      weights = TRUE)

  expect_equal(output, output2[[1]])
})




