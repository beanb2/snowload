# Test script to ensure that the new mapping functions are equivalent
# in results to the old ones.

test_that("prism predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- prism(variables = c("yr50", "ELEVATION"),
                  ut2017, utdem2, wbasin = c("HUC", 3),
                  bound = FALSE,
                  transform = "log")
  expect_equal(as.vector(output$prism), test_df$prism)
})

test_that("snlw predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- snlwf(utdem2)
  expect_equal(as.vector(output$snlw), test_df$snlw)
})

test_that("IDW predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- idw_snow(variables = c("yr50", "ELEVATION"),
                     ut2017, utdem2, bound_output = FALSE)
  expect_equal(as.vector(output$idw_snow), test_df$idw)
})

test_that("SKLM predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- rkriging(formula = log(yr50)~ELEVATION,
                     ut2017, utdem2, sklm = TRUE,
                     gstat::vgm(psill = .21, model = "Sph",
                                range = 200, nugget = .06))
  expect_equal(exp(as.vector(output$rkriging)), test_df$sklm)
})

test_that("UK predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- rkriging(formula = log(yr50)~ELEVATION,
                     ut2017, utdem2, sklm = FALSE,
                     gstat::vgm(psill = .21, model = "Sph",
                                range = 200, nugget = .06))
  expect_equal(exp(as.vector(output$rkriging)), test_df$uk)
})

test_that("tri_snow predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- tri_snow(c("yr50", "ELEVATION"),
                     ut2017, utdem2, density = c(100, 100),
                     NGSL = FALSE)
  expect_equal(as.vector(output$tri_snow), test_df$tri)
})

test_that("lm_snow predictions match parent function", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)

  output <- lm_snow(log(yr50) ~ ELEVATION, ut2017, utdem2)
  expect_equal(exp(as.vector(output$lm_snow)), test_df$lm)
})







