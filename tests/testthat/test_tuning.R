test_that("prism tuning works as expected", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)
  set.seed(123)

  output1 <- prism(formula = log(yr50 + 1) ~ ELEVATION,
                   ut2017, utdem2,
                   distImp = 0.8,
                   minRad = c(10, 20),
                   wdistance = 2,
                   welevRange = list(lwr = c(200, 500), upr = c(1500, 2500)),
                   welevation = c("ELEVATION", 1, 2),
                   wbasin = c("HUC", 3, 4),
                   bound = FALSE,
                   tune_score = "MAE",
                   tune_folds = 10)
  output2 <- prism(formula = log(yr50 + 1) ~ ELEVATION,
                   ut2017, utdem2,
                   distImp = 0.8,
                   minRad = 10,
                   wdistance = 2,
                   welevRange = list(lwr = 200, upr = 1500),
                   welevation = c("ELEVATION", 2),
                   wbasin = c("HUC", 4),
                   bound = FALSE)

  expect_equal(output1, output2)
})

test_that("prism tuning returns lowest score", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)
  set.seed(123)

  output <- prism_tune(formula = log(yr50 + 1) ~ ELEVATION,
                       ut2017,
                       distImp = 0.8,
                       minRad = 10,
                       wdistance = 2,
                       welevRange = list(lwr = c(200, 500), upr = c(1500, 2500)),
                       welevation = c("ELEVATION", 1),
                       wbasin = c("HUC", 3, 4),
                       bound = FALSE,
                       tune_score = "MAE",
                       tune_folds = 10)

  expect_equal(unname(output[1, 9]), min(unname(output[, 9])))
})

test_that("rkriging auto variogram fitting working as expected", {
  data(ut2017)
  sp::coordinates(ut2017) <- c("LONGITUDE", "LATITUDE")
  sp::proj4string(ut2017) <- sp::proj4string(utdem2)
  set.seed(123)

  g <- gstat::gstat(NULL, "vario", log(yr50 + 1) ~ ELEVATION, ut2017)
  g_vario <- gstat::variogram(g)
  g_varioFit <- gstat::fit.variogram(g_vario, model = gstat::vgm("Sph"))

  output <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                     model = g_varioFit)
  output2 <- rkriging(log(yr50 + 1) ~ ELEVATION, ut2017, utdem2,
                      model = gstat::vgm("Sph"))

  expect_equal(output, output2)
})


