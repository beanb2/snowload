load("data-raw/newUtah2.R")
utsub <- newUtah.final[newUtah.final$yr50 != newUtah.final$yr50.2, ]

mean((newUtah.final$yr50 - newUtah.final$yr50.2) / newUtah.final$yr50.2)

# I incorrectly listed this stat as the 50 year load.
mean((utsub$maxobs - utsub$maxobs2) / utsub$maxobs2)
mean((utsub$yr50 - utsub$yr50.2) / utsub$yr50.2)

utsub[utsub$yr50 < 21, ] <- 21
utsub[utsub$yr50.2 < 21, ] <- 21

mean((utsub$yr50 - utsub$yr50.2) / utsub$yr50.2)
mean((utsub$maxobs - utsub$maxobs2) / utsub$maxobs2)

utsub[utsub$yr50 < 25, ] <- 25
utsub[utsub$yr50.2 < 25, ] <- 25

mean((utsub$yr50 - utsub$yr50.2) / utsub$yr50.2)
mean((utsub$maxobs - utsub$maxobs2) / utsub$maxobs2)
