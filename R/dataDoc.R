#' The National Snow Load Dataset
#'
#' 50 year ground snow load estimates at 11,444 measurement locations
#' in North America.
#'
#' @format A data frame with 11,444 rows and 14 variables:
#'
#'  \describe{
#'  \item{STATION}{The numeric station identifier. Actual station names are
#'  not provided as each location may represent more than one station.}
#'  \item{STATION_NAME}{The character name of the first station
#'  (in alphabetical order) in the cluster.}
#'  \item{STATE}{Two character state abbreviation of the
#'  measurement location. Includes Canadian provinces.}
#'  \item{LONGITUDE}{Longitude coordinate position.}
#'  \item{LATITUDE}{Latitude coordinate position.}
#'  \item{ELEVATION}{Elevation of the measurement location (meters).}
#'  \item{REGIONIII}{U.S. Environmental Protection Agency level III
#'  numeric ecoregion.}
#'  \item{climate}{Climate class based on Sturm et al. (1995).}
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  \item{numST}{The number of stations comprising the measurement location.}
#'  \item{numMax}{The number of annual maximums that
#'  were retained from the period of record (i.e. the sample size).}
#'  \item{medMax}{The median of the annual maximum snow loads from the
#'  historical record.}
#'  \item{mMax}{The maximum snow load observation from the
#'  period of record (kpa).}
#'  \item{snowload}{The estimated design (i.e. 50 year) ground snow load (kpa).}
#'  }
#'
#' @details This dataset was created by scraping data from the global
#' historical climatological network. Values were subject to automatic
#' outlier detection and removal. 50 year estimates were obtained by
#' fitting 5 distributions (log-normal, gamma, extreme value distributions
#' types I, II, and III) and selecting the median 98th percentile
#' from the 5 distributions. Climate types were determined using
#' eco-region maps from the Environmental Protection Agency (EPA)
#' and climate classes developed by Sturm et al. (1995).
#'
#' Note that each measurement location may represent more than
#' one station. When two stations had overlapping period of
#' records, only the maximum snow load estimate was retained.
#' Preference was also always given to direct measurements
#' of snow load.
#'
#' @references
#' \insertRef{eco2018}{snowload}
#' \insertRef{Sturm1995}{snowload}
"snowloads"


#' The 2017 Utah Snow Load Study Dataset
#'
#' 50 year ground snow load estimates at 415 measurement locations in
#' locations in and around the state of Utah. Used as part of the
#' 2018 Utah snow load study.
#'
#' @format A data frame with 415 rows and 10 variables:
#'
#'  \describe{
#'  \item{STATION}{The National Weather Service COOP, or Natural Resource
#'  conservation service SNOTEL station identifier.}
#'  \item{STATION_NAME}{The character name of the station.}
#'  \item{STATE}{Two character state abbreviation of the
#'  measurement location.}
#'  \item{LONGITUDE}{Longitude coordinate position.}
#'  \item{LATITUDE}{Latitude coordinate position.}
#'  \item{ELEVATION}{Elevation of the measurement location (meters).}
#'  \item{YRS}{The number of years of recorded load maximums. This is
#'  effectively the sample size of the 50 year estimate.}
#'  \item{maxobs}{The maximum snow load observation from the
#'  period of record (kpa).}
#'  \item{yr50}{The estimated design (i.e. 50 year) ground snow load (kpa).}
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  }
#'
#' @details This dataset contains 279 (192 COOP, 87 SNOTEL) Utah stations
#' with an additional 136 stations (103 COOP, 33 SNOTEL), all located within
#' 100 km of the Utah border. 50 year estimates were obtained by fitting
#' a log-normal distribution to the annual maximum snow loads (by water
#' year, not calendar year) via maximum likelihood estimation and
#' extracting the 98th percentile.
#'
#' @references
#' \insertRef{Bean2018-report}{snowload}
"ut2017"

#' The 1992 Utah Snow Load Study Dataset
#'
#' 50 year snow load estimates at 413 measurement locations
#' (210 Snow Course, 203 COOP) in Utah. Used to determine Utah's
#' 1992 ground snow load requirements.
#'
#' @format A data frame with 413 rows and 10 variables:
#'
#'  \describe{
#'  \item{STATION}{The National Weather Service COOP, or Natural Resource
#'  conservation service Snow Course station identifier.}
#'  \item{STATION_NAME}{The character name of the station.}
#'  \item{LONGITUDE}{Longitude coordinate position.}
#'  \item{LATITUDE}{Latitude coordinate position.}
#'  \item{ELEVATION}{Elevation of the measurement location (meters).}
#'  \item{approx}{Indicates whether or not the station location was
#'  approximated using Google Earth (1-yes, 0-no).}
#'  \item{maxobs}{The maximum snow load observation from the
#'  period of record (kpa).}
#'  \item{yr50}{The estimated design (i.e. 50 year) ground snow load (kpa).}
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  }
#'
#' @details The data include
#' 50 year estimates were determined using a log-Pearson III
#' applied to annual maximum snow loads. Where required, snow loads were
#' estimated from snow depths using the rocky mountain conversion density.
#'
#' Station locations were not provided in the original data. Locations were
#' inferred from NWS station databases, personal correspondence with the
#' Utah Snow Survey Office in Salt Lake City, and Google Earth.
#'
#' @references
#' \insertRef{SEAU1992}{snowload}
"ut1992"


#' The 2015 Idaho Snow Load Study Dataset
#'
#' 50 year ground snow load estimates at 651 measurement locations
#' in and around the state of Idaho. Used as part of the
#' 2015 Idaho snow load study.
#'
#' @format A data frame with 651 rows and 9 variables:
#'
#'  \describe{
#'  \item{STATION}{The National Weather Service COOP, or Natural Resource
#'  conservation service Snow Course station identifier.}
#'  \item{STATION_NAME}{The character name of the station.}
#'  \item{STATE}{Two character state abbreviation of the
#'  measurement location.}
#'  \item{LONGITUDE}{Longitude coordinate position.}
#'  \item{LATITUDE}{Latitude coordinate position.}
#'  \item{ELEVATION}{Elevation of the measurement location (meters).}
#'  \item{YRS}{The number of years of recorded load maximums. This is
#'  effectively the sample size of the 50 year estimate.}
#'  \item{yr50}{The estimated design (i.e. 50 year) ground snow load (kpa).}
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  }
#'
#' @details These data consist of 394 (246 SC/SNOTEL, 148 COOP) Idaho stations
#' with an additional 257 (222 SC/SNOTEL, 35 COOP) located near the Idaho
#' border with the most recent measurements being taken in 2014. log-Pearson
#' type III distribution parameter estimates were determined using the sample
#' mean, standard deviation skew of annual of yearly maximums at each station
#' location (i.e. method of moments).
#'
#' @references
#' \insertRef{Hatailah2015}{snowload}
"id2015"

#' A 1km by 1km digital elevation model for Utah.
#'
#' These elevation values were obtained from the U.S.
#' Geological Survey National Map. They were aggregated from
#' their original 30m by 30m resolution to a 1km by 1km
#' resultion using \code{\link[raster]{aggregate}}.
#'
#'
#' @format A SpatialPixelsDataFrame with the following variables.
#'
#'  \describe{
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  \item{ELEVATION}{The elevation of the location (meters).}
#'  }
#'
"utdem"

#' Post office locations in Utah
#'
#' These 253 post office locations provide the ideal set
#' of prediction locations to test the snow load functions.
#' These locations are also the same locations comprising the
#' city tables provided in the 2018 Utah Snow Load Study.
#'
#' @format A data frame with 253 rows and 7 variables:
#'
#'  \describe{
#'  \item{CITY}{The character name of the post office city.}
#'  \item{COUNTY}{The character name of the county.}
#'  \item{STATE}{The two-character state abbrevation.}
#'  \item{LONGITUDE}{Longitude coordinate position.}
#'  \item{LATITUDE}{Latitude coordinate position.}
#'  \item{ELEVATION}{Elevation of the measurement location (meters).}
#'  \item{HUC}{The United States Geological Survey 8 digit
#'  hydrologic unit code.}
#'  }
#'
#' @references
#' \insertRef{Bean2018-report}{snowload}
"utpost"

