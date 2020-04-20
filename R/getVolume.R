#' Calculate whole stem volume and merchantable volume
#'
#' @param species Tree species code. See details for the list of accepted species
#' @param H Tree height in m
#' @param DBH Diameter at breast height in cm
#' @param h_DBH Breast height in m. Default is 1.3
#' @param stumpHeight Stump height in m. Is used as a lower limit for the calculation of merchantable volume
#' @param min_d Minimum stem diamter in cm. Used as upper limit for the calculation of merchantable volume
#' @param source Character. Name of published allometric equations. See details for list of accepted equations
#'

getVolume <- function(species,
                      H,
                      DBH,
                      h_dbh = 1.3,
                      stumpHeight,
                      min_d,
                      source) {



  }
