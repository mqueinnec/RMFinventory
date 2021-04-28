#' Define or retrieve strata from a layer with PCA values
#'
#' @param PCA_layer update
#' @param method update
#' @param nbreaks update
#' @param breaks update
#' @param summary update
#'
#' @importFrom dplyr n
#'
#' @export

getPCAstrata <- function(PCA_layer,
                         method = "equal_interval",
                         nbreaks,
                         breaks,
                         summary = FALSE) {

  returnRaster = FALSE #Logical used to indicate if a Raster object is returned

  if (!missing(breaks)) {
    message("breaks have been provided and won't be recalculated")
    if(!is(breaks,"list")) {
      stop("breaks must be a list")
    }
    nbreaks <- as.vector(sapply(breaks, length))

    if (any(c(is(PCA_layer,"RasterLayer"), is(PCA_layer,"RasterBrick"), is(PCA_layer,"RasterStack")))) {
      PCA_layer_r <- PCA_layer
      PCA_layer <- raster::as.data.frame(PCA_layer_r, na.rm = TRUE, xy = TRUE)
      cells <- raster::cellFromXY(raster::subset(PCA_layer_r,1), as.matrix(PCA_layer[,c("x","y")]))
      PCA_layer <- PCA_layer[,!colnames(PCA_layer) %in% c("x","y")]
      returnRaster <- TRUE
    }else if (!is(PCA_layer,"data.frame")) {
      stop("PCA_layer must be a data.frame or Raster* object")
    }

    nfeatures <- length(breaks)
    features_names <- names(breaks)

  }else{

    if (any(c(is(PCA_layer,"RasterLayer"), is(PCA_layer,"RasterBrick"), is(PCA_layer,"RasterStack")))) {
      PCA_layer_r <- PCA_layer
      PCA_layer <- raster::as.data.frame(PCA_layer_r, na.rm = TRUE, xy = TRUE)
      cells <- raster::cellFromXY(raster::subset(PCA_layer_r,1), as.matrix(PCA_layer[,c("x","y")]))
      PCA_layer <- PCA_layer[,!colnames(PCA_layer) %in% c("x","y")]
      returnRaster <- TRUE
    }else if (!is(PCA_layer,"data.frame")) {
      stop("PCA_layer must be a data.frame or Raster* object")
    }

    nfeatures <- NCOL(PCA_layer)
    features_names <- colnames(PCA_layer)

    if (method == "equal_interval") {

      if(missing(nbreaks)) {
        stop("Must provide nbreaks to get strata if method = equal interval")
      }else if (length(nbreaks) != nfeatures) {
        stop("nbreaks must have the same length as the number of PCA_raster layers")
      }

      breaks <- list()
      for (n in 1:nfeatures) {
        breaks[[features_names[n]]] <- seq(min(PCA_layer[,features_names[n]]),
                                           max(PCA_layer[,features_names[n]]),
                                           length.out = nbreaks[n])


      }
    }

  }

  # Assingn strata to df
  all_strata <- do.call(tidyr::crossing, lapply(breaks, function(x) seq(from = 1, length.out = length(x) - 1)))
  # Add strata ID
  all_strata$strata <- 1:NROW(all_strata)
  # all_strata <- do.call(paste0, as.list(all_strata))

  strata <- lapply(seq_along(breaks), function(x) {findInterval(PCA_layer[[names(breaks)[x]]], breaks[[x]], rightmost.closed = TRUE)})
  names(strata) <- names(breaks)
  strata <- strata[features_names] # Make sure that list is ordered

  PCA_strata <- strata

  strata_df <- as.data.frame(strata)
  strata_df$rowID <- 1:NROW(strata_df)

  strata_ID_df <- merge(strata_df, all_strata, sort = FALSE)

  strata_ID_df$strata <- factor(strata_ID_df$strata, levels = all_strata$strata, ordered = T)

  # Make sure that strata_df is ordered by rowID
  strata_ID_df <- dplyr::arrange(strata_ID_df, rowID)

  PCA_layer <- strata_ID_df

  if (summary) {
    PCA_layer_summary <- strata_ID_df %>%
      dplyr::group_by(strata) %>%
      dplyr::summarize(count = n(), frac = n()/NROW(strata_ID_df)) %>%
      tidyr::complete(strata, fill = list(count = 0, frac = 0))
  }

  if (returnRaster) {
    strata_ID_df$cells <- cells
    seq_cells <- seq(from = 1, to = raster::ncell(raster::subset(PCA_layer_r,1)))
    # Assign back NA values to cells
    if (NROW(strata_ID_df) != length(seq_cells)) {
      df_temp <- rbind(strata_ID_df[,c("strata", "cells")], data.frame(strata = NA, cells = seq_cells[!(seq_cells %in% strata_ID_df$cells)]))
    }
    df_temp <- dplyr::arrange(df_temp, cells)
    PCA_layer <- raster::setValues(raster::subset(PCA_layer_r,1), as.numeric(as.character(df_temp$strata)))
    names(PCA_layer) <- "strata"
  }

  if(summary) {
    out <- list(strata_layer = PCA_layer$strata,
                PCA_strata = PCA_strata,
                breaks = breaks,
                summary = PCA_layer_summary)
  }else{
    out <- list(strata_layer = PCA_layer$strata,
                PCA_strata = PCA_strata,
                breaks = breaks)
  }
  out$matrix <- makeStrataMatrix(nbreaks, all_strata)
  out
}
