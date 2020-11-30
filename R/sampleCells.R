#' Stratified random sampling
#'
#' This is the function implementing the stratified random sampling approach
#'
#' @param strata_layer RasterLayer. A layer of candidates cells that can be sampled. The value of each cell should be the strata ID
#' @param matrix_strata A matrix obtained from the getPCAstrata function
#' @param existing_sample Optional. A data.frame with a column strata, a column x and a column y giving the XY coordinates of sampled cells centroids or existing plots.
#' @param toSample A data.frame with the first column refering to strat ID and second column the number of cells to sample in each strata
#' @param mindist Minimum distance between sampled cells
#' @param message Logical. Should messages be printed to console while the function is running
#' @param wrow Number of row in the focal window (default is 3)
#' @param wcol Number of columns in the focal window (default is 3)
#' @param checkStrata Logical. If TRUE, will calculate unique values of strata layer and will check if they match the ones provided in toSample. Default is FALSE
#' @export
#'
#' @importFrom raster %in%
#'

sampleCells <- function(strata_layer,
                        matrix_strata,
                        existing_sample,
                        toSample,
                        mindist,
                        message = TRUE,
                        wrow = 3,
                        wcol = 3,
                        checkStrata = FALSE) {

  if(!is(strata_layer, "RasterLayer")) {
    stop("strata_layer must be a RasterLayer object (opened with the raster function)")
  }

  if(missing(existing_sample)) {
    addSamples <- data.frame()
    extraCols <- character(0)
  }else{
    if(!is(existing_sample, "data.frame")) {
      stop("existing_samples must be a data.frame")
    }
    if(any(! c("strata", "x", "y") %in% colnames(existing_sample))) {
      stop("existing_samples must have columns named strata, x and y")
    }
    addSamples <- existing_sample
    extraCols <- colnames(existing_sample)[!colnames(existing_sample) %in% c("x", "y", "strata")]

    # Transform strata to numeric if factor
    if(is(addSamples$strata, "factor")) {
      addSamples$strata <- as.numeric(as.character(addSamples$strata))
    }

  #   # Assign NA to cells already sampled
  #   if (NROW(addSamples) > 0) {
  #     cells_toRemove <- raster::cellFromXY(strata_layer, matrix(c(addSamples$x, addSamples$y), nc = 2))
  #     strata_layer[cells_toRemove] <- NA
  #   }
   }

  if(!is(toSample, "data.frame")) {
    stop("toSample must be a data.frame")
  }else if (NCOL(toSample) < 2){
    stop("toSample must have at least two columns: strata and number of cells to sample")
  }

  # Check that all strata to sample correspond to a valid PCA_layer strata
  unique_strata_toSample <- unique(toSample[,1])
  if (checkStrata) {
    unique_strata_layer <- raster::unique(strata_layer)
    if(!all(unique_strata_toSample %in% unique_strata_layer)) {
      stop("Not all strata to sample correspond to a strata in PCA_layer")
    }
  }

  # Main loop through each strata

  for(i in 1:NROW(toSample)) {
    s <- toSample[i, 1] #strata
    need <- toSample[i,2] #number of cells to sample

    if (message) message(sprintf("Strata %s (%d/%d): %d cells to sample",s,i,length(unique_strata_toSample), need))

    if(need > 0) {
      # Mask all cells that are not strata
      group_s <- raster::mask(strata_layer, mask = strata_layer, maskvalue = s, inverse = TRUE)
      names(group_s) <- "strata"

      #Assign NA to cells already sampled

      # RULE 1: select only cells surrounded by cells with same strata
      w <- matrix(1/(wrow*wcol), nr = wrow, nc = wcol) #Focal window

      suppressWarnings(group_s_cluster <- raster::focal(group_s, w = w, na.rm = FALSE, pad = FALSE, NAonly = FALSE))
      names(group_s_cluster) <- "strata"

      values_strata_cluster <- is.na(getValues(group_s_cluster))
      group_s_cluster_cells <- 1:ncell(group_s_cluster)
      validCells_s_cluster <- group_s_cluster_cells[!values_strata_cluster]

      #Initiate number of sampled cells
      add_strata <- dplyr::filter(addSamples, strata == s)
      if(NROW(add_strata) > 0) {
        add_strata$type <- "Existing"
        if (! "cluster" %in% colnames(add_strata)){
          add_strata$cluster <- NA
        }
      }
      nCount <- 0 #Number of sampled cells

      if (message) message("Selecting in clusters ...")

      # While loop for RULE 1
      while(length(validCells_s_cluster) > 0 & nCount < need) {
        # Temporary sampled cells
        smp <- sample(1:length(validCells_s_cluster), size = 1)
        smp_cell <- validCells_s_cluster[smp]
        add_temp <- data.frame(cell = smp_cell,
                               x = raster::xFromCell(group_s_cluster, smp_cell),
                               y = raster::yFromCell(group_s_cluster, smp_cell),
                               strata = group_s_cluster[smp_cell])

        # Remove sampled cell from validCells so that it cannot be sampled again later
        validCells_s_cluster <- validCells_s_cluster[-smp]

        add_temp$type <- "New"
        add_temp$cluster <- "cluster"
        add_temp[,extraCols] <- NA
        # If add_strata is empty, sampled cell accepted
        if (NROW(add_strata) == 0) {
          add_strata <- add_temp[,c("x","y","strata","cluster","type",extraCols)]
          nCount = nCount +1
        }else{ # Else, check distance with all other sampled cells in strata
          dist <- spatstat::crossdist(add_temp$x,add_temp$y,add_strata$x,add_strata$y)
          if(all(as.numeric(dist)>mindist)){ # Is all less than mindist, accept sampled cell otherwise reject
            add_strata <- rbind(add_strata, add_temp[,c("x","y","strata","cluster","type",extraCols)])
            nCount = nCount + 1
          }
        }
      }

      # OUT OF RULE 1

      if (nCount < need) {
        #We still don't have enough samples, select within extended strata clusters (RULE 2)

        # Get strata that are located in the "neighbourhood" of current strata
        select_strata <- getNeiMatrix(matrix_strata, s)

        if (message) message("Selecting in extended clusters ...")


        # Mask cells not in extended strata
        mask_lyr <- strata_layer %in% select_strata
        group_s_ext <- raster::mask(strata_layer, mask_lyr, maskvalue = 0)

        # Select only cells in extended cluster
        suppressWarnings(group_s_ext_cluster <- raster::focal(group_s_ext, w = w, na.rm = FALSE, pad = FALSE))
        suppressWarnings(group_s_ext_cluster <- raster::mask(group_s_ext_cluster,group_s))
        names(group_s_ext_cluster) <- "strata"

        values_strata_ext <- is.na(getValues(group_s_ext_cluster))
        group_s_ext_cluster_cells <- 1:ncell(group_s_ext_cluster)
        validCells_s_ext_cluster <- group_s_ext_cluster_cells[!values_strata_ext]

        # While loop for RULE 2

        while(length(validCells_s_ext_cluster) > 0 & nCount < need){

          smp <- sample(1:length(validCells_s_ext_cluster), size = 1)
          smp_cell <- validCells_s_ext_cluster[smp]
          add_temp <- data.frame(cell = smp_cell,
                                 x = raster::xFromCell(group_s_ext_cluster, smp_cell),
                                 y = raster::yFromCell(group_s_ext_cluster, smp_cell),
                                 strata = group_s_ext_cluster[smp_cell])
          validCells_s_ext_cluster <- validCells_s_ext_cluster[-smp]

          add_temp$cluster <- "cluster extended"
          add_temp$type <- "New"
          add_temp[,extraCols] <- NA
          if (NROW(add_strata) == 0) {
            add_strata <- add_temp[,c("x","y","strata","cluster","type",extraCols)]
            nCount = nCount +1
          }else{
            dist <- spatstat::crossdist(add_temp$x,add_temp$y,add_strata$x,add_strata$y)
            if(all(as.numeric(dist)>mindist)){
              add_strata <- rbind(add_strata, add_temp[,c("x","y","strata","cluster","type",extraCols)])
              nCount = nCount +1
            }
          }
        }
      }

      if (nCount < need) {

        val_strata <- is.na(raster::getValues(group_s))
        group_s_cells <- 1:ncell(group_s)
        validCells_s <- group_s_cells[!val_strata]
        #We still don't have enough samples, select isolated

        if (message) message("Selecting isolated cells ...")

        while(length(validCells_s) > 0 & nCount < need){

          smp <- sample(1:length(validCells_s), size = 1)
          smp_cell <- validCells_s[smp]
          add_temp <- data.frame(cell = smp_cell,
                                 x = raster::xFromCell(group_s, smp_cell),
                                 y = raster::yFromCell(group_s, smp_cell),
                                 strata = group_s[smp_cell])
          validCells_s <- validCells_s[-smp]

          add_temp$cluster <- "isolated"
          add_temp$type <- "New"
          add_temp[,extraCols] <- NA
          if (NROW(add_strata) == 0) {
            add_strata <- add_temp[,c("x","y","strata","cluster","type",extraCols)]
            nCount = nCount +1
          }else{
            dist <- spatstat::crossdist(add_temp$x,add_temp$y,add_strata$x,add_strata$y)
            if(all(as.numeric(dist)>mindist)){
              add_strata <- rbind(add_strata, add_temp[,c("x","y","strata","cluster","type",extraCols)])
              nCount = nCount +1
            }
          }
        }
      }

      # SAMPLING FINISHED

      if (nCount < need){
        warning(sprintf("Strata %s: couldn't select required number of samples: %i instead of %i \n",s,nCount, need))
      }


      # If there are sampled cells in strata then add them to the ouput
      if (nrow(add_strata) > 0) {
        #Make sure that correct strata number is assigned and only correct columns are selected
        add_strata$strata <- s
        add_strata <- add_strata[,c("x","y","strata","cluster","type",extraCols)]


        # Create out object if first iteration of loop
        # Else just rbind output with what has been processed in the loop
        if (i == 1) {
          out <- add_strata
        }else{
          out <- rbind(out,add_strata)
        }
      }
    }
  }
  out
}
