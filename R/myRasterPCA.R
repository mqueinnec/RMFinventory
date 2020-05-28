#' Calculate PCA of a raster from covariance matrix
#'
#' @param x Raster object
#' @param covMat Matrix of covariance with rownames and colnames correspondign to raster layer names
#' @param xmeans Mean of each layer of x, as a named numeric vector
#' @param threads Number of parallel threads
#' @param ... Additional arguments passed to writeRaster
#'
#' @export
#'

myRasterPCA <- function(x,
                        covMat,
                        xmeans,
                        nComp = 2,
                        threads = 1,
                        ...) {

  #Check names of covMat and x
  if(!any(names(x) %in% rownames(covMat))) stop("covMat rownames and x layers names not matching")
  if(!any(names(x) %in% colnames(covMat))) stop("covMat colnames and x layers names not matching")

  ellip <- list(...)

  model  <- stats::princomp(covmat = covMat, cor=TRUE)
  model$center <- xmeans
  model$n.obs  <- ncell(x)

  S <- diag(covMat)
  model$scale <- sqrt(S * (model$n.obs-1)/model$n.obs)

  #Predict
  if(threads > 1) {
    raster::beginCluster(threads)
    if("filename" %in% names(ellip)){
      out <- raster::clusterR(x, fun = raster::predict, args = c(list(model = model,
                                                                      na.rm = TRUE,
                                                                      index = 1:nComp)))

      writeRaster(out, ...)
    }else{
      out <- raster::clusterR(x, fun = raster::predict, args = c(list(model = model,
                                                                      na.rm = TRUE,
                                                                      index = 1:nComp),
                                                                 ellip))
    }

    raster::endCluster()

  }else{
    out <- raster::predict(x, model = model, na.rm = TRUE, index = 1:nComp, ...)
  }

  names(out) <- paste0("PC",1:nComp)

  return(list(map = out,
              model = model))

}
