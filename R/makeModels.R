#' ABA models and accuracy assessment
#'
#' Based on plot-level forest attributes and correpsonding ALS metrics, this function train predictive models based on a modeling approach (currently only random forest is implemented) and assess models accuracy using k-folds cross-validation
#'
#' Further details
#'
#'@param dat data.frame. Need to contain both forest attributes and ALS metrics
#'@param attNames Character. Name of attrbiutes to be modeled
#'@param preds Character. Name of ALS metrics to be included
#'@param k Numeric. Number of folds to be created for k-folds cross-validation. Default is 5
#'@param saveModel Logical. If TRUE, models will be saved in outdir in .rds format. Default if FALSE
#'@param saveFigure Logical. If TRUE, save scatterplots in outdir. Currently save a pdf file of 4 inches width and 4 inches height
#'@param titles List containing titles of scatterplots. Elements must be named according to attNames. If missing, scatterplots won't have any title
#'@param outdir Character. Path to exisiting directory where models and scatterplots will be saved if wanted
#'@return A list with one element per modeled forest attribute and one element containing accuracy measures (R2, RMSE and bias).
#'@export
#'@importFrom dplyr %>%

makeModels <- function(dat,
                       attNames,
                       preds,
                       k = 5,
                       titles,
                       saveModel,
                       saveFigure,
                       outdir) {

  if (missing(dat)) {
    stop("dat must be provided")
  }

  if (missing(attNames)) {
    stop("attNames needs to be provided")
  }else if (!all(attNames %in% colnames(dat))){
    stop("Not all attNames correspond to dat column names")
  }

  if (missing(preds)) {
    stop("preds needs to be provided")
  }else if (!all(preds %in% colnames(dat))){
    stop("Not all attNames correspond to dat column names")
  }

  if (missing(titles)) {
    titles <- NULL
  }else{
    if(!all(attNames %in% names(titles))) {
      stop("titles must be a list with elements names corresponding to attNames")
    }
  }

  if (missing(saveModel)) {
    saveModel <- FALSE
  }else if(!is.logical(saveModel)) {
    stop("saveModel must be a logical")
  }

  if (missing(saveFigure)) {
    saveFigure <- FALSE
  }else if(!is.logical(saveFigure)) {
    stop("saveFigure must be a logical")
  }

  if(missing(outdir)) {
    outdir = getwd()
  }else if(!is.character(outdir)){
    stop("outdir must be a character (path to directory)")
  }else if (!dir.exists(outdir)) {
    stop("outdir must be a valid path to directory")
  }

  # Prepare dat
  dat <- dat[complete.cases(dat),]
  if (NROW(dat) == 0) {stop("dat only contains non-complete cases (rows with at least one NA)")}

  # Variables to be returned
  out <- list() # List that will be returned
  all_preds_combined <- data.frame()

  # Loop through each attNames
  for (n in attNames) {
    # Data frame with reponse variable and ALS metrics
    df <- dat[,c(n, preds)]

    # Range of mtry
    nx <- length(preds)
    mtry_tune <- data.frame(mtry = seq(from = 1, to = nx, by = 2))

    # Set up cross validation and traning parameters
    folds <- caret::createFolds(df[,n], k = k, returnTrain = TRUE)

    trControl <- caret::trainControl(method = "cv",
                                     index = folds,
                                     number = k,
                                     savePredictions = "all")

    # Train RF model
    rf_folds <- caret::train(eval(parse(text = paste0(n," ~ ."))),
                             data = df,
                             method = "rf",
                             trControl = trControl,
                             importance = TRUE)

    out[[n]] <- rf_folds

    #Save model
    if(saveModel) {
      saveRDS(rf_folds, file = file.path(outdir,paste0("rf_5folds_",n,".rds")))
    }

    # Accuracy
    all_preds <- dplyr::filter(rf_folds$pred, mtry == as.numeric(rf_folds$bestTune))
    all_preds$variable <- n
    all_preds_combined <- rbind(all_preds_combined, all_preds)

    # Make scatterplots
    p <- RMFinventory::scatter(x = all_preds$obs, y = all_preds$pred, title = titles[[n]], label_text = c("bias","bias%","RMSE","RMSE%") )

    if (saveFigure) {
      pdf(file = file.path(outdir, paste0(n,"_scatter.pdf")),
          useDingbats = FALSE,
          height = 4,
          width = 4)
      print(p)
      dev.off()
    }


    out[[n]][["model"]] <- rf_folds
    out[[n]][["plot"]] <- p

  }
  accuracy_folds <- all_preds_combined %>%
    group_by(Resample, variable) %>%
    summarise(r2 = caret::R2(pred = pred, obs = obs, formula = "traditional"),
              rmse = caret::RMSE(obs = obs, pred = pred),
              rmse_per = rmse / mean(obs) * 100,
              bias = sum(pred - obs)/n(), # or mean(pred - obs)
              bias_per = bias / mean(obs) * 100)

  accuracy_summary <- accuracy_folds %>%
    group_by(variable) %>%
    summarise(r2_m = mean(r2),
              r2_sd = sd(r2),
              rmse_m = mean(rmse),
              rmse_sd = sd(rmse),
              rmse_per_m = mean(rmse_per),
              rmse_per_sd = sd(rmse_per),
              bias_m = mean(bias),
              bias_sd = sd(bias),
              bias_per_m = mean(bias_per),
              bias_per_sd = sd(bias_per))

  out[["accuracy"]][["all folds"]] <- accuracy_folds
  out[["accuracy"]][["summary"]] <- accuracy_summary

  return(out)
}
