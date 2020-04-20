#' Get random forest importance values and rank variables per order of importance
#'
#'@param rfList List of random forest models trained with the function \code{\link[RMFinventory]{makeModels}}
#'@param measure Character. Random forest importance measure. Can be either IncMSE (default) or IncNodePurity.
#'@param nrank Numeric. Maximum importance rank to be included. Default is 5
#'@param saveFigure Logical. If TRUE, save scatterplots in outdir. Currently save a pdf file of 7 inches width and 5 inches height. Default is FALSE
#'@param outdir Character. Path to exisiting directory where models and scatterplots will be saved if wanted
#'@return A list containing two elements:
#'\itemize{
#'\item{"rank"}{ a data.frame showing the number of time a variable is part of the nrank most important varaibles}
#'\item{"plot}{a barplot (ggplot2 object)}
#'}
#'@export



varImpRank <- function(rfList,
                       measure = "IncMSE",
                       nrank = 5,
                       saveFigure = FALSE,
                       outdir) {

  if (length(measure) != 1) {
    stop("measure must be either IncMSE or IncNodePurity")
  }

  if (!measure %in% c("IncMSE","IncNodePurity")) {
    stop("measure must be either IncMSE or IncNodePurity")
  }

  if (length(nrank) > 1) {
    stop("nrank must be a single value")
  }

  rf_imp <- lapply(rfList, function(x) as.data.frame(randomForest::importance(x$model$finalModel)))

  rf_imp <- lapply(rf_imp, function(x) {
    colnames(x) <- c("IncMSE","IncNodePurity")
    x$metric <- row.names(x)
    x <- dplyr::arrange(x,eval(parse(text = paste0("desc(",measure,")"))))
    x$rank = 1:NROW(x)
    x
  })

  nrank <- 1:round(nrank,0)

  rf_imp_summary <- dplyr::group_by(bind_rows(rf_imp),metric)
  rf_imp_summary <- eval(parse(text=paste0("dplyr::summarize(rf_imp_summary,",paste0("n",nrank," = sum(rank == ",nrank,")", collapse = ","),", tot_n = sum(",paste0("n",nrank, collapse = ","),"))")))

  rf_imp_summary <- dplyr::arrange(rf_imp_summary, desc(tot_n))
  rf_imp_summary$metric <- factor(rf_imp_summary$metric, levels = rf_imp_summary$metric, ordered = T)
  rf_imp_summary$tot_n <- NULL
  rf_imp_summary <- eval(parse(text = paste0("tidyr::gather(rf_imp_summary,\"rank\",\"count\",",paste0("c(",paste0("n",nrank, collapse = ","),"))")) ))
  rf_imp_summary$rank <- eval(parse(text=paste0("factor(rf_imp_summary$rank, levels = c(",paste0("\"","n",nrank,"\"" ,collapse = ","),"), ordered = T)")))



  rf_imp_summary_plot <- rf_imp_summary[rf_imp_summary$count > 0,]

  # Make plot
  p <- ggplot2::ggplot(rf_imp_summary_plot, aes(x = metric, y = count, fill = rank)) +
    geom_bar(stat = "identity", color = "black", position = position_stack(reverse = T)) +
    scale_fill_viridis_d() +
    #scale_fill_manual(breaks = c("n1", "n2", "n3", "n4", "n5"),labels = c("1", "2","3","4","5"), values = c("#00442A","#248443","#78C679","#D9EFA2","#FFFFE5")) +
    labs(fill = "Variable importance rank", x = "ALS metric", y = "Number of models") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "bottom")

  if (saveFigure) {
    pdf(file = file.path(outdir, paste0("varImpRank_top",nrank[length(nrank)],".pdf")),
        useDingbats = FALSE,
        height = 5,
        width = 7)
    print(p)
    dev.off()
  }

  return(list(rank = rf_imp_summary,plot = p))

}
