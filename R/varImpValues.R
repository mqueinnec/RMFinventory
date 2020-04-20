#' Retrieve importance of models
#'
#' @param rfList update
#' @param measure update
#' @param saveFigure update
#' @param outdir update
#'
#' @export
#'
 varImpValues <- function(rfList,
                          measure = "IncMSE",
                          saveFigure = FALSE,
                          outdir) {

   attNames <- names(rfList)

   out <- list()

   for(n in attNames) {
     importance_df <- as.data.frame(randomForest::importance(rfList[[n]][["model"]][["finalModel"]]))
     importance_df$metrics <- row.names(importance_df)
     colnames(importance_df) <- c("IncMSE", "IncNodePurity", "metrics")
     importance_df <- eval(parse(text = paste0("dplyr::arrange(importance_df,", measure,")")))

     importance_df$metrics <- factor(importance_df$metrics, levels = importance_df$metrics, ordered = TRUE)

    importance_df$variable <- n

    out[[n]][["importance"]] <- importance_df

    p <- ggplot2::ggplot(importance_df, aes(x = IncMSE, y = metrics)) +
      geom_point() +
      ggtitle(titles[[n]]) +
      theme_bw()

    out[[n]][["plot"]] <- p

    if(saveFigure) {
      pdf(file = file.path(outdir,paste0("importance_",n,".pdf")),
          height = 5,
          width = 7)
      print(p)
      dev.off()
    }


   }

   return(out)
 }
