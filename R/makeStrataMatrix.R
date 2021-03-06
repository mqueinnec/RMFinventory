#' Matrix of strata
#'
#' This function returns a matrix of strata base on the number of breaks
#'
#' @param nbreaks Number of breaks
#' @param all_strata Data frame with all possible strata (breaks and ID)
#'
#' @export
#'
makeStrataMatrix <- function(nbreaks, all_strata) {

  all_strata <- dplyr::arrange(all_strata, names(all_strata)[1], names(all_strata)[1])

  matrix_strata <- matrix(all_strata$strata,nr=nbreaks[2]-1,nc=nbreaks[1]-1, byrow = FALSE)
  # matrix_strata <- matrix(NA,nr=nbreaks[2]-1,nc=nbreaks[1]-1)
  # for(i in 1:(nbreaks[1]-1)){
  #   for(j in 1:(nbreaks[2]-1)){
  #     matrix_strata[j,i] <- all_strata[i,j]
  #       as.numeric(sprintf("%i%i",i,j))
  #   }
  # }
  return(matrix_strata)
  }

#' Neighbouring strata
#'
#' @param matrix A matrix
#' @param strata Strata number
#'
#' @export

getNeiMatrix <- function(matrix, strata) {

  ind_strata <- which(matrix==strata,arr.ind=T)
  if (ind_strata[1]==dim(matrix)[1]){ #If on edge
    select_rows <- c((ind_strata[1]-1):(ind_strata[1]))
  }else if (ind_strata[1] == 1){ #If on edge
    select_rows <- c((ind_strata[1]):(ind_strata[1]+1))
  }else{
    select_rows <- c((ind_strata[1]-1):(ind_strata[1]+1))
  }

  if (ind_strata[2]==dim(matrix)[2]){ #If on edge
    select_cols <- c((ind_strata[2]-1):(ind_strata[2]))
  }else if (ind_strata[1] == 1){ #If on edge
    select_cols <- c((ind_strata[2]):(ind_strata[2]+1))
  }else{
    select_cols <- c((ind_strata[2]-1):(ind_strata[2]+1))
  }

  as.vector(matrix[select_rows,select_cols])
}
