#' @title Prepare data for fitting
#'
#' @description  Parses data for the \code{pre_data} function.
#'    Validates if data is a summary data file & creates the tidy format.
#'    Tidy format is created by checking for a column name \code{variable}.
#'    Returns class of either type \code{"prepped_time"} or \code{""prepped_data"},
#'    depending if a \code{timeID} or \code{timeAfterExposure} column was found with multiple unique values.
#'
#'
#' @author Steven Wink
#'
#' @param xcol string, column name of independant variable (usually the concentration).
#'     Default: "dose_uM"
#' @param dataset dataframe containing the x, y and groups columns.
#'     Default: summary_data
#' @param groups string, column name of grouping variable which entries uniquely identify
#'     Default: "treatment"
#' the replicate set of concentration response curves for a certain condition.
#' @param minD integer, minimum number of unique concentrations
#' per concentration response curve, curves with less are removed
#'     Default: 7

#' @examples

#' parse_data(xcol = "dose_uM", groups = "treatment",
#'          minD = 7, dataset = test_data)
#'
#' @export



## TODO: maak een method voor het prepareren van time-lapse data voor de fm1 functie, en een
## method voor het prepareren van single time point data voor de fm1 functie



parse_data <- function(xcol = "dose_uM", groups = "treatment", minD = 7, dataset = summary_data) {

  # check if variable in column names, if so, make tidy
  if("variable" %in% colnames(dataset)){
    dataset <- dataset %>% tidyr::spread_(key = "variable", value = "value")
    message("Colum variable, tidy formatted table")
  }

  if(any(c("timeID", "timeAfterExposure") %in% colnames(dataset)) &
     (length(unique(dataset$time)) > 1 | length(unique(dataset$timeAfterExposure))) > 1) {

    class(dataset_list) <- append(class(dataset), "prepped_time")
    message("Found multiple time points, assigning class prepped_time")

    } else {

    class(dataset_list) <- append(class(dataset), "prepped_data")
    message("No multiple time points detected, assigning class prepped_data")
  }

}


#' @title generic for summarising data
#'
#' @description prepares \code{prepped_data} and \code{prepped_time} for model fitting.
#'    time summarisation methods are \code{max}, \code{auc} or \code{sel_time}
#'
#' @author Steven Wink

summarise_data <- function(input) {
  UseMethod("summarise_data", input)
}



summarise_data.prepped_data <- function(input) {
  dataset <- as.data.frame(
    dataset %>% dplyr::group_by_(groups, xcol) %>% dplyr::mutate(n=n()) # of replicates
  )
  dataset<-dataset[order(dataset[, groups], dataset[, xcol]),]

  dataset_list <- split(dataset, f = dataset[, groups]) # assuming 1 variable, 1 cell line, 1 time point
  ind_rm <- which(sapply(dataset_list, function(x) length(unique(x[, xcol]))) < (minD+1))
  dataset_list <- dataset_list[-ind_rm]


}


summarise_data.prepped_time <- function(input, method) {

  dataset <- as.data.frame(
    dataset %>% dplyr::group_by_(groups, xcol) %>% dplyr::mutate(n=n()) # of replicates
  )
  dataset<-dataset[order(dataset[, groups], dataset[, xcol]),]


  dataset_list <- split(dataset, f = dataset[, groups]) # assuming 1 variable, 1 cell line, 1 time point
  ind_rm <- which(sapply(dataset_list, function(x) length(unique(x[, xcol]))) < (minD+1))
  dataset_list <- dataset_list[-ind_rm]

}


#' @title fit loess regression
#'
#' @description dfdsf
#'
#' @author Steven Wink
#'
#' @param xcol string, column name of independant variable (usually the concentration).
#' @param ycol string, column name of dependent variable/ measure/ response variable.
#' @param dataset dataframe containing the x, y and groups columns.
#' @param respLev Numeric, A number representing threshold for point of departure,
#' strictly speaking this should be set to 0 for point of departures.
#' @param groups string, column name of grouping variable which entries uniquely identify
#' the replicate set of concentration response curves for a certain condition.
#' @param minD integer, minimum number of unique concentrations
#' per concentration response curve, curves with less are removed

#' @examples

#' model_pod(xcol = "GFP_int", ycol = "dose_uM", dataset = test_data,
#'       respL = 0.1,  groups = "treatment", minD = 7, span = 3/4, degree = 1)
#'
#' @export


model_pod <- function(xcol, ycol, dataset, respLev, groups, minD, ...) {


fm1_fun <- function(input, ...) {
  #input = dataset_list[[1]]
  fm1 <- loess(formula = get(ycol) ~ get(xcol), data = input, ...)


  conc.vec <- seq(min(input[, xcol]), max(input[, xcol]), length.out = 100)

  new_data <- data.frame( conc.vec )
  colnames(new_data) <- xcol
  predict(fm1, newdata = new_data, se = TRUE)

  prediction <- as.data.frame(
    predict(fm1, newdata =  new_data, se = TRUE )
  )



  list(model = fm1, conc.vec = conc.vec,
       prediction = prediction[, "fit"],
       se = prediction[, "se.fit"],
       n = input$n)
}

model_out <- lapply(dataset_list, fm1_fun, ... )





calc_pod <- function() {
  PoD_est = alist()
  TH_est = alist()
  for( i in seq_along(dataset_list)){

    ctrl_level <- model_out[[i]]$prediction[1]  # first predicted concentration level
    sd_regr <- model_out[[i]]$se * sqrt(model_out[[i]]$n[1]) # vector of regression standard errors * sqrt(n)
    TH_level = sd_regr + respLev + ctrl_level

    ind_predict <- model_out[[i]]$prediction > TH_level
    if(sum(ind_predict) > 0){
      ind_predict <- min(which(ind_predict))
    } else{
      ind_predict <- NA
    }
    if(!is.na(ind_predict)){
      PoD_est[[i]] <- model_out[[i]]$conc.vec[ind_predict]
      TH_est[[i]] <- model_out[[i]]$prediction[ind_predict]
    } else{
      PoD_est[[i]] <- NA
      TH_est[[i]] <- NA
    }

  }

  return(list(PoD_est = PoD_est, TH_est = TH_est))
}


est <- calc_pod()
TH_est <- est[["TH_est"]]
PoD_est <- est[["PoD_est"]]

write_plot <- function(dataset_list, model_out, PoD_est){

  if(!dir.exists("../results")){
    dir.create("../results")
  }
  pdf("../results/pod_fits_loess.pdf", height = 15, width = 10)
  par(mfrow = c(5,4))
  for(i in seq_along(dataset_list)){
    plot(log(dataset_list[[i]][, xcol]+1),dataset_list[[i]][, ycol], main = unique(dataset_list[[i]][, groups]
    ),
    xlab = paste("log(", xcol, "+1)"), ylab = ycol)
    lines(log(model_out[[i]]$conc.vec+1), model_out[[i]]$prediction)
    abline(v=log(PoD_est[[i]]+1), lwd = 3, lty = 1)
    abline(h=TH_est[[i]], lwd = 3, lty = 1)
  }
  dev.off()
}

write_plot(dataset_list, model_out, PoD_est)




result <- data.frame( treatment = unique(unlist(lapply(dataset_list, "[[", groups))),
                      PoD_est = unlist(PoD_est),
                      TH_est = unlist(TH_est)
)

#' @return dataframe with the point of departure estimates
result

}
#
# head(test_data)
#
# model_pod(xcol = "dose_uM", ycol = "GFP_int",
#           dataset = test_data, respLev = 0, groups = "treatment",
#           minD = 7, span = 3/4, degree = 1)


