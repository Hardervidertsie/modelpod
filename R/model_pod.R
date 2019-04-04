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
#' @param dataset dataframe containing the x, y and Groups columns.
#'     Default: summary_data
#' @param Groups string, column name of grouping variable which entries uniquely identify
#'     Default: "treatment"
#' the replicate set of concentration response curves for a certain condition.
#' @param minD integer, minimum number of unique concentrations
#' per concentration response curve, curves with less are removed
#'     Default: 7

#' @examples

#' parse_data(xcol = "dose_uM", Groups = "treatment",
#'          minD = 7, dataset = test_data)
#'
#' @export





# S3 class constructor for pod
new_pod <- function(x = data.frame(), type = "poddata") {
  stopifnot(is.data.frame(x))
  type = match.arg(type, c("poddata", "podtime"))

  structure(x,
            class = append(class(x), type)
            )
}

# pod validator function
validate_pod <- function(x) {
  if(inherits(x, "pod")){
    stop("Not class pod, please call parse_data() first ")
  }

  if(nrow(x) < 1){
    stop("No rows in data")
  }

    x
}


# parse_data: set class based on time lapse or no-time lapse

parse_data <- function(dataset ) {

  message("A summary data file:\n
          replID entries are biological replicates and plateID entries are unique
          for each individual plate.
          So for a dataset with 3 biological replicates each plateID should pair with 3 replID entries
          missing replicates are ok, for example:

          plateID:  replID:\n
          plate1    rep1
          plate2    rep2
          plate3    rep3
          plate4    rep1
          plate5    rep2
          plate6    rep1
          plate7    rep2
          plate8    rep3")
  # check if input is a summary data file with minimal set of required columns
  min_columnset <- c("treatment", "replID", "plateID", "dose_uM", "cell_line")
  ind_ok <- min_columnset %in% colnames(dataset)

  if (!all(ind_ok)){
    stop(message(paste("missing columns: ", min_columnset[!ind_ok], collapse = "\n")))
  }

  # check if variable in column names, if so, make tidy
  if ("variable" %in% colnames(dataset)){
    dataset <- dataset %>% tidyr::spread_(key = "variable", value = "value")
    message("Colum variable, tidy formatted table")
  }

  if (any(c("timeID", "timeAfterExposure") %in% colnames(dataset)) &
     (length(unique(dataset$time)) > 1 | length(unique(dataset$timeAfterExposure)) > 1)) {

    dataset <- new_pod(dataset, type = "podtime")
    message("Found multiple time points, assigning class podtime")

    } else {
      dataset <- new_pod(dataset, type = "poddata")
    message("No multiple time points detected, assigning class poddata")
  }
return(dataset)
}

#load("data/test_data.rda")
summary_data <- parse_data(dataset = summary_data)

head(summary_data)

#' @title generic for summarising data
#'
#' @description prepares \code{poddata} and \code{podtime} for model fitting.
#'    time summarisation methods are \code{max}, \code{auc} or \code{sel_time}
#'
#' @author Steven Wink

summarise_data <- function(dataset, ... ) {

   UseMethod("summarise_data")

  }


summarise_data.poddata <- function(dataset, Groups, xcol , method ) {

  Groups <- enquo(Groups)
  xcol <- enquo(xcol)

calc_tech_rep <- function(dataset, Groups, xcol ) {


   tmp <- dataset %>%
    dplyr::group_by(!! Groups, replID, !! xcol) %>%
    dplyr::mutate(all_rep = n()) %>% #of technical? replicates
    ungroup() %>%
    dplyr::group_by(!! Groups, !! xcol) %>%
    dplyr::mutate(bio_rep =  n_distinct(replID)) %>%
    mutate(tech_rep = all_rep > bio_rep)
if (any(tmp$tech_rep)){
    return(TRUE)
 } else {
    return(FALSE)
 }
}

tech_rep <- calc_tech_rep(dataset, Groups, xcol)
tech_rep


# if technical reps: send message, then calculate summary over Groups + reps + xcol

if (tech_rep) {

min_columnset <- c("treatment", "replID", "plateID", "dose_uM", "cell_line")

message( "Technical replicates detected, calculating mean over: \n")

message( paste(c(
                 min_columnset,
                 paste(quo_name(Groups), "(Groups)")
                 ),
                 collapse = "\n"
                 )
              )

  calc_summary <- function(dataset, min_columnset, Groups) {

  groups_var <- rlang::syms(min_columnset)

  dataset %>% dplyr::group_by(!!! groups_var, !! Groups) %>%
  dplyr::summarise_all(.funs = mean)
  }

  dataset <- calc_summary(dataset, min_columnset, Groups )

 }

class(dataset) <- append(class(dataset), "pod")
attributes(dataset)$sum_vars <- sum_vars
attributes(dataset)$all_vars_set <- all_vars_set
attributes(dataset)$Groups <- quo_text(Groups)

 return(dataset)

}

summarise_data.podtime <- function(dataset, Groups, xcol, method  ) {

  Groups <- rlang::enquo(Groups)
  xcol <- rlang::enquo(xcol)
  method = rlang::enquo(method)

  method <- match.arg(rlang::quo_text(method), choices = c("auc", "max", "mean"))

  #  calculate summary over time:

  min_columnset <- c("treatment", "replID", "plateID", "dose_uM", "cell_line")
  all_vars <- colnames(dataset)[!colnames(dataset) %in% c("timeID", "timeAfterExposure")]
  all_vars_set <- intersect(all_vars, min_columnset)

  sum_vars <- all_vars[ !all_vars %in% all_vars_set]


  message("summarising time lapse data per:\n")
  message(paste( all_vars_set, collapse = "\n"))
  message("   ========   ")
  message("calculating summary over:\n")
  message(paste(sum_vars, collapse = "\n"))


  calc_summary_t <- function(dataset, all_vars_set, method) {
    all_vars_set <- rlang::syms(all_vars_set)
    dataset %>% dplyr::group_by(!!! all_vars_set) %>%
        dplyr::summarise_at(vars(sum_vars), .funs = method)
    }


  auc_fun <- function() {
    timecol <-  dplyr::if_else("timeID" %in% colnames(dataset),  "timeID", "timeAfterExposure")
    dataset %>% dplyr::group_by(!!! rlang::syms(all_vars_set)) %>%
      dplyr::summarise_at(vars(sum_vars),
                          funs(MESS::auc(x = !! rlang::sym(timecol), y = .  , type = "spline")))
  }

if (method == "auc"){
  dataset <- auc_fun()
  } else {
  dataset <- calc_summary_t(dataset, all_vars_set, method)
  }


  class(dataset) <- append(class(dataset), "pod")
  attributes(dataset)$sum_vars <- sum_vars
  attributes(dataset)$all_vars_set <- all_vars_set
  attributes(dataset)$Groups <- quo_text(Groups)

  return(dataset)

}

summarized_data <- summarise_data(dataset = summary_data, Groups = treatment,  xcol = dose_uM, method = auc  )

attributes(summarized_data)



#' @title fit loess regression
#'
#' @description dfdsf
#'
#' @author Steven Wink
#'
#' @param xcol string, column name of independant variable (usually the concentration).
#' @param ycol string, column name of dependent variable/ measure/ response variable.
#' @param dataset dataframe containing the x, y and Groups columns.
#' @param respLev Numeric, A number representing threshold for point of departure,
#' strictly speaking this should be set to 0 for point of departures.
#' @param Groups string, column name of grouping variable which entries uniquely identify
#' the replicate set of concentration response curves for a certain condition.
#' @param minD integer, minimum number of unique concentrations
#' per concentration response curve, curves with less are removed

#' @examples

#' model_pod(xcol = "GFP_int", ycol = "dose_uM", dataset = test_data,
#'       respL = 0.1,  Groups = "treatment", minD = 7, span = 3/4, degree = 1)
#'
#' @export


model_pod <- function(dataset, ...) {
  UseMethod("model_pod")
}


model_pod.pod <- function(dataset = summarized_data, xcol = "dose_uM", ycol = "GFP_int", respLev = 0, minD = 7, ...) {

  all_vars_set <- attributes(dataset)$all_vars_set
  sum_vars <- attributes(dataset)$sum_vars
  Groups <- attributes(dataset)$Groups
  xcol <- rlang::quo_text(rlang::enquo(xcol))
  ycol <- rlang::quo_text(rlang::enquo(ycol))

  # calc number of replicates

  dataset <- dataset %>%
    dplyr::group_by(!! rlang::sym(Groups), !! rlang::sym(xcol)) %>%
    dplyr::mutate(n = n())


dataset_list <- split(dataset, f = dataset[ , Groups])
rmind <- which(lapply(dataset_list, nrow) < (minD + 1))
if(!length(rmind) == 0) {
dataset_list <- dataset_list[-rmind]
}


fm1_fun <- function(input, ...) {

  fm1 <- stats::loess(formula = get(ycol) ~ get(xcol), data = input)

  conc.vec <- seq(min(input[, xcol]), max(input[, xcol]), length.out = 100)

  new_data <- data.frame( conc.vec )
  colnames(new_data) <- xcol

  prediction <- as.data.frame(
    predict(fm1, newdata =  new_data, se = TRUE )
  )




  list(model = fm1, conc.vec = conc.vec,
        prediction = prediction[, "fit"],
        se = prediction[, "se.fit"],
        n = input$n)

  }

model_out <- lapply(dataset_list, fm1_fun )

class(model_out) <- append(class(model_out), "podmodel")

return(model_out)


} # ... eg span, degree


model_data <- model_pod(summarized_data, dose_uM, GFP_int, respLev = 0, minD = 7, degree = 1, span = 3/4)


calc_pod <- function(dataset, ...) {
  UseMethod("calc_pod")
}

calc_pod.podmodel <- function(model_data) {
  PoD_est = alist()
  TH_est = alist()

  for( i in seq_along(model_data)){

    ctrl_level <- model_data[[i]]$prediction[1]  # first predicted concentration level
    sd_regr <- model_data[[i]]$se * sqrt(model_data[[i]]$n[1]) # vector of regression standard errors * sqrt(n)
    TH_level = sd_regr + respLev + ctrl_level

    ind_predict <- model_data[[i]]$prediction > TH_level
    if(sum(ind_predict) > 0){
      ind_predict <- min(which(ind_predict))
    } else{
      ind_predict <- NA
    }
    if(!is.na(ind_predict)){
      PoD_est[[i]] <- model_data[[i]]$conc.vec[ind_predict]
      TH_est[[i]] <- model_data[[i]]$prediction[ind_predict]
    } else{
      PoD_est[[i]] <- NA
      TH_est[[i]] <- NA
    }

  }



  result <- data.frame( Groups = unlist(names(model_data)),
                        PoD_est = unlist(PoD_est),
                        TH_est = unlist(TH_est)


  )


  class(result) <- append(class(result), "podest")
  return(result)
}



pod_result <- calc_pod(model_data)



## todo next: write_plot
write_plot <- function(){

  if(!dir.exists("../results")){
    dir.create("../results")
  }
  pdf("../results/pod_fits_loess.pdf", height = 15, width = 10)
  par(mfrow = c(5,4))
  for(i in seq_along(dataset_list)){
    plot(log(dataset_list[[i]][, xcol]+1),dataset_list[[i]][, ycol], main = unique(dataset_list[[i]][, Groups]
    ),
    xlab = paste("log(", xcol, "+1)"), ylab = ycol)
    lines(log(model_out[[i]]$conc.vec+1), model_out[[i]]$prediction)
    abline(v=log(PoD_est[[i]]+1), lwd = 3, lty = 1)
    abline(h=TH_est[[i]], lwd = 3, lty = 1)
  }
  dev.off()
}

write_plot(dataset_list, model_out, PoD_est)




#
# head(test_data)
#
# model_pod(xcol = "dose_uM", ycol = "GFP_int",
#           dataset = test_data, respLev = 0, Groups = "treatment",
#           minD = 7, span = 3/4, degree = 1)


