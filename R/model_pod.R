#' fit piece wise splines
#'
#' @param xcol string, column name of dependant variable.
#' @param ycol string, column name of independent variable/ measure variable.
#' @param dataset dataframe containing the x, y and groups columns.
#' @param respLev Numeric, A number representing threshold for point of departure.
#' @param groups string, column name of grouping variable.
#'
#' @examples

#' model_pod(xcol = "GFP_int", ycol = "dose_uM", dataset = test_data,
#'       respL = 0.1,  groups = "treatment")
#'
#' @export
model_pod <- function(xcol, ycol, dataset, respLev, groups) {

dataset<-dataset[order(dataset[, groups], dataset[, xcol]),]

dataset_list <- split(dataset, f = dataset[, groups]) # assuming 1 variable, 1 cell line, 1 time point
ind_rm <- which(sapply(dataset_list, function(x) length(unique(x[, xcol]))) < 7)
dataset_list <- dataset_list[-ind_rm] #

fm1_fun <- function(input) {

  fm1 <- lm(input[,ycol] ~ splines::ns(input[, xcol], df = 3))
  conc.vec <- seq(min(input[, xcol]), max(input[, xcol]), length.out = nrow(input))
  new_data <- data.frame( conc.vec )
  colnames(new_data) <- xcol
  prediction <- as.data.frame(
    predict(fm1, newdata =  new_data, interval = "confidence" )

  )

  list(model = fm1, conc.vec = conc.vec, prediction = prediction[, "fit"], upr = prediction[,"upr"])

}


model_out <- lapply(dataset_list, fm1_fun )





# calculation POD, defined as 2xSD above control sample mean With SD the combined deviance and control SD


#control_model <- lm(GFP_pos ~ 1  , data = controlData) #
#control_est_plus_interval <- confint(control_model, level = 0.9)[2]

# POD_est POD estimate, mean


# use confidence interval of control and fit to calculate pod
calc_pod <- function(dataset_list) {
PoD_est = alist()
  for( i in seq_along(dataset_list)){
  #square root of the mean of the squared residual values
  #sd_total <- sqrt( (sqrt(mean(model_out[[i]]$model$residuals^2)))^2 + control_sd^2 ) #total standard deviation
  #TH_level <- control_est + respLev  #control mean plus 2x the standard deviation

  cl_regr <- model_out[[i]]$upr - model_out[[i]]$prediction
  TH_level = cl_regr + respLev

  ind_predict <- model_out[[i]]$prediction > TH_level
  if(sum(ind_predict) > 0){
    ind_predict <- min(which(ind_predict))
  } else{
    ind_predict <- NA
  }
  if(!is.na(ind_predict)){
    PoD_est[[i]] <- model_out[[i]]$conc.vec[ind_predict]
  } else{
    PoD_est[[i]] <- NA
  }

}

  return(PoD_est)
}


PoD_est <- calc_pod(dataset_list)


write_plot <- function(dataset_list, model_out, PoD_est){

  if(!dir.exists("../results")){
  dir.create("../results")
}
pdf("../results/pod_fits.pdf", height = 15, width = 10)
par(mfrow = c(5,4))
for(i in seq_along(dataset_list)){
  plot(log(dataset_list[[i]][, xcol]),dataset_list[[i]][, ycol], main = unique(dataset_list[[i]][, groups]
  ), ylim = c(-0.1,1.1),
  xlab = xcol, ylab = ycol)
  lines(log(model_out[[i]]$conc.vec), model_out[[i]]$prediction)
  abline(v=log(PoD_est[[i]]), lwd = 3, lty = 1)



}
dev.off()
}

write_plot(dataset_list, model_out, PoD_est)




result <- data.frame( treatment = unique(unlist(lapply(dataset_list, "[[", groups))),
                      PoD_est = unlist(PoD_est)
                      )

#' @return dataframe with the point of departure estimates
result

}
# model_pod(xcol = "dose_uM", ycol = "GFP_int", dataset = test_data, respLev = 0.1, groups = "treatment" )
