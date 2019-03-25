

model_pod <- function(x, y, dataset, respL, groups) {

#respL = 0.15 # absolute difference between control and response PoD
#dataset <- summarised_data



# fit piece wise splines
dataset<-dataset[order(dataset$treatment, dataset$dose_uM),]

dataset_list <- split(dataset, f = dataset[, groups]) # assuming 1 variable, 1 cell line, 1 time point
ind_rm <- which(sapply(dataset_list, function(x) length(unique(x$dose_uM))) < 7)
dataset_list <- dataset_list[-ind_rm] #


fm1_fun <- function(input) {

  fm1 <- lm(GFP_pos ~ splines::ns(dose_uM, df = 3), data = input)
  conc.vec <- with(input, seq(min(dose_uM), max(dose_uM), length.out = 100))
  prediction <- as.data.frame(
    predict(fm1, newdata =  data.frame(dose_uM = conc.vec), interval = "confidence" )
  )

  list(model = fm1, conc.vec = conc.vec, prediction = prediction[, "fit"], upr = prediction[,"upr"])

}

model_out <- lapply(dataset_list, fm1_fun )

# calculation POD, defined as 2xSD above control sample mean With SD the combined deviance and control SD

controlData <- dataset[dataset$treatment == "DMSO", ]

controlData$replID <- as.factor(controlData$replID)

#control_model <- lm(GFP_pos ~ 1  , data = controlData) #
#control_est_plus_interval <- confint(control_model, level = 0.9)[2]

# POD_est POD estimate, mean
PoD_est = alist()

# use confidence interval of control and fit to calculate pod

for( i in seq_along(dataset_list)){
  #square root of the mean of the squared residual values
  #sd_total <- sqrt( (sqrt(mean(model_out[[i]]$model$residuals^2)))^2 + control_sd^2 ) #total standard deviation
  #TH_level <- control_est + respL  #control mean plus 2x the standard deviation

  cl_regr <- model_out[[i]]$upr - model_out[[i]]$prediction
  TH_level = cl_regr + respL

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


if(!dir.exists("../results")){
  dir.create("../results")
}
pdf("../results/pod_fits.pdf", height = 15, width = 10)
par(mfrow = c(5,4))
for(i in seq_along(dataset_list)){
  plot(log(dataset_list[[i]]$dose_uM),dataset_list[[i]]$GFP_pos, main = unique(dataset_list[[i]]$treatment
  ), ylim = c(-0.1,1.1),
  xlab = "Dose", ylab = "GFP_pos")
  lines(log(model_out[[i]]$conc.vec), model_out[[i]]$prediction)
  abline(v=log(PoD_est[[i]]), lwd = 3, lty = 1)



}


dev.off()


result <- data.frame( treatment = unique(unlist(lapply(dataset_list, "[[", "treatment"))),
                      PoD_est = unlist(PoD_est)
                      )

write.table(result, file = "../results/PoD_data.txt", sep ="\t", col.names = NA)

}
