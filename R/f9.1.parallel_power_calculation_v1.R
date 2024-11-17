parallel_power_calculation_v1 <- function(fit,yournsim=10,OR_cutoff_up=NULL,OR_cutoff_down=NULL,object="rate",CCI=T){
  #library(simr)
  environment(attributes(fit@frame)$formula) <- globalenv()
  cell2df <<-fit@frame
  #celldf <<-fit@frame
  # 计算均值差值
  result_diff <- aggregate(rate ~ Phenotype, data = fit@frame, FUN = mean)
  colnames(result_diff) <- c("Phenotype", "mean_rate")  # 重命名列名

  # 是否考虑细胞变化量
  ifelse(CCI, mean_difference <<- diff(result_diff$mean_rate),mean_difference<<-1)



  if(!is.null(OR_cutoff_up)){

    if(fixef(fit)[object]>0){
      fixef(fit)[object] <- log(OR_cutoff_up)/mean_difference #指定固定效应x为-0.1

      print(paste(object,"risk were set as",OR_cutoff_up))
    }else{
      fixef(fit)[object] <- log(OR_cutoff_down)/mean_difference #指定固定效应x为-0.1
      print(paste(object,"risk were set as",OR_cutoff_down))
    }


  }

  power <- simr::powerSim(fit,nsim=yournsim)
  pow <- stats::confint(power)
  df <- data.frame(pow)
  # df$powerLow <- paste0(pow[, 1] * 100, "%")
  # df$powerHigh <- paste0(pow[, 2] * 100, "%")
  df$powerLow <- pow[, 1]
  df$powerHigh <- pow[, 2]
  df$power <- (pow[, 1] + pow[, 2]) / 2
  # df$power <- paste0(df$power * 100, "%")
  return(df[, 3:5])

}
