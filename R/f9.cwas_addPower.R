#' Add Power Information to Fitted Models
#'
#' This function adds power and confidence interval (CI) information to the given fitted models for analysis.
#'
#' @param mymodel A list of fitted models obtained from previous analyses.
#' @param yournsim Number of simulations to perform for power estimation (default is 10).
#' @param OR_cutoff_up Upper cutoff value for odds ratio (OR) (optional).
#' @param OR_cutoff_down Lower cutoff value for odds ratio (OR) (optional).
#' @param CCI Logical indicating whether to consider cell change index (CCI) (default is TRUE).
#' @param interaction Logical indicating whether to consider interaction effects (default is FALSE).
#' @param Tablename Name of the table to be created (default is "Tablename").
#'
#' @return A list of fitted models with added power and CI information.
#'
#' @export
#'
#' @examples cwas_addPower(cwas.test.data.single.model.best)
#'
cwas_addPower <- function(mymodel,yournsim=10,OR_cutoff_up=NULL,OR_cutoff_down=NULL,CCI=T,interaction=F,Tablename="Tablename"){

  if (!requireNamespace("simr", quietly = TRUE)) {
    stop("Package \"simr\" needed for this function to work.\n         Please install it by install.packages('simr')",
         call. = FALSE)
  }

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")
  #usethis::use_package(package="lme4",type="Import")
  #usethis::use_package(package="stats",type="Import")
  #usethis::use_package(package="ggrepel",type="Import")
  #usethis::use_package(package="gtsummary",type="Import")
  #usethis::use_package(package="grid",type="Import")
  #usethis::use_package(package="forestploter",type="Import")
  #usethis::use_package(package="magrittr",type="Import")
  #usethis::use_package(package="MuMIn",type="Import")
  #usethis::use_package(package="simr",type="Import")

  temp0 <-data.frame(0,0,0)
  colnames(temp0) <- c("powerLow","powerHigh","power")


  alltemp <- c()
  for (i in 1:length(mymodel[[1]])) {

    if(interaction){


      # 尝试执行混合模型
      tryCatch({
        temp <- parallel_power_calculation_v1(fit = mymodel[[1]][[i]],yournsim =yournsim,OR_cutoff_up = OR_cutoff_up,OR_cutoff_down = OR_cutoff_down,CCI=CCI,object="rate:rate2")

      }, error = function(e) {
        # 如果出现错误，提示无法完成混合模型计算,执行固定效应模型
        cat("Error: Mixed effect model was not applicable for this data.\n")
        temp <<- temp0
      })




    }else{

      # 尝试执行混合模型
      tryCatch({
        temp <- parallel_power_calculation_v1(fit = mymodel[[1]][[i]],yournsim =yournsim,OR_cutoff_up = OR_cutoff_up,OR_cutoff_down = OR_cutoff_down,CCI=CCI)

      }, error = function(e) {
        # 如果出现错误，提示无法完成混合模型计算,执行固定效应模型
        cat("Error: Mixed effect model was not applicable for this data.\n")
        temp <<- temp0
      })

    }

    print(i)
    alltemp <- rbind(alltemp,temp)
    print(i)


  }

  mymodel[[2]] <-  cbind(mymodel[[2]],alltemp)

  data <-  mymodel[[2]]

  data[c("powerHigh","power")]  <- lapply(data[c("powerHigh","power")], function(x) format(x, scientific = F, digits = 2))
  data$powerLow <- round(data$powerLow, 2)
  data$powerLow <- ifelse(data$powerLow<0.01,"0",data$powerLow)

  data$"Power(95%CI)" <- paste0(data$power, "(",data$powerLow,",",data$powerHigh,")")

  mymodel[[2]] <- data

  OR_plot(data = mymodel,show_variable=c("Celltype","Model","OR(95%CI)","P-value","P-value adjusted","Power(95%CI)"),plot.position=4,Tablename=Tablename)

  ### 画森林图


  return(mymodel)
}
