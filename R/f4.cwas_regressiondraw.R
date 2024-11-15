#' Perform Subset Regression Analysis for All Cell Types
#'
#' This function performs subset regression analysis for all cell types based on the provided data.
#'
#' @param subcelldata The subset of data containing sample IDs, cell types, cell rates, and phenotype groups.
#' @param sampleID The column name in `subcelldata` representing sample IDs (default is "Sample").
#' @param celltype The column name in `subcelldata` representing cell types (default is "Celltype").
#' @param cellrate The column name in `subcelldata` representing cell rates (default is "rate").
#' @param div_group The column name in `subcelldata` representing phenotype groups (default is "Phenotype").
#'
#' @return None
#'
#' @export
#'
#' @examples cwas_regressiondraw(cwas.test.data)
cwas_regressiondraw <- function(subcelldata,sampleID=c("Sample"),celltype=c("Celltype"),cellrate=c("rate"),div_group=c("Phenotype")){

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")


  df <- subcelldata[,c(sampleID,celltype,cellrate,div_group)]
  names(df) <- c("sampleID","celltype","cellrate","div_group")
  wide_df <- tidyr::spread(df, key = celltype, value = cellrate)
  wide_df[is.na(wide_df)] <- 0
  dd2 <- wide_df[,-1]
  #dd2 <-  dd2%>% dplyr::select_all() %>% dplyr::filter(div_group%in%c("controlcontrol","severe/criticalprogression"))
  #dd2$div_group <- ifelse(dd2$div_group=="controlcontrol",0,1)
  #allphno <- unique(wide_df$div_group)
  # 使用 ifelse() 函数进行条件判断和赋值
  dd2$div_group <- ifelse(dd2$div_group == "Control", 0, 1)

  sub.fit <- leaps::regsubsets(div_group ~ ., data = dd2)
  #sub.fit<- summary(sub.fit)
  best_summary <- summary(sub.fit)
  #par(mfrow = c(3, 2))

  # 定义一个函数来执行绘图并处理可能的错误
  plot_with_error_handling <- function(plot_function, ...) {
    tryCatch(
      {
        plot_function(...)
      },  error = function(e) {
        # 如果出现错误，则输出错误信息
        cat("An error occurred:", conditionMessage(e), "\n")
      }
    )
  }

  # 调用 plot_with_error_handling 函数执行绘图并处理错误
  plot_with_error_handling(plot, best_summary$adjr2, ylab = "adjr2", xlab = "number of features")
  plot_with_error_handling(plot, sub.fit, scale = "adjr2")

  plot_with_error_handling(plot, best_summary$cp, ylab = "cp", xlab = "number of features")
  plot_with_error_handling(plot, sub.fit, scale = "Cp")

  plot_with_error_handling(plot, best_summary$bic, ylab = "bic", xlab = "number of features")
  plot_with_error_handling(plot, sub.fit, scale = "bic")




}




