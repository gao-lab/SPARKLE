#' Generate a heatmap of AIC values for multiple models
#'
#' This function takes a list of model dataframes and generates a heatmap
#' to visualize the AIC values across different models.
#'
#' @param mdlist A list of model dataframes. Each element in the list corresponds
#'               to a different cell type or model type.
#' @param method Character string specifying the method to use for ranking.
#'               Currently supports "rank" (default) and "AIC".
#'
#' @details
#' The function processes each model dataframe in the input list (\code{mdlist}),
#' extracts relevant columns (formula, AIC, Pvalue, rank), and combines them into
#' a single dataframe (\code{alldata}). It then cleans unnecessary columns,
#' orders the data, and prepares it for heatmap plotting. P-values are used to
#' annotate the heatmap cells based on significance levels.
#'
#' @examples
#' # Example usage:
#' # cwas_allmodel_heatmap(list(model_data1, model_data2))
#'
#' @import pheatmap
#' @export
#'

cwas_allmodel_heatmap <- function (mdlist,method=c("rank")) {
  i = 1
  celltype <- names(mdlist)
  mydata <- mdlist[[i]][[2]]
  mydata$rank <- seq(1:length(mydata$AIC))
  mydata <- mydata[, c("formula", "AIC", "Pvalue", "rank")]
  colnames(mydata) <- paste0(celltype[1], ".", colnames(mydata))
  mydata$name <- mydata[, 1]
  alldata <- mydata
  for (i in 2:length(mdlist)) {
    mydata <- mdlist[[i]][[2]]
    mydata$rank <- seq(1:length(mydata$AIC))
    mydata <- mydata[, c("formula", "AIC", "Pvalue", "rank")]
    colnames(mydata) <- paste0(celltype[i], ".", colnames(mydata))
    mydata$name <- mydata[, 1]
    alldata <- merge(alldata, mydata, by = "name", all = T)
    #print(i)
  }
  alldata_cleaned <- alldata[, !grepl("celltype", names(alldata))]
  alldata_cleaned <- alldata_cleaned[, !grepl("model", names(alldata_cleaned))]
  alldata_cleaned <- alldata_cleaned[, !grepl("formula", names(alldata_cleaned))]
  alldata_cleaned <- alldata_cleaned[, !grepl("id", names(alldata_cleaned))]
  alldata_cleaned$formula <- alldata[, 2]
  alldata_cleaned$modelID <- alldata[, 3]
  alldata_cleaned <- alldata_cleaned[order(alldata_cleaned$modelID),
  ]
  alldata <- alldata[order(alldata[, 3]), ]
  allformula <- alldata$name
  pvalue_columns <- grep("Pvalue", names(alldata_cleaned),
                         value = TRUE)
  AICrank_columns <- grep(method, names(alldata_cleaned), value = TRUE)
  pvalue_data <- alldata_cleaned[c("name", pvalue_columns)]
  heatmap_data <- alldata_cleaned[c("name", AICrank_columns)]

  ## 绘制热图

  heatmap_data2 <- heatmap_data
  heatmap_data2[, -1] <- lapply(heatmap_data2[, -1], as.numeric)
  max_value <- max(heatmap_data2[, -1], na.rm = TRUE)

  # 将 heatmap_data2 中的所有 NA 值替换成最大值
  heatmap_data2[is.na(heatmap_data2)] <- max_value
  heatmap_data2$name <- gsub("binomial link function: ", "", heatmap_data2$name)
  row.names(heatmap_data2) <- heatmap_data2$name
  colnames(heatmap_data2) <- gsub(".rank", "", colnames(heatmap_data2))
  colnames(heatmap_data2) <- gsub(".AIC", "", colnames(heatmap_data2))
  heatmap_data2$name=NULL

  pvalue_matrix <- pvalue_data
  pvalue_matrix$name=NULL
  rownames(pvalue_matrix) <- rownames(heatmap_data2)
  colnames(pvalue_matrix) <- colnames(heatmap_data2)

  # 定义p值注释
  annotation_matrix <- ifelse(pvalue_matrix < 0.001, "***",
                              ifelse(pvalue_matrix < 0.01, "**",
                                     ifelse(pvalue_matrix < 0.05, "*", "")))
  annotation_matrix[is.na(annotation_matrix)] <- "NA"


  pheatmap::pheatmap(heatmap_data2,
                     display_numbers = annotation_matrix,
                     color = colorRampPalette(c("red", "white"))(100),
                     cluster_rows = T,
                     cluster_cols = T,
                     angle_col = 90,
                     fontsize_row = 10,
                     fontsize_col = 10,
                     annotation_colors = list(display_numbers = c("black")),
                     main = paste0("Heatmap of AIC Values for All Models","(",method,")")
  )

}

