#' Plot Heatmap of AIC Values and P-values for All Models
#'
#' This function generates a heatmap displaying AIC (Akaike Information Criterion) values and P-values for all models.
#'
#' @param mdlist A list of models containing AIC values and P-values.
#'
#' @export
#'
#' @examples cwas_model_heatmap(cwas.test.data.single.model)
#'
cwas_model_heatmap <- function(mdlist){

  #library(ggplot2)
  #library(reshape2)

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")
  #usethis::use_package(package="lme4",type="Import")
  #usethis::use_package(package="stats",type="Import")
  #usethis::use_package(package="ggrepel",type="Import")
  #usethis::use_package(package="reshape2",type="Import")



  i=1
  celltype <- names(mdlist)

  mydata <- mdlist[[i]][[2]]
  mydata$rank <- seq(1:length(mydata$AIC))
  mydata <- mydata[,c("formula","AIC","Pvalue","rank")]

  colnames(mydata) <- paste0(celltype[1],".",colnames(mydata))
  mydata$name <- mydata[,1]

  alldata <- mydata


  for(i in 2:length(mdlist)){

    mydata <- mdlist[[i]][[2]]
    mydata$rank <- seq(1:length(mydata$AIC))
    mydata <- mydata[,c("formula","AIC","Pvalue","rank")]
    colnames(mydata) <- paste0(celltype[i],".",colnames(mydata))
    mydata$name <- mydata[,1]
    alldata <- merge(alldata,mydata,by="name",all=T)
    #plist[[i]] <-plot_aic_pvalue(mydata)
    print(i)
  }


  alldata_cleaned <- alldata[, !grepl("celltype", names(alldata))]
  alldata_cleaned <- alldata_cleaned[, !grepl("model", names(alldata_cleaned))]
  alldata_cleaned <- alldata_cleaned[, !grepl("formula", names(alldata_cleaned))]
  alldata_cleaned <- alldata_cleaned[, !grepl("id", names(alldata_cleaned))]
  alldata_cleaned$formula <- alldata[,2]#formula
  alldata_cleaned$modelID <- alldata[,3]#id

  alldata_cleaned <- alldata_cleaned[order(alldata_cleaned$modelID),]

  alldata <- alldata[order(alldata[,3]),]
  allformula <- alldata$name



  pvalue_columns <- grep("Pvalue", names(alldata_cleaned), value = TRUE)
  AICrank_columns <- grep("AIC", names(alldata_cleaned), value = TRUE)

  # 创建一个新的数据框，仅包含上述列和 name 列
  pvalue_data <- alldata_cleaned[c("name", pvalue_columns)]
  heatmap_data <- alldata_cleaned[c("name", AICrank_columns)]

  # 使用 reshape2 包的 melt 函数将数据转换为适合绘制热图的格式
  melted_data <- reshape2::melt(heatmap_data, id.vars = "name", variable.name = "variable", value.name = "value")
  pvalue_melted_data <- reshape2::melt(pvalue_data, id.vars = "name", variable.name = "variable", value.name = "value")
  pvalue_melted_data$variable <- gsub("\\..*", "", pvalue_melted_data$variable)
  melted_data$variable <- gsub("\\..*", "", melted_data$variable)

  #melted_data$name <- factor(melted_data$name, levels = rev(unique(melted_data$name)))
  library(scales)


  p1 <-ggplot2::ggplot(melted_data, ggplot2::aes(x = variable, y = name, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "red", high = "white") +  # 设置颜色渐变
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),  # 旋转 x 轴文本
          axis.text.y = ggplot2::element_text(angle = 0, hjust = 0)) +  # y 轴标签水平显示
    #scale_y_discrete(labels = function(y) paste("model", y, allformula)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),  # 去除背景
          panel.grid.major = ggplot2::element_blank(),  # 去除主要网格线
          panel.grid.minor = ggplot2::element_blank(),  # 去除次要网格线
          axis.text.y.right = ggplot2::element_text(angle = 0, hjust = 1)) +  # 右侧 y 轴文本水平显示
    ggplot2::labs(x = "celltype", y = "model") # y轴标签放在右侧

  # 根据 Pvalue 添加显著性标记
  # 删除包含 NA 值的行
  pvalue_melted_data <- pvalue_melted_data[complete.cases(pvalue_melted_data$value), ]
  print(p1)
  p1 <- p1 +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.001, ],
                       ggplot2::aes(label = "***"), color = "black") +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.01 & pvalue_melted_data$value >= 0.001, ],
                       ggplot2::aes(label = "**"), color = "black") +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.05 & pvalue_melted_data$value >= 0.01, ],
                       ggplot2::aes(label = "*"), color = "black")

  # 打印热图
  print(p1)



  # 提取所有包含 "Pvalue" 的列
  pvalue_columns <- grep("Pvalue", names(alldata_cleaned), value = TRUE)
  AICrank_columns <- grep("rank", names(alldata_cleaned), value = TRUE)

  # 创建一个新的数据框，仅包含上述列和 name 列
  pvalue_data <- alldata_cleaned[c("name", pvalue_columns)]
  heatmap_data <- alldata_cleaned[c("name", AICrank_columns)]

  # 使用 reshape2 包的 melt 函数将数据转换为适合绘制热图的格式
  melted_data <- reshape2::melt(heatmap_data, id.vars = "name", variable.name = "variable", value.name = "value")
  pvalue_melted_data <- reshape2::melt(pvalue_data, id.vars = "name", variable.name = "variable", value.name = "value")
  pvalue_melted_data$variable <- gsub("\\..*", "", pvalue_melted_data$variable)
  melted_data$variable <- gsub("\\..*", "", melted_data$variable)

  #melted_data$name <- factor(melted_data$name, levels = rev(unique(melted_data$name)))
  library(scales)

  #unique(pvalue_melted_data$name)

  p2 <- ggplot2::ggplot(melted_data,  ggplot2::aes(x = variable, y = name, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "red", high = "white", trans = "sqrt") +  # 设置颜色渐变
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),  # 旋转 x 轴文本
          axis.text.y =  ggplot2::element_text(angle = 0, hjust = 0)) +  # y 轴标签水平显示
    #scale_y_discrete(labels = function(y) paste("model", y, allformula)) +
    ggplot2::theme(panel.background =  ggplot2::element_blank(),  # 去除背景
          panel.grid.major = ggplot2::element_blank(),  # 去除主要网格线
          panel.grid.minor = ggplot2::element_blank(),  # 去除次要网格线
          axis.text.y.right = ggplot2::element_text(angle = 0, hjust = 1)) +  # 右侧 y 轴文本水平显示
    ggplot2::labs(x = "celltype", y = "model") # y轴标签放在右侧

  # 根据 Pvalue 添加显著性标记
  # 删除包含 NA 值的行
  pvalue_melted_data <- pvalue_melted_data[complete.cases(pvalue_melted_data$value), ]
  print(p2)
  p2 <- p2 +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.001, ],
                       ggplot2::aes(label = "***"), color = "black") +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.01 & pvalue_melted_data$value >= 0.001, ],
                       ggplot2::aes(label = "**"), color = "black") +
    ggplot2::geom_text(data = pvalue_melted_data[pvalue_melted_data$value < 0.05 & pvalue_melted_data$value >= 0.01, ],
                       ggplot2::aes(label = "*"), color = "black")

  # 打印热图
  print(p2)

}
