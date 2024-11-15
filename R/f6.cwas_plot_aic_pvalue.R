#' Plot AIC and P-value Scatter Plots for Each Cell Type
#'
#' This function generates separate scatter plots of AIC (Akaike Information Criterion) values against P-values for each cell type.
#'
#' @param mdlist A list of models containing AIC values and P-values for each cell type.
#' @param color Optional. A vector of color codes to use for differentiating models in the scatter plots.
#'
#' @return A list of AIC and P-value Scatter Plots
#' @export
#'
#' @examples cwas_plot_aic_pvalue(cwas.test.data.single.model)
#'
cwas_plot_aic_pvalue <- function(mdlist,color=NULL){

  if(is.null(color)){

    default_color <- c(
      "#396AB1", "#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#8B3D47", "#196D96", "#A3A3A3", "#FF9309", "#8B9BCA", "#A3A3A3", "#67BF5C", "#B5B5B5", "#D8B83F", "#D86375", "#6ED641", "#A3A3A3", "#D89841", "#8CACC6", "#DB6D00", "#97DFFC", "#B5B5B5", "#AB9F45", "#94BFA2", "#E17C05", "#71588F", "#DE6A61", "#6D6E71", "#D8B83F", "#9278A8", "#969696", "#000000", "#B5B5B5", "#D88E2D", "#76A233", "#6B4C9A", "#6D9CCF", "#8CACC6", "#8B9BCA", "#70AD47", "#C0504D", "#778899", "#B5B5B5", "#B8860B", "#008B8B", "#8B4513", "#A0522D", "#FFA07A", "#DC143C", "#00BFFF", "#B8860B", "#483D8B", "#2F4F4F", "#00CED1", "#9400D3", "#008000", "#BC8F8F", "#006400", "#FF6347", "#DAA520", "#FF69B4", "#4682B4", "#32CD32", "#F4A460", "#87CEEB", "#FF8C00", "#FA8072", "#8B008B", "#2F4F4F", "#00FF00", "#008080", "#B22222", "#000080", "#FF4500", "#DAA520", "#A52A2A", "#8A2BE2", "#6495ED", "#00FA9A", "#9932CC", "#FFD700", "#20B2AA", "#708090", "#CD5C5C", "#F0E68C", "#7B68EE", "#FF00FF", "#F08080", "#87CEFA", "#7FFF00", "#4682B4", "#00FF00", "#4169E1", "#5F9EA0", "#0000CD", "#FA8072", "#AFEEEE", "#20B2AA", "#778899", "#FF0000", "#00FFFF", "#7FFFD4", "#396AB1", "#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#8B3D47", "#196D96", "#A3A3A3", "#FF9309", "#8B9BCA", "#A3A3A3", "#67BF5C", "#B5B5B5", "#D8B83F", "#D86375", "#6ED641", "#A3A3A3", "#D89841", "#8CACC6", "#DB6D00", "#97DFFC", "#B5B5B5", "#AB9F45", "#94BFA2", "#E17C05", "#71588F", "#DE6A61", "#6D6E71", "#D8B83F", "#9278A8", "#969696", "#000000", "#B5B5B5", "#D88E2D", "#76A233", "#6B4C9A", "#6D9CCF", "#8CACC6", "#8B9BCA", "#70AD47", "#C0504D", "#778899", "#B5B5B5", "#B8860B", "#008B8B", "#8B4513", "#A0522D", "#FFA07A", "#DC143C", "#00BFFF", "#B8860B", "#483D8B", "#2F4F4F", "#00CED1", "#9400D3", "#008000", "#BC8F8F", "#006400", "#FF6347", "#DAA520", "#FF69B4", "#4682B4", "#32CD32", "#F4A460", "#87CEEB", "#FF8C00", "#FA8072", "#8B008B", "#2F4F4F", "#00FF00", "#008080", "#B22222", "#000080", "#FF4500", "#DAA520", "#A52A2A", "#8A2BE2", "#6495ED", "#00FA9A", "#9932CC", "#FFD700", "#20B2AA", "#708090", "#CD5C5C", "#F0E68C", "#7B68EE", "#FF00FF", "#F08080", "#87CEFA", "#7FFF00", "#4682B4", "#00FF00", "#4169E1", "#5F9EA0", "#0000CD", "#FA8072", "#AFEEEE", "#20B2AA", "#778899", "#FF0000", "#00FFFF", "#7FFFD4"
    )

  }



  plot_aic_pvalue <- function(data,colors) {


    # 手动生成200种不同的颜色

    #usethis::use_package(package="tidyr",type="Import")
    #usethis::use_package(package="FactoMineR",type="Import")
    #usethis::use_package(package="leaps",type="Import")
    #usethis::use_package(package="lme4",type="Import")
    #usethis::use_package(package="stats",type="Import")
    #usethis::use_package(package="ggrepel",type="Import")




    # library(ggplot2)
    # library(ggrepel)

    # 根据 formula 列中是否包含 cov1 和 cov2 字符串设置形状
    data$shape <- ifelse(grepl("Cov1", data$formula) & grepl("Cov2", data$formula), "Covariates1+Covariates2 adjusted",
                         ifelse(grepl("Cov1", data$formula), "Covariates1 adjusted",
                                ifelse(grepl("Cov2", data$formula), "Covariates2 adjusted", "Non-adjusted")))

    data$formula <- paste(data$model ,":", data$formula)

    minAIC <- min(data$AIC)


    p <- ggplot2::ggplot(data, ggplot2::aes(x = Pvalue, y = AIC, group = formula, color = formula, size = 1, shape = shape)) +
      ggplot2::geom_point() +  # 绘制散点图
      ggplot2::scale_size_continuous(range = c(1, 5)) +  # 设置大小范围
      ggplot2::labs(x = "Pvalue", y = "AIC") +  # 设置坐标轴标签
      ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +  # 添加 x=0.05 的线
      ggplot2::annotate("text", x = 0.05, y = max(data$AIC), label = "Pvalue =0.05", vjust = -1, hjust = 0.5, color = "red") +  # 在图上标注 x=0.05
      ggplot2::scale_color_manual(values = colors) +  # 设置颜色
      ggplot2::scale_shape_manual(values = c("square", "triangle", "star", "circle")) +  # 设置形状
      # guides(color = FALSE, size = FALSE, shape = FALSE) +  # 隐藏图例
      ggplot2::theme_minimal()+  # 设置主题
      ggrepel::geom_text_repel(data=subset(data, Pvalue < 0.05),
                               ggplot2::aes(label=model),col="red",alpha = 0.8)+
      ggrepel::geom_text_repel(data=subset(data, AIC <= minAIC),
                               ggplot2::aes(label=model),col="black",alpha = 0.8)+
      ggplot2::theme(panel.grid = ggplot2::element_blank(),legend.position = "none")+  # 去掉网格线+隐藏图例（这一步不是必要的，但是有时候会因为绘图的图例太大而报错，所以可以隐藏图例）
      ggplot2::ggtitle(paste(data$celltype,"cell AIC-pvalue sactter plot"))

    print(p)
  }



  plist <- list()

  for(i in 1:length(mdlist)){

    mydata <- mdlist[[i]][[2]]
    plist[[i]] <-plot_aic_pvalue(mydata,colors=default_color)
    print(i)
  }

  return(plist)

}

