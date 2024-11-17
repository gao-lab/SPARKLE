#' Generate a Mediation Forest Plot
#'
#' This function generates a mediation forest plot to visualize mediation analysis results.
#'
#' @param data A data frame containing the necessary information for plotting the mediation forest plot.
#'
#' @return A mediation forest plot.
#'
#' @export
#'


Mediation_forest_plot <- function(data){


  tm <- forestploter::forest_theme(
    base_size = 10,        # 设置文本的基础大小
    # 设置可信区间的外观
    ci_pch = 15,           # 可信区间点的形状
    ci_col = "blue4",    # 可信区间的边框颜色
    ci_fill = "blue4",      # 可信区间的填充颜色
    ci_alpha = 0.8,        # 可信区间的透明度
    ci_lty = 1,            # 可信区间的线型
    ci_lwd = 1.5,          # 可信区间的线宽
    ci_Theight = 0.2,      # 设置T字在可信区间末端的高度，默认是NULL

    # 设置参考线的外观
    refline_lwd = 2,         # 参考线的线宽
    refline_lty = "dashed",  # 参考线的线型
    refline_col = "red4",  # 参考线的颜色

    # 设置垂直线的外观
    vertline_lwd = 1,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
    vertline_lty = "dashed",  # 垂直线的线型
    vertline_col = "grey20",  # 垂直线的颜色

    # 设置脚注的字体大小、字体样式和颜色
    footnote_cex = 0.6,            # 脚注字体大小
    footnote_fontface = "italic",  # 脚注字体样式
    footnote_col = "red4"          # 脚注文本的颜色
  )


  data$CI_lower <- c()

  data$CI_upper <- c()

  for(i in 1:length(data$CI)){
    matches <- regmatches(data$CI[i], gregexpr("-?\\d+\\.\\d+", data$CI[i]))
    extracted_numbers <- as.numeric(unlist(matches))
    data$CI_lower[i] <- extracted_numbers[1]
    data$CI_upper[i] <- extracted_numbers[2]

  }


  data$Pvalue <-  data$Indirect
  data$FDR <- p.adjust(data$Pvalue, method = "BH")
  data$Mediation <-gsub("^.*Mediation:","",data$Model)
  data$Celltype <-gsub("^.*X:","",data$Model)
  data$Celltype <-gsub("Mediation.*","",  data$Celltype)

  numeric_cols <- sapply(data, is.numeric)
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) round(x, 2))
  data$OR<- c(data$CI_lower+data$CI_upper)
  data$OR<- data$OR/2

  data$sig <- ifelse(data$Pvalue<0.001,"***",
                       ifelse(data$Pvalue<0.01,"**",
                              ifelse(data$Pvalue<0.05,"*","")))

  data$Pvalue <- ifelse(data$Pvalue<0.01,"<0.01",data$Pvalue )
  data$"P-value" <- paste0(data$Pvalue,data$sig)

  data1 <- data
  numeric_cols <- sapply(data1, is.numeric)
  data1[numeric_cols]  <- lapply(data[numeric_cols], function(x) format(x, scientific = TRUE, digits = 2))

  data$"OR(95%CI)" <- paste0(  data1$OR," ","[",data1$CI_lower,",",data1$CI_upper,"]")
  data$OR <- ifelse(data$OR > 20, 10, data$OR)
  data$"OR(95%CI)"[grepl("Inf", data1$OR)] <- "Not reliable"

  show_variable2 <- c("Celltype","Mediation","Effect","                        ","P-value")

  data$"                        " <- paste(rep("    ", length(data[,1])), collapse = " ")
  pp <- forestploter::forest(data[,show_variable2],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
                             est = data$Effect,          # 效应值，也就是HR列
                             lower = data$CI_lower,  # 置信区间下限
                             upper = data$CI_upper,  # 置信区间上限
                             sizes = 0.8,        # 黑框框的大小
                             ci_column = 4,             # 在第3列（可信区间列）绘制森林图
                             ref_line = 0,              # 添加参考线
                             arrow_lab = c("Down regulate", "Up regulate"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
                             #xlim = c(-1, 1),          # 设置x轴范围
                             #ticks_at = c( -1,0,1 ),  # 在指定位置添加刻度
                             theme = tm,                # 添加自定义主题
                             footnote = "Red line:Adjusted.Pvaule<0.05")  # 添加脚注信息

  pp <- forestploter::edit_plot(pp, row = which(data$Pvalue<0.0500001), gp = grid::gpar(col = "red4", fontface = "italic"))
  pp <- forestploter::insert_text(pp,
                                  text ="CWAS Mediation Analysis",
                                  col = 1:5,
                                  part = "header",
                                  gp =  grid::gpar(fontface = "bold"))
  pp <- forestploter::add_border(pp, part = "header", where = "bottom")
  return(pp)
}


