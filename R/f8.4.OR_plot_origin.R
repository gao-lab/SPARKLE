OR_plot_origin <- function(data,show_variable=c("Celltype","Model","OR(95%CI)","P-value","P-value adjusted"),plot.position=4,Tablename="Tablename"){
  
  cwas.data <- data[["All_models"]][[1]][["Rawdata"]]
  atri_info <- attributes(cwas.data)
  data <- data[[2]]
  show_variable2 <- show_variable[1:plot.position-1]
  show_variable2[plot.position] <- c(" ")
  num <- length(show_variable)
  num1 <- plot.position+1
  num2 <- length(show_variable)+1
  show_variable2[num1:num2] <- show_variable[plot.position:num]
  
  # 导入所需的包
  # library(metafor)
  #library(grid)
  #library(forestploter)
  
  # 创建数据框，包括研究名称、OR值和95%置信区间的上下限，以及事件数量和总样本量
  
  numeric_cols <- sapply(data, is.numeric)
  
  data1 <- data
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) round(x, 2))
  data1[numeric_cols]  <- lapply(data[numeric_cols], function(x) format(x, scientific = TRUE, digits = 2))
  
  data$" " <- paste(rep(" ", length(data[,1])), collapse = " ")
  data$"  " <- paste(rep(" ", length(data[,1])), collapse = " ")
  
  
  #data$CI_upper_table <- ifelse(data$CI_upper > 10, ">10", data$CI_upper)
  
  data$"OR(95%CI)" <- paste0(data1$OR," ","[",data1$CI_lower,",",data1$CI_upper,"]")
  data$OR <- ifelse(data$OR > 20, 10, data$OR)
  data$"OR(95%CI)"[grepl("Inf", data1$OR)] <- "Not reliable"
  # 对所有class是numeric的列保留小数点后两位数字
  data$CI_upper <- ifelse(data$CI_upper > 10, 10, data$CI_upper)
  
  data$"P-value" <-ifelse(data$Pvalue < 0.01, "<0.01", data$Pvalue)
  data$"P-value" <- paste(data$"P-value",data$significance)
  data$"P-value" <- paste0("  ",data$"P-value")
  data$"P-value" [grepl("Inf", data1$OR)] <- "Not reliable"
  
  tryCatch({
    data$"P-value adjusted" <-ifelse(data$adjustedPvalue < 0.01, "<0.01", data$adjustedPvalue)
    data$"P-value adjusted" <- paste(data$"P-value adjusted",data$significance.adjusted)
    data$"P-value adjusted" <- paste0("  ",data$"P-value adjusted")
    data$"P-value adjusted"[grepl("Inf", data1$OR)] <- "Not reliable"
  }, error = function(e) {
    # 捕获到错误时执行的代码
    # 打印出错信息
    print(paste("Warning Message:","No Adjusted p value in this data"))
    
    show_variable2<<-  show_variable2[show_variable2 != "P-value adjusted"]
    
  })
  
  data$Pvalue[grepl("Inf", data1$OR)] <- "Not reliable"
  
  if(is.null(data$formula)){ data$formula <- c("Integrated model") }
  
  data$Model <- data$formula
  
  data$Celltype <- data$celltype
  data$Significance <- data$significance
  
  
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
  
  # 绘制森林图
  
  
  footnote1 <- c("* p<0.05 \ ** p<0.01 \ *** p<0.001 ")
  
  info_parts <- list(
    paste("Phenotype:", ifelse(is.null(atri_info$Phenotype), "Not Provided", atri_info$Phenotype),"  ",
          "Control_label:", ifelse(is.null(atri_info$Control_label), "Not Provided", atri_info$Control_label),"  ",
          "Disease_label:", ifelse(is.null(atri_info$Disease_label), "Not Provided", atri_info$Disease_label)),
    paste("Group:", ifelse(is.null(atri_info$Group), "Not Provided", atri_info$Group),"  ",
          "Subgroup:", ifelse(is.null(atri_info$Subgroup), "Not Provided", atri_info$Subgroup),"  ",
          "Covariate 1:", ifelse(is.null(atri_info$Covariate1), "Not Provided", atri_info$Covariate1),"  ",
          "Covariate 2:", ifelse(is.null(atri_info$Covariate2), "Not Provided", atri_info$Covariate2) )
  )
  
  
  info <- paste(info_parts, collapse = "\n")
  footnote2 <- paste0(footnote1, info)
  
  
  
  pp <- forestploter::forest(data[,show_variable2],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
                             est = data$OR,          # 效应值，也就是HR列
                             lower = data$CI_lower,  # 置信区间下限
                             upper = data$CI_upper,  # 置信区间上限
                             sizes = 0.8,        # 黑框框的大小
                             ci_column = 4,             # 在第3列（可信区间列）绘制森林图
                             ref_line = 1,              # 添加参考线
                             arrow_lab = c("Low risk", "High Risk"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
                             xlim = c(-1, 10),          # 设置x轴范围
                             ticks_at = c(0,1, 3, 5),  # 在指定位置添加刻度
                             theme = tm,                # 添加自定义主题
                             footnote = footnote2)  # 添加脚注信息
  
  pp <- forestploter::edit_plot(pp, row = which(data$Pvalue<0.0500001), gp = grid::gpar(col = "red4", fontface = "italic"))
  
  pp <- forestploter::insert_text(pp,
                                  text =Tablename,
                                  col = 1:5,
                                  part = "header",
                                  gp =  grid::gpar(fontface = "bold"))
  pp <- forestploter::add_border(pp, part = "header", where = "bottom")
  print(pp)
  
  return(pp)
  
}
