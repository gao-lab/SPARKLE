#' Perform Classic Comparison of Cell Proportions and Phenotypes
#'
#' This function calculates classic comparisons of cell proportions and phenotypes based on the provided CWAS data.
#'
#' @param cwas.data The CWAS data containing information on cell types, phenotypes, and rates.
#'
#' @return A list containing the results of statistical tests and a list of plots for each cell type.
#'
#' @export
#'
#' @examples cwas_classic_comparision(cwas.test.data)

cwas_classic_comparision <- function(cwas.data){

  #usethis::use_package(package="ggpubr",type="Import")

  #usethis::use_package(package="ggplot2",type="Import")
  #usethis::use_package(package="tidyverse",type="Import")
  #usethis::use_package(package="tidyverse",type="Import")
  #usethis::use_package(package="gtools",type="Import")



  ## 自定义主题
     #自定义主题2；
  #隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
  mytheme2 <- ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", colour = "gray30", size = 12),
      axis.line = ggplot2::element_line(size = 0.6, colour = "gray30"),
      axis.ticks = ggplot2::element_line(size = 0.6, colour = "gray30"),
      axis.ticks.length = ggplot2::unit(1.5, units = "mm"),
      plot.margin = ggplot2::unit(c(0.2, 0.2, 0.2, 0.2), units = "inches")
    )
  myattr <- attributes(cwas.data)

  mydata <- cwas.data
  mydata$Phenotype <- ifelse(mydata$Phenotype=="Disease",myattr$Disease_label,myattr$Control_label)
  cellnameall <- unique(mydata$Celltype)

  plotlist <- list()
  result.all <- c()
  for(i in 1:length(cellnameall)){

    cellname <- cellnameall[i]
    df <- mydata
    df <- subset(df, Celltype == cellname)
    plotlist[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = Phenotype, y = rate)) +
      ggplot2::stat_summary(mapping = ggplot2::aes(fill = Phenotype), fun = mean, geom = "bar", fun.args = list(mult = 1), width = 0.7) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
      ggplot2::geom_jitter(ggplot2::aes(), position =ggplot2::position_jitter(0.2), shape = 19, size = 2, alpha = 0.9) +
      mytheme2 +
      ggpubr::stat_compare_means(comparisons = list(c(myattr$Control_label,myattr$Disease_label))) +
      ggplot2::ggtitle(cellname)

    plotlist[[i]]

    group1 <- df[df$Phenotype == myattr$Control_label, "rate"]
    group2 <- df[df$Phenotype == myattr$Disease_label, "rate"]

    # 执行 Wilcoxon 符号秩检验（或 Mann-Whitney U 检验）
    test_result <- wilcox.test(group1, group2)
    test_result$p.value

    result.df <- data.frame(cellname,test_result$p.value)
    result.all <- rbind(result.all,result.df)

  }

  result.all$p.value_Non_adjusted=result.all$test_result.p.value
  result.all$adjusted.p.vaule_Bonferroni <-  p.adjust(result.all$test_result.p.value, method = "bonferroni")
  result.all$adjusted.p.vaule_Benjamini_Hochberg <-  p.adjust(result.all$test_result.p.value, method = "BH")
  result.all$adjusted.p.vaule_Holm <-  p.adjust(result.all$test_result.p.value, method = "holm")
  result.all$adjusted.p.vaule_Benjamini_Yekutieli <-  p.adjust(result.all$test_result.p.value, method = "BY")
  result.all$adjusted.p.vaule_Hommel <-  p.adjust(result.all$test_result.p.value, method = "hommel")
  result.all$adjusted.p.vaule_Hochberg <-  p.adjust(result.all$test_result.p.value, method = "hochberg")

  result.all$test_result.p.value=NULL


  return(list(result.all,plotlist))

}

