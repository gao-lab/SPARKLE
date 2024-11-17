Simpson_plot <- function(fibro.cwas.select,cellname=cellname){

  fibro.cwas.select$Phenotype2 <- ifelse(fibro.cwas.select$Phenotype=="Disease",1,0)
  fibro.cwas.select$Subgroup <- as.factor(fibro.cwas.select$Subgroup )
  library(ggplot2)
  top.mar = 0.2
  right.mar = 0.2
  botton.mar = 0.2
  left.mar = 0.2
  ## 自定义主题
  ## 合并上面的参数，将其合并成mythemel
  mytheme1 <- theme(panel.background = element_blank(),
                    axis.ticks.length=unit(1.6,"mm"),
                    plot.margin=unit(x=c(top.mar,right.mar,botton.mar,left.mar),
                                     units="inches"))
  #自定义主题2；
  #隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
  mytheme2<-theme_classic()+
    theme(text=element_text(family = "sans",colour ="gray30",size = 12),
          axis.line = element_line(size = 0.6,colour = "gray30"),
          axis.ticks = element_line(size = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"),
          plot.margin=unit(x=c(top.mar,right.mar,botton.mar,left.mar),
                           units="inches"))

  label_mapping <- function(x) {
    ifelse(x == 0, "Control", ifelse(x == 1, "Disease", ""))
  }

  p1 <- ggplot(fibro.cwas.select, aes(x =rate , y = Phenotype2, color = Subgroup)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_y_continuous(labels = label_mapping) +  # 使用自定义函数进行标签映射
    labs(title = "Group-wise",
         x = paste0(cellname," Cell Proportion"),
         y = "Phenotype")+mytheme2+ggsci::scale_color_aaas()+ggsci::scale_color_aaas()

  # 合并数据的关系




  p2 <-ggplot(fibro.cwas.select, aes(x = rate, y = Phenotype2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(title = "Overall",
         x = "Cell Proportion",
         y = "Phenotype") +
    scale_y_continuous(labels = label_mapping) +  # 使用自定义函数进行标签映射
    theme(axis.text.y = element_text(size = 12),  # 设置纵坐标文本大小
          axis.ticks.y = element_line(size = 0.5))+mytheme2
  library(patchwork)
  print(p2+p1)

}


Simpson_cal <- function(fibro.cwas.select,subgroup) {

  # 创建二元表型变量
  fibro.cwas.select$Phenotype2 <- ifelse(fibro.cwas.select$Phenotype == "Disease", 1, 0)

  fibro.cwas.select1 <-  fibro.cwas.select%>% dplyr::select_all() %>% dplyr::filter(Subgroup==subgroup[1])
  model1 <- lm(Phenotype2 ~ rate, data = fibro.cwas.select1)
  beta1 <- coef(model1)["rate"]

  fibro.cwas.select2 <-  fibro.cwas.select%>% dplyr::select_all() %>% dplyr::filter(Subgroup==subgroup[2])
  model2 <- lm(Phenotype2 ~ rate, data = fibro.cwas.select2)
  beta2 <- coef(model2)["rate"]


  # 计算整体数据的回归系数
  overall_model <- lm(Phenotype2 ~ rate, data = fibro.cwas.select)
  overall_beta <- coef(overall_model)["rate"]

  df <- data.frame(group1=subgroup[1],group2=subgroup[2],model1=beta1,model2=beta2,allmodel=overall_beta)

  # 返回结果
  return(df)
}


sf.rds <- function(variable,vn = "viarablename"){
  dir.create(vn)
  namefile <- gsub(" ","",paste0(vn,Sys.time(),".rds"))
  namefile <- gsub(":","",namefile)
  #print(variable)
  saveRDS(variable,file = paste0(vn,"/",namefile),compress = F)
  output <- paste0("文件 ",namefile,"存储在",vn,"文件夹")
  return(output)
}
