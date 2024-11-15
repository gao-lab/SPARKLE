 
library(ggplot2)

# 计算Jaccard Index的函数
calculate_jaccard_index <- function(setA, setB) {
  
  
  
  intersection_size <- length(intersect(setA, setB))
  union_size <- length(union(setA, setB))
  jaccard_index <- intersection_size / union_size
  return(jaccard_index)
}

# 可视化Jaccard Index的函数
visualize_jaccard_index <- function(DEG.MTF,EMTF.mediation,com.num=NULL) {
  
  
  
  cwas.data <- separate(EMTF.mediation[["Summary"]], Model, into = c("X", "Mediation"), sep = " Mediation:")
  cwas.data$Proportion <- cwas.data$Indirect/cwas.data$Total
  cwas.data_sorted <- cwas.data %>% arrange(Indirect)
  
  setA <- cwas.data_sorted$Mediation
  
  
  if(!is.null(DEG.MTF$gene[1]) ){
    DEG.MTF <-  DEG.MTF%>% arrange(p_val_adj)
    DEGgenes <-DEG.MTF$gene }else{
      DEG.MTF$gene<-rownames(DEG.MTF)
      DEG.MTF <-  DEG.MTF%>% arrange(p_val_adj)
      DEGgenes <-DEG.MTF$gene
    }
  
  setB <- DEGgenes
  
  if(is.null(com.num)){
    com.num <- min(length(setB),length(setA))
  } 
  setA <- setA[1:com.num]
  setB <- setB[1:com.num]
  
   
  jaccard_indices <- numeric(com.num)
  
  for (i in 1:com.num) {
    prefixA <- setA[1:i]
    prefixB <- setB[1:i]
    jaccard_indices[i] <- calculate_jaccard_index(prefixA, prefixB)
  }
  
  data <- data.frame(
    Index = 1:com.num,
    JaccardIndex = jaccard_indices
  )
  
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
  
  P <- ggplot(data, aes(x = Index, y = JaccardIndex)) +
    geom_line(color = "gray") +
    geom_point(color = "black") +
    ggtitle("Jaccard Index Over Prefixes") +
    xlab("Prefix Length") +
    ylab("Jaccard Index") +
    ylim(c(0,1))+
    mytheme2+ggsci::scale_color_aaas()+ggsci::scale_fill_aaas()
  
  print(P)
  
  return(data)
}
 

