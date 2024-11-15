
diff_gene <- function(fibro.cwas.data){
  fibro.cwas.data.mtf <- fibro.cwas.data %>% dplyr::select_all() %>% dplyr::filter(Celltype=="MTF")
  difsample <- setdiff(fibro.cwas.data$Sample,fibro.cwas.data.mtf$Sample)
  fibro.cwas.data.mtf
  fibro.cwas.data.mtf <- fibro.cwas.data.mtf[,-c(1,2,3,5:6)]
  # 加载必要的库
  library(dplyr)
  library(broom)

  # 假设你的数据已经加载到一个名为fibro_cwas_data的数据框中
  # fibro_cwas_data <- read.csv('path_to_your_data.csv') # 读取数据的示例

  # 创建一个空的数据框来存储结果
  results <- data.frame(Gene = character(), Estimate = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

  names(fibro.cwas.data.mtf) <- gsub("-","_",names(fibro.cwas.data.mtf))

  # 获取所有的列名，假设除了'rate'列以外，其他都是基因列
  genes <- setdiff(names(fibro.cwas.data.mtf), c('rate'))
  fibro.cwas.data.mtf$Phenotype <- ifelse(fibro.cwas.data.mtf$Phenotype=="Disease",1,0)
  # 对每个基因进行线性回归分析
  for (i in 2:length(genes)) {
    gene <- genes[i]
    formula <- as.formula(paste("Phenotype ~", gene))
    model <- glm(formula, data = fibro.cwas.data.mtf,family = "binomial")
    summary_model <- summary(model)

    # 提取回归系数和p值
    estimate <- summary_model$coefficients[2, "Estimate"]
    p_value <- summary_model$coefficients[2, "Pr(>|z|)"]

    # 将结果添加到数据框中
    results <- rbind(results, data.frame(Gene = gene, Estimate = estimate, P.Value = p_value, stringsAsFactors = FALSE))
  }

  # 按照p值对结果进行排序
  results <- results %>% arrange(P.Value)
  results$p_val_adj=results$P.Value
  results$gene <- results$Gene

  return(results)

}

cellrate_gene <- function(fibro.cwas.data){

  fibro.cwas.data.mtf <- fibro.cwas.data %>% dplyr::select_all() %>% dplyr::filter(Celltype=="MTF")
  difsample <- setdiff(fibro.cwas.data$Sample,fibro.cwas.data.mtf$Sample)
  fibro.cwas.data.mtf
  fibro.cwas.data.mtf <- fibro.cwas.data.mtf[,-c(1,2,4:6)]
  # 加载必要的库
  library(dplyr)
  library(broom)

  # 假设你的数据已经加载到一个名为fibro_cwas_data的数据框中
  # fibro_cwas_data <- read.csv('path_to_your_data.csv') # 读取数据的示例

  # 创建一个空的数据框来存储结果
  results <- data.frame(Gene = character(), Estimate = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

  names(fibro.cwas.data.mtf) <- gsub("-","_",names(fibro.cwas.data.mtf))

  # 获取所有的列名，假设除了'rate'列以外，其他都是基因列
  genes <- setdiff(names(fibro.cwas.data.mtf), c('rate'))

  # 对每个基因进行线性回归分析
  for (i in 1:length(genes)) {
    gene <- genes[i]
    formula <- as.formula(paste("rate ~", gene))
    model <- lm(formula, data = fibro.cwas.data.mtf)
    summary_model <- summary(model)

    # 提取回归系数和p值
    estimate <- summary_model$coefficients[2, "Estimate"]
    p_value <- summary_model$coefficients[2, "Pr(>|t|)"]

    # 将结果添加到数据框中
    results <- rbind(results, data.frame(Gene = gene, Estimate = estimate, P.Value = p_value, stringsAsFactors = FALSE))
  }

  # 按照p值对结果进行排序
  results <- results %>% arrange(P.Value)
  results$p_val_adj=results$P.Value
  results$gene <- results$Gene

  return(results)
}

target_test_4group  <- function(rawadata,DEG.MTF,EMTF.mediation,Drug.target=NULL,com.num=NULL,NamesofDEG="DEG",figname="Target Potential Analysis of Top genes"){

  results.diff <- diff_gene(fibro.cwas.data)
  results.cell <- cellrate_gene(fibro.cwas.data)

  library(tidyr)
  library(dplyr)
  if(is.null(Drug.target)){
    fpath <- system.file("extdata", "TTDdatabase240523.rdata", package="SPARKLE")
    load(fpath)
    Drug.target <- unique(TTDdatabase[["Target_info"]][["Target_gene"]])
  }

  #DEG.MTF <-  DEG.MTF%>% arrange(-abs(avg_log2FC))


  if(!is.null(DEG.MTF$gene[1]) ){
    DEG.MTF <-  DEG.MTF%>% arrange(p_val_adj)
    DEGgenes <-DEG.MTF$gene }else{
      DEG.MTF$gene<-rownames(DEG.MTF)
      DEG.MTF <-  DEG.MTF%>% arrange(p_val_adj)
      DEGgenes <-DEG.MTF$gene
    }



  cwas.data <- separate(EMTF.mediation[["Summary"]], Model, into = c("X", "Mediation"), sep = " Mediation:")
  cwas.data$Proportion <- cwas.data$Indirect/cwas.data$Total

  cwas.data_sorted <- cwas.data %>% arrange(Indirect)
  cwas.data_sorted3 <- cwas.data %>% arrange(Proportion)
  cwas.data_sorted4 <- cwas.data %>% arrange(-abs(Effect))

  cwas.data_sorted2 <- cwas.data_sorted %>% dplyr::select_all() %>% dplyr::filter(Effect>0)


  DrugGene <- Drug.target
  # 定义函数
  compute_overlap_statistics <- function(Indirect, DrugGene,total_genes=65859) {
    # 初始化数据框
    result <- data.frame(
      K = integer(),
      OverlapCount = integer(),
      OverlapRatio = numeric(),
      OR = numeric()
    )

    # 生成数据表格
    for (k in 1:length(Indirect)) {
      # 取前K个字符
      substring_indirect <- Indirect[1:k]

      # 计算重叠数量
      overlap_count <- sum(substring_indirect %in% DrugGene)

      # 计算重叠比例
      overlap_ratio <- overlap_count / k

      # 计算 OR 值
      a <- overlap_count
      b <- k - overlap_count
      c <- length(DrugGene ) - overlap_count
      d <-  total_genes- a - b - c

      # 确保没有零值以避免分母为零的问题

      if (a == 0 | b == 0 | c == 0 | d == 0) {
        or_value <- 1  # 这里你可以选择如何处理零值
        ci_lower <- 1
        ci_upper <- 1
      }else {
        or_value <- (a * d) / (b * c)
        # 计算标准误差
        se <- sqrt(1/a + 1/b + 1/c + 1/d)
        # 计算log OR的95%置信区间
        log_or <- log(or_value)
        ci_lower <- exp(log_or - 1.96 * se)
        ci_upper <- exp(log_or + 1.96 * se)
      }

      # 添加到结果数据框
      result <- rbind(result, data.frame(
        K = k,
        OverlapCount = overlap_count,
        OverlapRatio = overlap_ratio,
        OR = or_value,
        CI_lower = ci_lower,
        CI_upper = ci_upper
      ))
    }

    # 返回结果数据框
    return(result)
  }
  # 示例字符向量
  # result0 <- compute_overlap_statistics(Marker.MTF$gene, Drug.target)
  # # result1$DEG <- result1$OverlapRatio
  # result0$Group <- "MarkerGene"
  # 调用函数

  result5 <- compute_overlap_statistics(Indirect=cwas.data_sorted2$Mediation, DrugGene=Drug.target)
  # result2$PMG <- result2$OverlapRatio
  result5$Group <- "PMG.up"


  result1 <- compute_overlap_statistics(DEGgenes, Drug.target)
  # result1$DEG <- result1$OverlapRatio
  result1$Group <- NamesofDEG
  result2 <- compute_overlap_statistics(cwas.data_sorted$Mediation, Drug.target)
  # result2$PMG <- result2$OverlapRatio
  result2$Group <- "PMG"

  result3 <- compute_overlap_statistics(cwas.data_sorted3$Mediation, Drug.target)
  # result2$PMG <- result2$OverlapRatio
  result3$Group <- "PMG.Prop"

  result4 <- compute_overlap_statistics(cwas.data_sorted4$Mediation, Drug.target)
  # result2$PMG <- result2$OverlapRatio
  result4$Group <- "PMG.Effect"

  results.diff <- compute_overlap_statistics(results.diff$gene, Drug.target)
  results.diff$Group <- "Pseudobulk.DEG"

  results.cell <- compute_overlap_statistics(results.cell$gene, Drug.target)
  results.cell$Group <- "Pseudobulk.CMG"

  if(is.null(com.num)){
    com.num <- length(intersect(result1$K,result2$K))
  }


  result <- rbind(result1[1:com.num,],result2[1:com.num,])
  result <- rbind(result,results.cell[1:com.num,])
  result <- rbind(result,results.diff[1:com.num,])
  #result <- rbind(result,result5[1:com.num,])
  #result <- rbind(result,result0[1:com.num,])
  #  result <- rbind(result,result3[1:com.num,])
  # result <- rbind(result,result4[1:com.num,])
  #
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

  library(ggsci)
  # 查看结果
  print(result)

  p <- ggplot(result, aes(x = K, y = OverlapRatio, color = Group)) +
    geom_line() +
    geom_point() +
    labs(title = figname,
         x = "Top target gene count",
         y = "Hit Ratio") +
    mytheme2+ggsci::scale_color_aaas()+ggsci::scale_fill_aaas()
  print(p)

  p1  <-  ggplot(result, aes(x = K, y = OR, color = Group)) +
    geom_line() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_point() +
    labs(title = figname,
         x = "Top target gene count",
         y = "Odds Ratio") +

    ggsci::scale_color_aaas() +
    ggsci::scale_fill_aaas() +
    mytheme2  # 你可以替换成你自己的主题 mytheme2
  print(p1)

  p2 <-  ggplot(result, aes(x = K, y = OR, color = Group)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    labs(title = figname,
         x = "Top target gene count",
         y = "Odds Ratio") +

    ggsci::scale_color_aaas() +
    ggsci::scale_fill_aaas() +
    mytheme2  # 你可以替换成你自己的主题 mytheme2

  print(p2)

  p3<- ggplot(result, aes(x = K, y = OR, color = Group, fill = Group)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Group), alpha = 0.2, color = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    labs(title = figname,
         x = "Top target gene",
         y = "Odds Ratio") +
    ggsci::scale_color_aaas() +
    ggsci::scale_fill_aaas() +

    mytheme2   # 你可以替换成你自己的主题 mytheme2

  print(p3)

  return(result)

}
