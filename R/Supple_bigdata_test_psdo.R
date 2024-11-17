# library(tidyverse)
# library(ggplot2)
# library(dplyr)
pvalue <- function(covid.data.selected){

  if (covid.data.selected$Phenotype[1] %in% c(0, 1)) {
  }else {
    covid.data.selected$Phenotype <- ifelse(covid.data.selected$Phenotype == "Control",
                                            0, 1)
  }

  allcelltype <- unique(covid.data.selected$Celltype)

  table(covid.data.selected$Celltype)

  dfall <- c()
  for(i in 1:length(allcelltype)){

    covid.data.selected2 <-  covid.data.selected%>% dplyr::select_all() %>% dplyr::filter(Celltype==allcelltype[i])
    #print(allcelltype[i])
    grouped1 <- covid.data.selected2 %>%  dplyr::select_all() %>% dplyr::filter(Phenotype=="0")
    grouped2 <- covid.data.selected2 %>%  dplyr::select_all() %>% dplyr::filter(Phenotype=="1")
    # 提取两个分组的rate
    group1 <- grouped1$rate
    group2 <- grouped2$rate
    if(length(group1)>0&length(group2)>0){
    # 进行Wilcoxon秩和检验
    wilcox_test_result <- wilcox.test(group1, group2, paired = F)
    fixed_model1 <- stats::glm(Phenotype ~ rate, family = "binomial", data = covid.data.selected2)
    fixed_cov <- summary(fixed_model1)
    #fixed_model2 <- stats::glm(Phenotype ~ rate , family = "binomial", data = covid.data.selected)
    #random_model <- lme4::glmer(Phenotype ~ rate + (1 | Subgroup), family = "binomial", data = covid.data.selected)

    pvalue <- fixed_cov$coefficients["rate",4]

    df <- data.frame(cellname=allcelltype[i],wilcox=wilcox_test_result$p.value,fixed=pvalue)
    dfall <- rbind(df,dfall)
    }

  }

  return(dfall)

}


alltest <- function(covid.cwas.data,num){

  dfall <- c()
  num.control <- round(num/2,0)
  covid.cwas.data.control <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="Control")
  alldataset.control <- unique(covid.cwas.data.control$Subgroup)
  selectedgroup.control <-as.data.frame( t(combn(alldataset.control,num.control)))
  random_sample.control <- sample_n(selectedgroup.control, 5,replace = TRUE)

  num.disease <- num-num.control
  covid.cwas.data.disease <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="Disease")
  alldataset.disease <- unique(covid.cwas.data.disease$Subgroup)
  selectedgroup.disease <-as.data.frame( t(combn(alldataset.disease,num.disease)))
  random_sample.disease <- sample_n(selectedgroup.disease, 5,replace = TRUE)

  random_sample <- cbind(random_sample.control,random_sample.disease)

  alldiff <- unique(covid.cwas.data$Subgroup)

  for(i in 1:length(random_sample[,1])){

    sub <- unique(as.character( random_sample[i,]))
    difnum <- num-length(sub)
    difsetsub <- setdiff(alldiff,sub)
    covid.data.selected <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Subgroup%in%c(sub,difsetsub[1:difnum]))
    #mydata <- pvalue(covid.data.selected)

    df <- tryCatch({
      pvalue(covid.data.selected)
    }, error = function(e) {
      message("Error and Skipping this celltype.")
      return(NULL) # 返回 NULL 以便后续处理
    })

    # 如果 df 不为 NULL，则将其绑定到 dfall
    if (!is.null(df)) {
      #print(df)
      dfall <- rbind(dfall, df)


    }



  }


  return(dfall)
}



alltest_acc <- function(covid.cwas.data,num){


  generate_random_combinations <- function(alldataset, num, n_samples) {
    n <- length(alldataset)
    combinations <- matrix(NA, nrow = n_samples, ncol = num)
    for (i in 1:n_samples) {
      combinations[i, ] <- sample(alldataset, num)
    }
    as.data.frame(combinations)
  }

  dfall <- c()
  num.control <- round(num/2,0)
  covid.cwas.data.control <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="Control")
  alldataset.control <- unique(covid.cwas.data.control$Subgroup)
  # selectedgroup.control <-as.data.frame( t(combn(alldataset.control,num.control)))
  # random_sample.control <- sample_n(selectedgroup.control, 5,replace = TRUE)
  random_sample.control <-generate_random_combinations(alldataset.control,num.control,5)

  num.disease <- num-num.control
  covid.cwas.data.disease <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="Disease")
  alldataset.disease <- unique(covid.cwas.data.disease$Subgroup)
  # selectedgroup.disease <-as.data.frame( t(combn(alldataset.disease,num.disease)))
  # random_sample.disease <- sample_n(selectedgroup.disease, 5,replace = TRUE)


  random_sample.disease <- tryCatch({
    generate_random_combinations(alldataset.disease,num.disease,5)
  }, error = function(e) {
    return(NULL) # 返回 NULL 以便后续处理
  })

  if(is.null(random_sample.disease)){
    random_sample.disease <- as.data.frame(t(replicate(5, alldataset.disease)))
  }


  random_sample <- cbind(random_sample.control,random_sample.disease)

  alldiff <- unique(covid.cwas.data$Subgroup)

  for(i in 1:length(random_sample[,1])){

    sub <- unique(as.character( random_sample[i,]))
    difnum <- num-length(sub)
    difsetsub <- setdiff(alldiff,sub)
    covid.data.selected <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Subgroup%in%c(sub,difsetsub[1:difnum]))
    #mydata <- pvalue(covid.data.selected)

    df <- tryCatch({
      pvalue(covid.data.selected)
    }, error = function(e) {
      message("Error and Skipping this celltype.")
      return(NULL) # 返回 NULL 以便后续处理
    })

    # 如果 df 不为 NULL，则将其绑定到 dfall
    if (!is.null(df)) {
      #print(df)
      dfall <- rbind(dfall, df)


    }



  }


  return(dfall)
}


alltest1 <- function(covid.cwas.data){


  covid.cwas.data.control <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="Control")

  covid.cwas.data.disease <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Phenotype=="disease")

  random_sample  <- intersect(covid.cwas.data.control$Subgroup,covid.cwas.data.control$Subgroup)

  dfall <- c()


  for(i in 1:length(random_sample)){
    covid.data.selected <-  covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(Subgroup%in%c(random_sample[i]))
    #mydata <- pvalue(covid.data.selected)

    df <- tryCatch({
      pvalue(covid.data.selected)
    }, error = function(e) {
      message("Error and Skipping this celltype.")
      return(NULL) # 返回 NULL 以便后续处理
    })

    # 如果 df 不为 NULL，则将其绑定到 dfall
    if (!is.null(df)) {
      #print(df)
      dfall <- rbind(dfall, df)


    }



  }


  return(dfall)
}

alltest_cal <- function(covid.cwas.data){


  alldataset <- unique(covid.cwas.data$Subgroup)
  #covid.cwas.data2 <- covid.cwas.data%>% dplyr::select_all() %>% dplyr::filter(!Celltype=="Epi")
  #myd <- alltest(covid.cwas.data2,2)

  mydfall <- alltest1(covid.cwas.data)
  mydfall$groupnum=1

  n <- length(alldataset)
  print(n)
  for(i in 2:n){


    myd <- tryCatch({
      alltest_acc(covid.cwas.data,i)
    }, error = function(e) {
      message("Error and Skipping this celltype.")
      return(mydfall) # 返回 NULL 以便后续处理
    })

    if(!is.null(myd)){
      myd$groupnum <- i
      mydfall <- rbind(mydfall, myd)
      print(i)

      # p <- ggplot(mydfall, aes(x = groupnum, y = wilcox, color = cellname, group = cellname)) +
      #   geom_point() +  # 添加散点图层
      #   geom_line() +   # 添加折线图层
      #   labs(title = "Scatter and Line Plot", x = "Group Number", y = "Wilcox Value") +  # 添加标题和坐标轴标签
      #   theme_minimal()  # 使用简洁主题
      #
      # # 显示图形
      # print(p)

    }


  }

  return(mydfall)

}


plot_bigdata<- function(covid.cwas.data,mydfall,cutoff=0.5){

 mean_values_all <- c()
 allcelltype <- unique(mydfall$cellname)
  for(i in 1:length(allcelltype)){
    mydfall.selected <- mydfall%>% dplyr::select_all() %>% dplyr::filter(cellname==allcelltype[i])
    covid.cwas.data2 <- covid.cwas.data %>% dplyr::select_all() %>% dplyr::filter(Celltype==allcelltype[i])
    df <- pvalue(covid.cwas.data2)
    df$sig <- ifelse(df$wilcox<0.05,1,0)
    mydfall.selected$sig <- ifelse(mydfall.selected$wilcox<0.05,1,0)
    mydfall.selected$consist <- ifelse(mydfall.selected$sig== df$sig,1,0)

    mean_values <- mydfall.selected %>%
      group_by(groupnum) %>%
      summarise(mean_wilcox = mean(consist, na.rm = TRUE))
    maxgroup <- max(mean_values$groupnum)

    mean_values <- as.data.frame(mean_values)

    if(mean_values[which(mean_values$groupnum==maxgroup),"mean_wilcox"]<cutoff){
      mean_values$mean_wilcox <- -(mean_values$mean_wilcox-1)
    }
    mean_values$cellname <- allcelltype[i]
    mean_values_all <- rbind(mean_values,mean_values_all)

  }


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


  colors <- c(ggsci::pal_aaas("default")(10),"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#1a1a1a", "#ffcc00", "#a6cee3", "#b2df8a", "#fb9a99",
              "#fdbf6f", "#cab2d6", "#ffff99", "#1f78b4", "#33a02c",
              "#6a3d9a", "#b15928", "#1b9e77", "#ff8e00", "#33a02c",
              "#a6761d", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928",
              "#1b9e77", "#ff8e00", "#33a02c", "#a6761d", "#e31a1c",
              "#ff7f00", "#6a3d9a", "#b15928", "#1b9e77", "#ff8e00", "#2b22a4")


  p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_wilcox, color = cellname, linetype = cellname)) +
    geom_line(size = 1.2, alpha = 0.8, linetype = "dashed") +  # 添加平均值折线图层，设置为点线
    geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
    scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
    labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
    theme_minimal() +
    scale_color_manual(values = colors) +  # 设置颜色
    mytheme2
print(p)
p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_wilcox, color = cellname)) +
  geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
  scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
  labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
  theme_minimal() +
  scale_color_manual(values = colors) +  # 设置颜色
  geom_smooth(method = "lm", se = FALSE)+mytheme2  # 添加直线拟合层，使用线性模型拟合数据，不显示标准误差带
print(p)

p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_wilcox, color = cellname)) +
  geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
  scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
  labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
  theme_minimal() +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  geom_smooth(method = "loess", se = FALSE) +mytheme2 # 使用 LOESS 方法进行曲线拟合，不显示标准误差带
print(p)

library(ggplot2)
groupnum <- max(mean_values_all$groupnum)
# Function to fit linear model and ensure line passes through (12, 1)
fit_line_through_point <- function(df, x0 = groupnum, y0 = 1) {
  # 创建一个新的变量，用于中心化数据
  df$groupnum_centered <- df$groupnum - x0
  df$mean_wilcox_centered <- df$mean_wilcox - y0

  # 拟合中心化后的数据，不需要截距项
  lm_fit <- lm(mean_wilcox_centered ~ groupnum_centered + 0, data = df)

  # 斜率
  slope <- coef(lm_fit)[1]

  # 计算截距，使得直线通过 (12, 1)
  intercept <- y0 - slope * x0

  return(c(intercept, slope))
}

# Fit linear model for each cellname group
models <- by(mean_values_all, mean_values_all$cellname, fit_line_through_point)

# Extract intercepts and slopes from models
intercepts <- sapply(models, function(m) m[1])
slopes <- sapply(models, function(m) m[2])

# Create a data frame with intercepts, slopes, and cellnames
line_data <- data.frame(
  intercepts = intercepts,
  slopes = slopes,
  cellname = unique(mean_values_all$cellname)
)



p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_wilcox, color = cellname)) +
  geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
  scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
  labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
  geom_abline(data = line_data, aes(intercept = intercepts, slope = slopes, color = cellname), size = 1) +
  theme_minimal() +
  scale_color_manual(values = colors) +  # 设置颜色
  mytheme2 # 使用 LOESS 方法进行曲线拟合，不显示标准误差带
print(p)


return(mean_values_all)

}

plot_bigdata_SPARKLE<- function(covid.cwas.data,mydfall,cutoff=0.5){

  mean_values_all <- c()
  allcelltype <- unique(mydfall$cellname)
  for(i in 1:length(allcelltype)){
    mydfall.selected <- mydfall%>% dplyr::select_all() %>% dplyr::filter(cellname==allcelltype[i])
    covid.cwas.data2 <- covid.cwas.data %>% dplyr::select_all() %>% dplyr::filter(Celltype==allcelltype[i])
    df <- pvalue(covid.cwas.data2)
    df$sig <- ifelse(df$fixed<0.05,1,0)
    mydfall.selected$sig <- ifelse(mydfall.selected$fixed<0.05,1,0)
    mydfall.selected$consist <- ifelse(mydfall.selected$sig== df$sig,1,0)

    mean_values <- mydfall.selected %>%
      group_by(groupnum) %>%
      summarise(mean_fixed = mean(consist, na.rm = TRUE))
    maxgroup <- max(mean_values$groupnum)

    mean_values <- as.data.frame(mean_values)

    if(mean_values[which(mean_values$groupnum==maxgroup),"mean_fixed"]<cutoff){
      mean_values$mean_fixed <- -(mean_values$mean_fixed-1)
    }
    mean_values$cellname <- allcelltype[i]
    mean_values_all <- rbind(mean_values,mean_values_all)

  }


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


  colors <- c(ggsci::pal_aaas("default")(10),"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
              "#1a1a1a", "#ffcc00", "#a6cee3", "#b2df8a", "#fb9a99",
              "#fdbf6f", "#cab2d6", "#ffff99", "#1f78b4", "#33a02c",
              "#6a3d9a", "#b15928", "#1b9e77", "#ff8e00", "#33a02c",
              "#a6761d", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928",
              "#1b9e77", "#ff8e00", "#33a02c", "#a6761d", "#e31a1c",
              "#ff7f00", "#6a3d9a", "#b15928", "#1b9e77", "#ff8e00", "#2b22a4")


  p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_fixed, color = cellname, linetype = cellname)) +
    geom_line(size = 1.2, alpha = 0.8, linetype = "dashed") +  # 添加平均值折线图层，设置为点线
    geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
    scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
    labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
    theme_minimal() +
    scale_color_manual(values = colors) +  # 设置颜色
    mytheme2
  print(p)
  p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_fixed, color = cellname)) +
    geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
    scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
    labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
    theme_minimal() +
    scale_color_manual(values = colors) +  # 设置颜色
    geom_smooth(method = "lm", se = FALSE)+mytheme2  # 添加直线拟合层，使用线性模型拟合数据，不显示标准误差带
  print(p)

  p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_fixed, color = cellname)) +
    geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
    scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
    labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
    theme_minimal() +
    ggsci::scale_fill_aaas() +
    ggsci::scale_color_aaas() +
    geom_smooth(method = "loess", se = FALSE) +mytheme2 # 使用 LOESS 方法进行曲线拟合，不显示标准误差带
  print(p)

  library(ggplot2)
  groupnum <- max(mean_values_all$groupnum)
  # Function to fit linear model and ensure line passes through (12, 1)
  fit_line_through_point <- function(df, x0 = groupnum, y0 = 1) {
    # 创建一个新的变量，用于中心化数据
    df$groupnum_centered <- df$groupnum - x0
    df$mean_fixed_centered <- df$mean_fixed - y0

    # 拟合中心化后的数据，不需要截距项
    lm_fit <- lm(mean_fixed_centered ~ groupnum_centered + 0, data = df)

    # 斜率
    slope <- coef(lm_fit)[1]

    # 计算截距，使得直线通过 (12, 1)
    intercept <- y0 - slope * x0

    return(c(intercept, slope))
  }

  # Fit linear model for each cellname group
  models <- by(mean_values_all, mean_values_all$cellname, fit_line_through_point)

  # Extract intercepts and slopes from models
  intercepts <- sapply(models, function(m) m[1])
  slopes <- sapply(models, function(m) m[2])

  # Create a data frame with intercepts, slopes, and cellnames
  line_data <- data.frame(
    intercepts = intercepts,
    slopes = slopes,
    cellname = unique(mean_values_all$cellname)
  )



  p <- ggplot(mean_values_all, aes(x = groupnum, y = mean_fixed, color = cellname)) +
    geom_point(position = position_jitter(width = 0.2,height = 0.001), size = 3) +  # 添加抖动点
    scale_x_continuous(breaks = seq(min(mean_values_all$groupnum), max(mean_values_all$groupnum), by = 1)) +  # 设置横坐标以1为单位
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +  # 设置纵坐标以0.1为单位，并从0开始
    labs(title = "", x = "Dataset Number", y = "Accuracy") +  # 添加标题和坐标轴标签
    geom_abline(data = line_data, aes(intercept = intercepts, slope = slopes, color = cellname), size = 1) +
    theme_minimal() +
    scale_color_manual(values = colors) +  # 设置颜色
    mytheme2 # 使用 LOESS 方法进行曲线拟合，不显示标准误差带
  print(p)


  return(mean_values_all)

}

