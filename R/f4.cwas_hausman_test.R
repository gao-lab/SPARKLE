#' Conducts Hausman Test for each Celltype in  data
#'
#' @param cwas.data Data frame containing  data
#' @param variable Name of the variable to test (dependent variable)
#' @param family Type of family for GLM (default is "binomial")
#'
#' @return A data frame with results of Hausman test for each Celltype
#'
#' @export
#'
cwas_hausman_test_all <- function(cwas.data, variable, family= "binomial"){

  hausman_test <- function(cwas.data, variable, family= "binomial",cellname="cellnames") {

    cwas.data$variable <- cwas.data[,variable]

    if (cwas.data$Phenotype[1] %in% c(0, 1)) {
    }
    else {
      cwas.data$Phenotype <- ifelse(cwas.data$Phenotype == "Control",
                                 0, 1)
    }

    # 拟合固定效应模型
    fixed_model <- stats::glm(Phenotype ~ rate + variable, family = family, data = cwas.data)

    # 拟合随机效应模型
    random_model <- lme4::glmer(Phenotype ~ rate + (1 | variable), family = family, data = cwas.data)

    # 提取系数和协方差矩阵
    fixed_coef <- coef(summary(fixed_model))[, "Estimate"]
    random_coef <- lme4::fixef(random_model)
    fixed_cov <- stats::vcov(fixed_model)
    random_cov <- as.matrix(stats::vcov(random_model))

    # 计算差异
    diff_coef <- random_coef - fixed_coef[1:length(random_coef)]
    cov_diff <- random_cov - fixed_cov[1:length(random_coef), 1:length(random_coef)]

    # 计算豪斯曼统计量
    hausman_stat <- t(diff_coef) %*% solve(cov_diff) %*% diff_coef
    p_value <- pchisq(hausman_stat, df = length(diff_coef), lower.tail = FALSE)

    # 输出结果
    if (p_value < 0.05) {
      model_choice <- "Fixed Effects Model"
    } else {
      model_choice <- "Random Effects Model"
    }

    # cat("Hausman Statistic: ", hausman_stat, "\n")
    # cat("P-value: ", p_value, "\n")
    # cat("Model Selection: ", model_choice, "\n")

    result <- data.frame(Celltype=cellname,Hausman_Statistic=hausman_stat,Pvalue=p_value,Model_Selection=model_choice)

    return(result)
  }


  cellnamesall <- unique(cwas.data$Celltype)
  dfall <- c()
  for(i in 1:length(cellnamesall)){

    #cwas.data.selected <-  cwas.data%>% dplyr::select_all() %>% dplyr::filter(Celltype==cellnamesall[i])
    cwas.data.selected <- subset(cwas.data, Celltype == cellnamesall[i])

    df <- tryCatch({
      hausman_test(cwas.data = cwas.data.selected, family = family, variable = variable, cellname = cellnamesall[i])
    }, error = function(e) {
      message("Error in hausman_test for celltype: ", cellnamesall[i], ". Skipping this celltype.")
      return(NULL) # 返回 NULL 以便后续处理
    })

    # 如果 df 不为 NULL，则将其绑定到 dfall
    if (!is.null(df)) {
      dfall <- rbind(dfall, df)
    }


  }
  myattr <- attributes(cwas.data)

  dfall$Varible <- paste0(variable," : ",myattr[[which(names(myattr)==variable)]])
  #print(myattr[[which(names(myattr)==variable)]])
  return(dfall)

}

