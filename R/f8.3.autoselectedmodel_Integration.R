
autoselectedmodel_Integration <- function(mdlist,effect_cell_rate="rate"){

  p_value_for_model_selection=0.05

  calculate_OR_CI <- function(df, conf_level = 0.95) {

    coef_estimate <- df$Estimate
    std_error <-df$`Std. Error`
    pvalue <- df[, grepl("Pr", colnames(df))]
    # 计算对数几率（log odds）的估计值
    log_odds <- coef_estimate

    # 计算对数几率的标准误差
    se_log_odds <- std_error

    # 计算置信水平对应的Z分位数
    z_value <- qnorm(1 - (1 - conf_level) / 2)

    # 计算对数几率的置信区间
    lower_bound <- log_odds - z_value * se_log_odds
    upper_bound <- log_odds + z_value * se_log_odds

    # 计算对应的比率（OR）值
    OR <- exp(log_odds)

    # 返回结果
    result <- data.frame(
      Beta=coef_estimate,
      SD= std_error,
      OR = OR,
      CI_lower = exp(lower_bound),
      CI_upper = exp(upper_bound),
      Pvalue =pvalue
    )

    return(result)
  }


  #library("gtsummary")
  finalchosenmodel <- list()
  alltemp <<- c()


  for(i in 1:length(mdlist)){

    dfres <- mdlist[[i]][[2]]

    indices <- c()
    # 去除带有||的模型，无法进行后续的计算
    for (j in 1:length(mdlist[[i]][[1]])) {
      variables <- attr(terms(as.formula(formula(mdlist[[i]][[1]][[j]]))), "variables")[-1L]
      result <- any( grepl("\\|\\|", variables))
      # 检查当前变量是否包含 ||
      if (result) {
        # 如果包含，则将索引添加到 indices 中
        indices <- c(indices, j)
      }
    }



    dfres$difAIC <- dfres$AIC-min(dfres$AIC)
    dfres$AIC_wt <- exp(0-dfres$difAIC/2)
    dfres$AIC_prob <- dfres$AIC_wt/sum(dfres$AIC_wt)
    cumulative_sum <- cumsum(dfres$AIC_prob)

    # 3. 确定所需的元素数量n，使得累加和达到或超过0.95
    nformultiplemodel <- which.max(cumulative_sum >= 1-p_value_for_model_selection)
    id <- dfres[1:nformultiplemodel,]$id
    result <- id[!id %in% indices]
    selectedmodel <- mdlist[[i]][[1]][result]

    #selectedmodel.res <- MuMIn::model.avg(selectedmodel,beta = "partial.sd")
    print(paste(length(result),"Models were integrated"))
    dup=0


    if(length(selectedmodel)>1){

      tryCatch({
        # 尝试执行 model.avg 函数(具体函数见：model_avg_debug2.R)
        model_avg <- MuMIn::model.avg(selectedmodel, beta = "partial.sd")
        # 调用 summary 函数查看模型平均结果
        selectedmodel.res  <- summary(model_avg)

      }, error = function(e) {
        # 捕获到错误时执行的代码
        # 打印出错信息
        print(paste("Warning Message:", conditionMessage(e)))
        error_message <<-  conditionMessage(e)

        matches <- gregexpr("(?<= = )\\d+",error_message, perl = TRUE)
        dup <- unlist(regmatches(error_message, matches))

        dup <- as.numeric(dup)

        print(paste(dup,"Models were removed"))
        result <- result[!result %in% dup]

        selectedmodel <-selectedmodel[-dup]
        # 再次尝试执行 model.avg 函数
       #selectedmodel.res <- MuMIn::model.avg(selectedmodel, beta = "partial.sd")%>%   summary()

        model_avg <- MuMIn::model.avg(selectedmodel, beta = "partial.sd")
        selectedmodel.res  <- summary(model_avg)

      })


      selectedmodel.df <- as.data.frame(selectedmodel.res$coefmat.subset)
      selectedmodel.df <-selectedmodel.df[effect_cell_rate,]
      #########SUMMARY table construct
      or.df <- calculate_OR_CI(selectedmodel.df)
      # 使用基本R语法重写 mutate 和嵌套的 ifelse 条件
      or.df$significance <- ifelse(or.df$Pvalue <= 0.001, "***",
                                   ifelse(or.df$Pvalue <= 0.01, "**",
                                          ifelse(or.df$Pvalue <= 0.05, "*", " ")))


      or.df$Modelnum <- length(result)-length(dup)
      or.df$celltype <- dfres$celltype[1]
      alltemp <- rbind(alltemp,or.df)
      print(paste("The",i,"celltype intergation calculation were finished (",or.df$celltype[1] ,")"))



    }else{

      print(paste("The",i,"celltype calculation were not needed "))

      selectedmodel <- mdlist[[i]][[1]][id[1]] #直接取出来AIC值排序第一个模型

      selectedmodel.res <-   summary(selectedmodel[[1]])
      selectedmodel.df <- as.data.frame(selectedmodel.res$coefficients)
      selectedmodel.df <-selectedmodel.df[effect_cell_rate,]

      # 检查列名中是否包含 "p" 的字符串
      if (!any(grepl("P", colnames(selectedmodel.df)))) {
        # 将最后一列的值赋给 $pvalue 列
        selectedmodel.df$Pr <- selectedmodel.df[, ncol(selectedmodel.df)]
      }

      #########SUMMARY table construct #
      or.df <- calculate_OR_CI(selectedmodel.df)
      # 使用基本R语法重写 mutate 和嵌套的 ifelse 条件
      or.df$significance <- ifelse(or.df$Pvalue <= 0.001, "***",
                                   ifelse(or.df$Pvalue <= 0.01, "**",
                                          ifelse(or.df$Pvalue <= 0.05, "*", " ")))

      or.df$Modelnum <- length(result)-length(dup)
      or.df$celltype <- dfres$celltype[1]
      alltemp <- rbind(alltemp,or.df)



    }

    # 循环结束
  }
  return(alltemp)


  finalchosenmodel <- "Integrated model"

  cwas.data <- mdlist
  finalresult <- list(finalchosenmodel,alltemp,cwas.data)
  names(finalresult) <- c("Chosen_model","Chosen_model_info","All_models")

  return(finalresult)



}

