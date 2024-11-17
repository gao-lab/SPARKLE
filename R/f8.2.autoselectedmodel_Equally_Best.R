
autoselectedmodel_Equally_Best <- function(mdlist,effect_cell_rate="rate"){

  p_value_for_signifcance=0.05

  p_value_for_model_selection=0.05

  calculate_OR_CI <- function(df, conf_level = 0.95) {

    coef_estimate <- df$Estimate
    std_error <-df$`Std. Error`
    pvalue <-  df[,grepl("Pr", names(df))]
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

    dfres$difAIC <- dfres$AIC-min(dfres$AIC)
    dfres$AIC_wt <- exp(0-dfres$difAIC/2)
    dfres$AIC_prob <- dfres$AIC_wt/sum(dfres$AIC_wt)
    cumulative_sum <- cumsum(dfres$AIC_prob)

    # 3. 确定所需的元素数量n，使得累加和达到或超过0.95
    nformultiplemodel <- which.max(cumulative_sum >= 1-p_value_for_model_selection)


    dfpvalue <- dfres[dfres$AIC_prob > 0.05,]
    dfpvalue$adjustedPvauel <- dfpvalue$Pvalue*length(dfpvalue$AIC)
    #p.adjust(dfpvalue$Pvalue,method =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") )
   # dfpvalue <- dfpvalue %>% dplyr::select_all() %>% dplyr::filter(adjustedPvauel<p_value_for_signifcance)
    # 使用基本R语法重写筛选操作
    dfpvalue <- subset(dfpvalue, adjustedPvauel < p_value_for_signifcance)


    if(length(dfpvalue$AIC)==0){##如果是空，则说明给定minAICnum中都不显著
      #outputmodelnum <- mdlist[[i]][[2]] %>%  dplyr::select_all() %>%  dplyr::arrange(AIC) %>%  dplyr::slice_min(order_by = AIC, n =1)
      # 使用基本R语法重写管道操作
      model_list <- mdlist[[i]][[2]]  # 从列表中提取模型对象
      # 按照 AIC 排序
      model_list <- model_list[order(model_list$AIC), ]
      # 提取 AIC 最小的行
      outputmodelnum <- model_list[which.min(model_list$AIC), , drop = FALSE]
    outputmodelnum <- outputmodelnum$id
      print("no significant model")
    }else{
      #outputmodelnum <- dfpvalue %>%  dplyr::select_all() %>%  dplyr::arrange(AIC) %>%  dplyr::slice_min(order_by = AIC, n =1)

      dfpvalue <- dfpvalue[order(dfpvalue$AIC), ]

      # 选择 AIC 列最小值所在的行
      outputmodelnum <- dfpvalue[which.min(dfpvalue$AIC), ]
      outputmodelnum <- outputmodelnum$id
      print("significant model selected")

    }


    finalchosenmodel[[i]] <- mdlist[[i]][[1]][[outputmodelnum]]


    #########SUMMARY table construct ##
    selectedmodel.res <- summary(mdlist[[i]][[1]][[outputmodelnum]])
    selectedmodel.df <- as.data.frame(selectedmodel.res$coefficients)
    selectedmodel.df <-selectedmodel.df[effect_cell_rate,]
    or.df <- calculate_OR_CI(selectedmodel.df)

    df <- or.df[1,]
    temp <- mdlist[[i]][[2]]
    temp <- temp[which(temp$id==outputmodelnum),]
    temp$OR <- df$OR
    temp$CI_lower <- df$CI_lower
    temp$CI_upper <- df$CI_upper
    temp$SD <- df$SD
    temp$Beta <- df$Beta
    temp$adjustedPvalue <- ifelse(length(dfpvalue$AIC)!=0,temp$Pvalue*length(dfpvalue$AIC),temp$Pvalue)

    print(temp$adjustedPvalue)

    # 使用基本R语法重写 mutate 和 ifelse 条件
    temp$significance <- ifelse(temp$adjustedPvalue <= 0.001, "***",
                                ifelse(temp$adjustedPvalue <= 0.01, "**",
                                       ifelse(temp$adjustedPvalue <= 0.05, "*", " ")))
    alltemp <- rbind(alltemp,temp)

    print(i)

  }

  cwas.data <- mdlist

  finalresult <- list(finalchosenmodel,alltemp,cwas.data)

  names(finalresult) <- c("Chosen_model","Chosen_model_info","All_models")

  return(finalresult)




}
