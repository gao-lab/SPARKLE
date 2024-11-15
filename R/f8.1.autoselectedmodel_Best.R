

autoselectedmodel_Best <- function(mdlist,minAICnum=1,effect_cell_rate="rate"){

  pvalue <- 0.05
  # library(tidyverse)
  # library(lme4)
  # library(lmerTest)

  calculate_OR_CI <- function(df, conf_level = 0.95) {

    coef_estimate <- df[,1]
    std_error <-df[,2]
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

    # 从列表中提取模型数据框，并按照 AIC 排序，取出前 minAICnum 个模型
    dfres <- mdlist[[i]][[2]]
    dfres <- dfres[order(dfres$AIC), ]
    dfres <- dfres[1:minAICnum, ]

    # 从排序后的模型中筛选出 Pvalue 小于 pvalue 的模型
    dfpvalue <- dfres[dfres$Pvalue < pvalue, , drop = FALSE]

    # 检查 dfpvalue 是否为空
    if (nrow(dfpvalue) == 0) {
      # 如果为空，则选择 AIC 最小的模型作为输出模型
      outputmodelnum <- dfres[which.min(dfres$AIC), "id"]
      cat("no significant model\n")
    } else {
      # 如果不为空，则选择符合条件的 AIC 最小的模型作为输出模型
      dfpvalue <- dfpvalue[order(dfpvalue$AIC), ]
      outputmodelnum <- dfpvalue[1, "id"]
      cat("significant model selected\n")
    }

    # 输出选择的模型的 id


    finalchosenmodel[[i]] <- mdlist[[i]][[1]][[outputmodelnum]]


    #SUMMARY table construct ##


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
    temp$significance <- ifelse(temp$Pvalue <= 0.001, "***",
                                ifelse(temp$Pvalue <= 0.01, "**",
                                       ifelse(temp$Pvalue <= 0.05, "*", " ")))


    alltemp <- rbind(alltemp,temp)


    print(i)

  }

  cwas.data <- mdlist

  finalresult <- list(finalchosenmodel,alltemp,cwas.data)
  names(finalresult) <- c("Chosen_model","Chosen_model_info","All_models")

  return(finalresult)

}



