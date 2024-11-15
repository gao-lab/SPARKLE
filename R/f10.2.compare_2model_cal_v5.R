

compare_2model_cal_v5<- function(cell2df,selectedcelltype="cellname",link="binomial"){
  comparefitres <- NULL

  comparefitall <- list()

  if(cell2df$Phenotype[1]%in%c(0,1)){}else{
    cell2df$Phenotype <- ifelse(cell2df$Phenotype=="Control",0,1)
  }




  tryCatch({

    # library(tidyverse)
    # library(lme4)
    # library(lmerTest)

    print(paste0(length(cell2df$Phenotype)," Sample were enronlled in the model"))

    # 固定效应模型

    # 单纯线性模型，放弃加入了，因为明显效果都不好
    # comparefit01 <- lm(Phenotype~rate,data = cell2df)
    # comparefit02 <- lm(Phenotype~rate+Subgroup,data = cell2df)
    # comparefit03 <- lm(Phenotype~rate+Group,data = cell2df)
    # comparefit03 <- lm(Phenotype~rate+Subgroup+Group,data = cell2df)

    # comparefit01 <- glm(Phenotype~rate,data = cell2df,family = link)
    # comparefit02 <- glm(Phenotype~rate+Subgroup,data = cell2df,family = link)
    # comparefit03 <- glm(Phenotype~rate+Group,data = cell2df,family = link)
    # comparefit04 <- glm(Phenotype~rate+Subgroup+Group,data = cell2df,family = link)

    # 混合效应模型

    #group1##

    ## 固定斜率+随机截距模型

    comparefit11 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit12 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit13 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit14 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit15 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit16 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit17 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit18 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    comparefit21 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit22 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit23 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit24 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit25 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit26 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit27 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit28 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    comparefit31 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit32 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit33 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit34 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit35 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit36 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit37 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit38 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    print(paste0("group1 ",selectedcelltype," comparefitting models have been calculated successful",""))

    #group2##
    ## 固定斜率+随机截距模型

    comparefit11_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit12_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit13_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit14_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit15_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit16_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit17_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit18_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    comparefit21_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit22_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit23_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit24_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit25_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit26_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit27_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit28_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    comparefit31_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit32_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit33_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit34_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit35_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit36_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit37_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit38_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)





    print(paste0("group2 ",selectedcelltype," comparefitting models have been calculated successful",""))
    #group3##

    ## 固定斜率+随机截距模型

    comparefit11_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit12_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit13_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit14_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit15_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit16_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit17_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit18_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    comparefit21_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit22_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit23_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit24_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit25_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit26_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit27_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit28_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    comparefit31_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit32_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit33_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit34_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit35_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit36_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit37_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit38_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    print(paste0("group3 ",selectedcelltype," comparefitting models have been calculated successful",""))

    #group4###
    #交叉随机效应补充


    ## 固定斜率+随机截距模型

    comparefit11_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit12_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit13_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit14_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit15_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit16_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit17_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit18_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    comparefit21_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit22_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit23_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit24_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit25_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit26_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit27_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit28_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    comparefit31_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit32_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit33_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit34_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit35_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit36_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    comparefit37_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    comparefit38_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+Cov1+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    #result###



    print(paste0("group4 ",selectedcelltype," comparefitting models have been calculated successful",""))

    comparefitall <- list(comparefit11   ,comparefit12   ,comparefit13  ,comparefit14   ,comparefit15   ,comparefit16    ,comparefit17  ,comparefit18   ,comparefit21   ,comparefit22   ,comparefit23  ,comparefit24   ,comparefit25   ,comparefit26    ,comparefit27  ,comparefit28   ,comparefit31   ,comparefit32   ,comparefit33  ,comparefit34   ,comparefit35   ,comparefit36    ,comparefit37  ,comparefit38   ,comparefit11_1   ,comparefit12_1   ,comparefit13_1  ,comparefit14_1   ,comparefit15_1   ,comparefit16_1    ,comparefit17_1  ,comparefit18_1   ,comparefit21_1   ,comparefit22_1   ,comparefit23_1  ,comparefit24_1   ,comparefit25_1   ,comparefit26_1    ,comparefit27_1  ,comparefit28_1   ,comparefit31_1   ,comparefit32_1   ,comparefit33_1  ,comparefit34_1   ,comparefit35_1   ,comparefit36_1    ,comparefit37_1  ,comparefit38_1 ,comparefit11_2   ,comparefit12_2   ,comparefit13_2  ,comparefit14_2   ,comparefit15_2   ,comparefit16_2    ,comparefit17_2  ,comparefit18_2   ,comparefit21_2   ,comparefit22_2   ,comparefit23_2  ,comparefit24_2   ,comparefit25_2   ,comparefit26_2    ,comparefit27_2  ,comparefit28_2   ,comparefit31_2   ,comparefit32_2   ,comparefit33_2  ,comparefit34_2   ,comparefit35_2   ,comparefit36_2    ,comparefit37_2  ,comparefit38_2   ,comparefit11_3   ,comparefit12_3   ,comparefit13_3  ,comparefit14_3   ,comparefit15_3   ,comparefit16_3    ,comparefit17_3  ,comparefit18_3   ,comparefit21_3   ,comparefit22_3   ,comparefit23_3  ,comparefit24_3   ,comparefit25_3   ,comparefit26_3    ,comparefit27_3  ,comparefit28_3   ,comparefit31_3   ,comparefit32_3   ,comparefit33_3  ,comparefit34_3   ,comparefit35_3   ,comparefit36_3    ,comparefit37_3  ,comparefit38_3  )


    comparefitall <- Filter(Negate(is.null), comparefitall)

    comparefitres <- do.call(anova, comparefitall)
    # comparefitres <- anova(comparefit11   ,comparefit12   ,comparefit13  ,comparefit14   ,comparefit15   ,comparefit16    ,comparefit17  ,comparefit18   ,comparefit21   ,comparefit22   ,comparefit23  ,comparefit24   ,comparefit25   ,comparefit26    ,comparefit27  ,comparefit28   ,comparefit31   ,comparefit32   ,comparefit33  ,comparefit34   ,comparefit35   ,comparefit36    ,comparefit37  ,comparefit38   ,comparefit11_1   ,comparefit12_1   ,comparefit13_1  ,comparefit14_1   ,comparefit15_1   ,comparefit16_1    ,comparefit17_1  ,comparefit18_1   ,comparefit21_1   ,comparefit22_1   ,comparefit23_1  ,comparefit24_1   ,comparefit25_1   ,comparefit26_1    ,comparefit27_1  ,comparefit28_1   ,comparefit31_1   ,comparefit32_1   ,comparefit33_1  ,comparefit34_1   ,comparefit35_1   ,comparefit36_1    ,comparefit37_1  ,comparefit38_1 ,comparefit11_2   ,comparefit12_2   ,comparefit13_2  ,comparefit14_2   ,comparefit15_2   ,comparefit16_2    ,comparefit17_2  ,comparefit18_2   ,comparefit21_2   ,comparefit22_2   ,comparefit23_2  ,comparefit24_2   ,comparefit25_2   ,comparefit26_2    ,comparefit27_2  ,comparefit28_2   ,comparefit31_2   ,comparefit32_2   ,comparefit33_2  ,comparefit34_2   ,comparefit35_2   ,comparefit36_2    ,comparefit37_2  ,comparefit38_2   ,comparefit11_3   ,comparefit12_3   ,comparefit13_3  ,comparefit14_3   ,comparefit15_3   ,comparefit16_3    ,comparefit17_3  ,comparefit18_3   ,comparefit21_3   ,comparefit22_3   ,comparefit23_3  ,comparefit24_3   ,comparefit25_3   ,comparefit26_3    ,comparefit27_3  ,comparefit28_3   ,comparefit31_3   ,comparefit32_3   ,comparefit33_3  ,comparefit34_3   ,comparefit35_3   ,comparefit36_3    ,comparefit37_3  ,comparefit38_3   )


    #result##

    comparefitres <- as.data.frame(comparefitres)
    i=1
    comparefitres$Pvalue=0
    comparefitres$formula=0
    comparefitres$id=seq(1:length(comparefitall))

    for(i in 1:length(comparefitall)){
      sum <- summary(comparefitall[[i]])
      coff <- as.data.frame(sum$coefficients)
      comparefitres$Pvalue[i] <- coff["rate",grepl("Pr", names(coff))]
      comparefitres$formula[i] <- as.character( sum[["call"]][["formula"]])[3]

    }

    comparefitres$model <- paste("model",seq(1:length(comparefitres$formula)))
    comparefitres$celltype <- selectedcelltype

    comparefitres_sorted <- comparefitres[order(comparefitres$AIC), ]



    # 使用 ifelse 函数根据 Pvalue 的值设置 significance 列的值
    comparefitres_sorted$significance <- ifelse(comparefitres_sorted$Pvalue <= 0.001, "***",
                                                ifelse(comparefitres_sorted$Pvalue <= 0.01, "**",
                                                       ifelse(comparefitres_sorted$Pvalue <= 0.05, "*", "-")))



    print(paste0(length(comparefitres$AIC),selectedcelltype," comparefitting models have been calculated successful",""))


  }, error = function(e) {
    # 捕获异常并输出错误信息
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
  })


  finalresult <- list(comparefitall,comparefitres_sorted,cell2df)


  names(finalresult) <- c("Models","Models_info","Rawdata")



  return(finalresult)

}

