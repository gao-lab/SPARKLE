
interaction_2model_cal_v5<- function(cell2df,selectedcelltype="cellname",link="binomial"){
  interactionfitres <- NULL

  interactionfitall <- list()

  if(cell2df$Phenotype[1]%in%c(0,1)){}else{
    cell2df$Phenotype <- ifelse(cell2df$Phenotype=="Control",0,1)
  }




  tryCatch({

    #lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (1|Subgroup), family=link,data=cell2df)
    # library(tidyverse)
    # library(lme4)
    # library(lmerTest)

    print(paste0(length(cell2df$Phenotype)," Sample were enronlled in the model"))

    # 固定效应模型

    # 单纯线性模型，放弃加入了，因为明显效果都不好
    # interactionfit01 <- lm(Phenotype~rate,data = cell2df)
    # interactionfit02 <- lm(Phenotype~rate+Subgroup,data = cell2df)
    # interactionfit03 <- lm(Phenotype~rate+Group,data = cell2df)
    # interactionfit03 <- lm(Phenotype~rate+Subgroup+Group,data = cell2df)

    # interactionfit01 <- glm(Phenotype~rate,data = cell2df,family = link)
    # interactionfit02 <- glm(Phenotype~rate+Subgroup,data = cell2df,family = link)
    # interactionfit03 <- glm(Phenotype~rate+Group,data = cell2df,family = link)
    # interactionfit04 <- glm(Phenotype~rate+Subgroup+Group,data = cell2df,family = link)

    # 混合效应模型

    #group1##

    ## 固定斜率+随机截距模型

    interactionfit11 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit12 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit13 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit14 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit15 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit16 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit17 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit18 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    interactionfit21 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit22 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit23 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit24 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit25 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit26 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit27 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit28 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    interactionfit31 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit32 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit33 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit34 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit35 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit36 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit37 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit38 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    print(paste0("group1 ",selectedcelltype," interactionfitting models have been calculated successful",""))

    #group2###
    ## 固定斜率+随机截距模型

    interactionfit11_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit12_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit13_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit14_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit15_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit16_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit17_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit18_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    interactionfit21_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit22_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit23_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit24_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit25_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit26_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit27_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit28_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    interactionfit31_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit32_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit33_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit34_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit35_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit36_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit37_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit38_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)





    print(paste0("group2 ",selectedcelltype," interactionfitting models have been calculated successful",""))
    #group3###

    ## 固定斜率+随机截距模型

    interactionfit11_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit12_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit13_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit14_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit15_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit16_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit17_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit18_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    interactionfit21_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit22_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit23_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit24_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit25_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit26_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit27_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit28_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    interactionfit31_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit32_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit33_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit34_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit35_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit36_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit37_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit38_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    print(paste0("group3 ",selectedcelltype," interactionfitting models have been calculated successful",""))

    #group4###

    #交叉随机效应补充


    ## 固定斜率+随机截距模型

    interactionfit11_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit12_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit13_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Group+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit14_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit15_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 +Group+ (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit16_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit17_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup+ (1|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit18_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    interactionfit21_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit22_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit23_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Group+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit24_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit25_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit26_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit27_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup+ (rate|Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit28_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    interactionfit31_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit32_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit33_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Group+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit34_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit35_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit36_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=cell2df), error = function(e) NULL)
    interactionfit37_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup+ (rate||Group), family=link,data=cell2df), error = function(e) NULL)
    interactionfit38_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+rate2+rate*rate2+Cov1+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=cell2df), error = function(e) NULL)



    #result##



    print(paste0("group4 ",selectedcelltype," interactionfitting models have been calculated successful",""))

    interactionfitall <- list(interactionfit11   ,interactionfit12   ,interactionfit13  ,interactionfit14   ,interactionfit15   ,interactionfit16    ,interactionfit17  ,interactionfit18   ,interactionfit21   ,interactionfit22   ,interactionfit23  ,interactionfit24   ,interactionfit25   ,interactionfit26    ,interactionfit27  ,interactionfit28   ,interactionfit31   ,interactionfit32   ,interactionfit33  ,interactionfit34   ,interactionfit35   ,interactionfit36    ,interactionfit37  ,interactionfit38   ,interactionfit11_1   ,interactionfit12_1   ,interactionfit13_1  ,interactionfit14_1   ,interactionfit15_1   ,interactionfit16_1    ,interactionfit17_1  ,interactionfit18_1   ,interactionfit21_1   ,interactionfit22_1   ,interactionfit23_1  ,interactionfit24_1   ,interactionfit25_1   ,interactionfit26_1    ,interactionfit27_1  ,interactionfit28_1   ,interactionfit31_1   ,interactionfit32_1   ,interactionfit33_1  ,interactionfit34_1   ,interactionfit35_1   ,interactionfit36_1    ,interactionfit37_1  ,interactionfit38_1 ,interactionfit11_2   ,interactionfit12_2   ,interactionfit13_2  ,interactionfit14_2   ,interactionfit15_2   ,interactionfit16_2    ,interactionfit17_2  ,interactionfit18_2   ,interactionfit21_2   ,interactionfit22_2   ,interactionfit23_2  ,interactionfit24_2   ,interactionfit25_2   ,interactionfit26_2    ,interactionfit27_2  ,interactionfit28_2   ,interactionfit31_2   ,interactionfit32_2   ,interactionfit33_2  ,interactionfit34_2   ,interactionfit35_2   ,interactionfit36_2    ,interactionfit37_2  ,interactionfit38_2   ,interactionfit11_3   ,interactionfit12_3   ,interactionfit13_3  ,interactionfit14_3   ,interactionfit15_3   ,interactionfit16_3    ,interactionfit17_3  ,interactionfit18_3   ,interactionfit21_3   ,interactionfit22_3   ,interactionfit23_3  ,interactionfit24_3   ,interactionfit25_3   ,interactionfit26_3    ,interactionfit27_3  ,interactionfit28_3   ,interactionfit31_3   ,interactionfit32_3   ,interactionfit33_3  ,interactionfit34_3   ,interactionfit35_3   ,interactionfit36_3    ,interactionfit37_3  ,interactionfit38_3  )


    interactionfitall <- Filter(Negate(is.null), interactionfitall)

    interactionfitres <- do.call(anova, interactionfitall)
    # interactionfitres <- anova(interactionfit11   ,interactionfit12   ,interactionfit13  ,interactionfit14   ,interactionfit15   ,interactionfit16    ,interactionfit17  ,interactionfit18   ,interactionfit21   ,interactionfit22   ,interactionfit23  ,interactionfit24   ,interactionfit25   ,interactionfit26    ,interactionfit27  ,interactionfit28   ,interactionfit31   ,interactionfit32   ,interactionfit33  ,interactionfit34   ,interactionfit35   ,interactionfit36    ,interactionfit37  ,interactionfit38   ,interactionfit11_1   ,interactionfit12_1   ,interactionfit13_1  ,interactionfit14_1   ,interactionfit15_1   ,interactionfit16_1    ,interactionfit17_1  ,interactionfit18_1   ,interactionfit21_1   ,interactionfit22_1   ,interactionfit23_1  ,interactionfit24_1   ,interactionfit25_1   ,interactionfit26_1    ,interactionfit27_1  ,interactionfit28_1   ,interactionfit31_1   ,interactionfit32_1   ,interactionfit33_1  ,interactionfit34_1   ,interactionfit35_1   ,interactionfit36_1    ,interactionfit37_1  ,interactionfit38_1 ,interactionfit11_2   ,interactionfit12_2   ,interactionfit13_2  ,interactionfit14_2   ,interactionfit15_2   ,interactionfit16_2    ,interactionfit17_2  ,interactionfit18_2   ,interactionfit21_2   ,interactionfit22_2   ,interactionfit23_2  ,interactionfit24_2   ,interactionfit25_2   ,interactionfit26_2    ,interactionfit27_2  ,interactionfit28_2   ,interactionfit31_2   ,interactionfit32_2   ,interactionfit33_2  ,interactionfit34_2   ,interactionfit35_2   ,interactionfit36_2    ,interactionfit37_2  ,interactionfit38_2   ,interactionfit11_3   ,interactionfit12_3   ,interactionfit13_3  ,interactionfit14_3   ,interactionfit15_3   ,interactionfit16_3    ,interactionfit17_3  ,interactionfit18_3   ,interactionfit21_3   ,interactionfit22_3   ,interactionfit23_3  ,interactionfit24_3   ,interactionfit25_3   ,interactionfit26_3    ,interactionfit27_3  ,interactionfit28_3   ,interactionfit31_3   ,interactionfit32_3   ,interactionfit33_3  ,interactionfit34_3   ,interactionfit35_3   ,interactionfit36_3    ,interactionfit37_3  ,interactionfit38_3   )


    #result##
    interactionfitres <- as.data.frame(interactionfitres)
    i=1
    interactionfitres$Pvalue=0
    interactionfitres$formula=0
    interactionfitres$id=seq(1:length(interactionfitall))

    for(i in 1:length(interactionfitall)){
      sum <- summary(interactionfitall[[i]])
      coff <- as.data.frame(sum$coefficients)
      interactionfitres$Pvalue[i] <- coff["rate:rate2",grepl("Pr", names(coff))]
      interactionfitres$formula[i] <- as.character( sum[["call"]][["formula"]])[3]

    }

    interactionfitres$model <- paste("model",seq(1:length(interactionfitres$formula)))
    interactionfitres$celltype <- selectedcelltype

    interactionfitres_sorted <- interactionfitres[order(interactionfitres$AIC), ]

    interactionfitres_sorted$significance <- ifelse(interactionfitres_sorted$Pvalue <= 0.001, "***",
                                                    ifelse(interactionfitres_sorted$Pvalue <= 0.01, "**",
                                                           ifelse(interactionfitres_sorted$Pvalue <= 0.05, "*", "-")))




    print(paste0(length(interactionfitres$AIC),selectedcelltype," interactionfitting models have been calculated successful",""))


  }, error = function(e) {
    # 捕获异常并输出错误信息
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
  })

  finalresult <- list(interactionfitall,interactionfitres_sorted,cell2df)


  names(finalresult) <- c("Models","Models_info","Rawdata")



  return(finalresult)

}


