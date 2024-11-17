#' CWAS Model Calculation for Mixed Effects
#'
#' This function computes mixed effects models for phenotype data using various combinations of fixed effects and random intercepts/slopes.
#'
#' @param cwas.data A data frame containing the phenotype data and covariates.
#' @param selectedcelltype A character specifying the selected cell type for modeling.
#' @param link The link function to use in the model. Default is "binomial" for logistic regression.
#'
#' @return A list containing two elements:
#'   \item{fitall}{A list of all fitted model objects that were successfully computed.}
#'   \item{fitres_sorted}{A data frame containing summary statistics of the fitted models, sorted by AIC (Akaike Information Criterion).}
#'
#' @details
#' This function fits a series of mixed effects logistic regression models using the glmer function from the lme4 package. The models include combinations of fixed effects (rate, covariates) and random intercepts/slopes for different grouping structures (Subgroup, Group, Group/Subgroup). The function attempts to fit each model and captures successful results for further analysis.
#'
#' The function outputs a list of fitted models (`fitall`) and a summary data frame (`fitres_sorted`) containing information about each model including model ID, cell type, formula, AIC, and p-value (based on the rate coefficient).
#'
#' @export
#'
#' @examples cwas_model_cal_mix(cwas.test.data)
#'
cwas_model_cal_mix<- function(cwas.data,selectedcelltype="selectedcelltype",link="binomial"){

  #library(tidyverse)
  #library(lme4)
  #library(lmerTest)

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")
  #usethis::use_package(package="lme4",type="Import")



  fitres <- NULL
  celldf <- cwas.data

  if(celldf$Phenotype[1]%in%c(0,1)){}else{
    celldf$Phenotype <- ifelse(celldf$Phenotype=="Control",0,1)
  }




  fitall <- list()
  tryCatch({



    print(paste0(length(celldf$Phenotype)," Sample were enronlled in the model"))



    # 混合效应模型



    ## 固定斜率+随机截距模型

    fit11 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit12 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit13 <- tryCatch(lme4::glmer(Phenotype ~ rate+Group+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit14 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit15 <- tryCatch(lme4::glmer(Phenotype ~ rate +Group+ (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit16 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1|Group)+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit17 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup+ (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit18 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    fit21 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit22 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit23 <- tryCatch(lme4::glmer(Phenotype ~ rate+Group+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit24 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit25 <- tryCatch(lme4::glmer(Phenotype ~ rate +Group+ (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit26 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate|Group)+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit27 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup+ (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit28 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    fit31 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit32 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit33 <- tryCatch(lme4::glmer(Phenotype ~ rate+Group+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit34 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit35 <- tryCatch(lme4::glmer(Phenotype ~ rate +Group+ (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit36 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate||Group)+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit37 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup+ (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit38 <- tryCatch(lme4::glmer(Phenotype ~ rate+Subgroup + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)

    print(paste0("group1 ",selectedcelltype," fitting models have been calculated successful",""))

    #group2#


    ## 固定斜率+随机截距模型

    fit11_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit12_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit13_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Group+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit14_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit15_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 +Group+ (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit16_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (1|Group)+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit17_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup+ (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit18_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    fit21_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit22_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit23_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Group+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit24_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit25_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 +Group+ (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit26_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate|Group)+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit27_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup+ (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit28_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    fit31_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit32_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit33_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Group+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit34_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit35_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 +Group+ (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit36_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1 + (rate||Group)+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit37_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup+ (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit38_1 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Subgroup + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)





    print(paste0("group2 ",selectedcelltype," fitting models have been calculated successful",""))
    #group3#


    ## 固定斜率+随机截距模型

    fit11_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit12_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit13_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Group+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit14_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit15_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 +Group+ (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit16_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit17_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup+ (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit18_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    fit21_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit22_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit23_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Group+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit24_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit25_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit26_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit27_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup+ (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit28_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    fit31_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit32_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit33_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Group+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit34_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit35_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit36_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit37_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup+ (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit38_2 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    print(paste0("group3 ",selectedcelltype," fitting models have been calculated successful",""))

    #group4##



    ## 固定斜率+随机截距模型

    fit11_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit12_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit13_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Group+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit14_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit15_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 +Group+ (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit16_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (1|Group)+ (1|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit17_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup+ (1|Group), family=link,data=celldf), error = function(e) NULL)
    fit18_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup + (1|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)



    ## 随机斜率+随机截距模型

    fit21_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit22_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit23_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Group+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit24_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit25_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 +Group+ (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit26_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate|Group)+ (rate|Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit27_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup+ (rate|Group), family=link,data=celldf), error = function(e) NULL)
    fit28_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup + (rate|Group/Subgroup), family=link,data=celldf), error = function(e) NULL)

    ## 随机斜率+随机截距模型 不相关

    fit31_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit32_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit33_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Group+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit34_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit35_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 +Group+ (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit36_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2 + (rate||Group)+ (rate||Subgroup), family=link,data=celldf), error = function(e) NULL)
    fit37_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup+ (rate||Group), family=link,data=celldf), error = function(e) NULL)
    fit38_3 <- tryCatch(lme4::glmer(Phenotype ~ rate+Cov1+Cov2+Subgroup + (rate||Group/Subgroup), family=link,data=celldf), error = function(e) NULL)







    print(paste0("group4 ",selectedcelltype," fitting models have been calculated successful",""))

    fitall <- list(fit11  ,fit12 ,fit13  ,fit14  ,fit15  ,fit16 ,fit17  ,fit18  ,fit21 ,fit22 ,fit23 ,fit24 ,fit25 ,fit26  ,fit27 ,fit28   ,fit31  ,fit32 ,fit33 ,fit34  ,fit35  ,fit36   ,fit37 ,fit38   ,fit11_1 ,fit12_1 ,fit13_1 ,fit14_1 ,fit15_1 ,fit16_1 ,fit17_1  ,fit18_1  ,fit21_1 ,fit22_1 ,fit23_1 ,fit24_1 ,fit25_1  ,fit26_1  ,fit27_1 ,fit28_1 ,fit31_1 ,fit32_1 ,fit33_1 ,fit34_1 ,fit35_1 ,fit36_1  ,fit37_1 ,fit38_1 ,fit11_2 ,fit12_2 ,fit13_2 ,fit14_2 ,fit15_2 ,fit16_2   ,fit17_2  ,fit18_2  ,fit21_2 ,fit22_2 ,fit23_2 ,fit24_2 ,fit25_2  ,fit26_2  ,fit27_2 ,fit28_2 ,fit31_2 ,fit32_2 ,fit33_2 ,fit34_2 ,fit35_2  ,fit36_2  ,fit37_2 ,fit38_2 ,fit11_3 ,fit12_3 ,fit13_3   ,fit14_3 ,fit15_3   ,fit16_3   ,fit17_3   ,fit18_3 ,fit21_3  ,fit22_3 ,fit23_3 ,fit24_3  ,fit25_3  ,fit26_3 ,fit27_3 ,fit28_3 ,fit31_3  ,fit32_3 ,fit33_3 ,fit34_3 ,fit35_3  ,fit36_3 ,fit37_3 ,fit38_3)


    fitall <- Filter(Negate(is.null), fitall)
    fitres <- do.call(anova, fitall)
    fitres <- as.data.frame(fitres)
    i=1
    fitres$Pvalue=0
    fitres$formula=0
    fitres$id=seq(1:length(fitall))

    for(i in 1:length(fitall)){
      sum <- summary(fitall[[i]])
      coff <- as.data.frame(sum$coefficients)

      tryCatch({fitres$Pvalue[i] <- coff["rate", grepl("Pr", names(coff))]}, error = function(e) {
        # 如果出现错误，则执行备用操作
        fitres$Pvalue[i]  <<- coff["rate", grepl("value", names(coff))]
      })
      fitres$formula[i] <- paste(link,"link function:",as.character( sum[["call"]][["formula"]])[3])

    }

    fitres$model <- paste("model",seq(1:length(fitres$formula)))

    fitres$celltype <- selectedcelltype

    fitres_sorted <- fitres[order(fitres$AIC), ]

    fitres_sorted <- fitres_sorted[,c("id","model","celltype","formula","AIC","Pvalue")]

    print(paste0(length(fitres$AIC),selectedcelltype," fitting models have been calculated successful",""))

    finalresult <- list(fitall,fitres_sorted,cwas.data)
    names(finalresult) <- c("Models","Models_info","Rawdata")




    return(finalresult)  }, error = function(e) {
      # 捕获异常并输出错误信息
      cat("Error in iteration", i, ":", conditionMessage(e), "\n")
      return("NotSuccess")


    })



}
