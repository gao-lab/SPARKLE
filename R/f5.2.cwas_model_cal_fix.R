#' CWAS Model Calculation for Fixed Effects
#'
#' This function computes fixed effects models for phenotype data using generalized linear models (GLMs) with different link functions.
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
#' This function fits a series of fixed effects models using generalized linear models (GLMs) from the stats package. The models include combinations of fixed effects (rate, covariates) with different grouping structures (Subgroup, Group, Subgroup/Group). The function attempts to fit each model and captures successful results for further analysis.
#'
#' The function outputs a list of fitted models (`fitall`) and a summary data frame (`fitres_sorted`) containing information about each model including model ID, cell type, formula, AIC, and p-value (based on the rate coefficient).
#'
#' @export
#'
#' @examples cwas_model_cal_fix(cwas.test.data)
cwas_model_cal_fix<- function(cwas.data,selectedcelltype="selectedcelltype_name",link="binomial"){

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")
  #usethis::use_package(package="lme4",type="Import")
  #usethis::use_package(package="stats",type="Import")


  fitres <- NULL

  celldf <- cwas.data
  celldf$Phenotype <- ifelse(celldf$Phenotype=="Control",0,1)
  fitall <- list()


  # library(tidyverse)
  # library(lme4)
  # library(lmerTest)


  print(paste0(length(celldf$Phenotype)," Sample were enronlled in the model"))

  # 固定效应模型
  if(link%in%c("binomial","binomial","poisson","Gamma","inverse.binomial","quasibinomial")){


    fit0001 <-   tryCatch(stats::glm(Phenotype~rate,data = celldf,family = link), error = function(e) NULL)
    fit0002 <-   tryCatch(stats::glm(Phenotype~rate+Subgroup,data = celldf,family = link), error = function(e) NULL)
    fit0003 <-   tryCatch(stats::glm(Phenotype~rate+Group,data = celldf,family = link), error = function(e) NULL)
    fit0004 <-   tryCatch(stats::glm(Phenotype~rate+Subgroup+Group,data = celldf,family = link), error = function(e) NULL)
    fit0001_1 <- tryCatch(stats::glm(Phenotype~rate+cov1,data = celldf,family =link), error = function(e) NULL)
    fit0002_1 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+cov1,data = celldf,family = link), error = function(e) NULL)
    fit0003_1 <- tryCatch(stats::glm(Phenotype~rate+Group+cov1,data = celldf,family = link), error = function(e) NULL)
    fit0004_1 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+Group+cov1,data = celldf,family = link), error = function(e) NULL)
    fit0001_2 <- tryCatch(stats::glm(Phenotype~rate+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0002_2 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0003_2 <- tryCatch(stats::glm(Phenotype~rate+Group+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0004_2 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+Group+cov2,data = celldf,family =link), error = function(e) NULL)
    fit0001_3 <- tryCatch(stats::glm(Phenotype~rate+cov1+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0002_3 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+cov1+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0003_3 <- tryCatch(stats::glm(Phenotype~rate+Group+cov1+cov2,data = celldf,family = link), error = function(e) NULL)
    fit0004_3 <- tryCatch(stats::glm(Phenotype~rate+Subgroup+Group+cov1+cov2,data = celldf,family =link), error = function(e) NULL)



  }else{


    print("No link function , please check stats::glm package")

  }





  fitall <- list(fit0001  ,fit0002  ,fit0003  ,fit0004  ,fit0001_1   ,fit0002_1    ,fit0003_1   ,fit0004_1   ,fit0001_2   ,fit0002_2   ,fit0003_2    ,fit0004_2   ,fit0001_3    ,fit0002_3   ,fit0003_3    ,fit0004_3)

  fitall <- Filter(Negate(is.null), fitall)

  fitres <- anova(fit0001  ,fit0002  ,fit0003  ,fit0004  ,fit0001_1   ,fit0002_1    ,fit0003_1   ,fit0004_1   ,fit0001_2   ,fit0002_2   ,fit0003_2    ,fit0004_2   ,fit0001_3    ,fit0002_3   ,fit0003_3    ,fit0004_3)

  fitres <- as.data.frame(fitres)


  if(length(fitall)==1){fitres <- fitres[2,]}

  i=1
  fitres$Pvalue=0
  fitres$formula=0
  fitres$id=seq(1:length(fitall))

  for(i in 1:length(fitall)){
    sum <- summary(fitall[[i]])
    coff <- as.data.frame(sum$coefficients)
    fitres$Pvalue[i] <- coff["rate",grepl("Pr", names(coff))]
    fitres$formula[i] <- paste(link,"link function:",as.character( sum[["call"]][["formula"]])[3])

  }

  fitres$model <- paste("model",seq(1:length(fitres$formula)))
  fitres$celltype <- selectedcelltype


  AIC_values <- sapply(fitall, function(model) AIC(model))
  fitres$AIC <- AIC_values
  fitres_sorted <- fitres[order(fitres$AIC), ]

  fitres_sorted <- fitres_sorted[,c("id","model","celltype","formula","AIC","Pvalue")]

  print(paste0(length(fitres$AIC),selectedcelltype," fitting models have been calculated successful",""))

  finalresult <- list(fitall,fitres_sorted,cwas.data)
  names(finalresult) <- c("Models","Models_info","Rawdata")



  return(finalresult)

}

