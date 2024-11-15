cwas_model_cal_mix2 <- function (cwas.data, selectedcelltype = "selectedcelltype",
          link = "binomial")
{
  fitres <- NULL
  celldf <- cwas.data
  if (celldf$Phenotype[1] %in% c(0, 1)) {
  }
  else {
    celldf$Phenotype <- ifelse(celldf$Phenotype == "Control",
                               0, 1)
  }
  fitall <- list()
  tryCatch({
    print(paste0(length(celldf$Phenotype), " Sample were enronlled in the model"))

    fit01 <- tryCatch(stats::glm(Phenotype ~ rate , family = link, data = celldf), error = function(e) NULL)

    fit02 <- tryCatch(stats::glm(Phenotype ~ rate+ Group, family = link, data = celldf), error = function(e) NULL)

    fit03 <- tryCatch(stats::glm(Phenotype ~ rate +Subgroup, family = link, data = celldf), error = function(e) NULL)

    fit04 <- tryCatch(stats::glm(Phenotype ~ rate+ Group+Subgroup, family = link, data = celldf), error = function(e) NULL)

    fit11 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1 |
                                                        Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit12 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1 |
                                                        Group), family = link, data = celldf), error = function(e) NULL)
    fit13 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (1 | Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit14 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1 |
                                                        Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit15 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (1 | Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit16 <- tryCatch(lme4::glmer(Phenotype ~ rate + (1 |
                                                        Group) + (1 | Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit17 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (1 | Group), family = link, data = celldf), error = function(e) NULL)
    fit18 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (1 | Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit21 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate |
                                                        Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit22 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate |
                                                        Group), family = link, data = celldf), error = function(e) NULL)
    fit23 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (rate | Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit24 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate |
                                                        Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit25 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (rate | Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit26 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate |
                                                        Group) + (rate | Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit27 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (rate | Group), family = link, data = celldf), error = function(e) NULL)
    fit28 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (rate | Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit31 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate ||
                                                        Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit32 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate ||
                                                        Group), family = link, data = celldf), error = function(e) NULL)
    fit33 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (rate || Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit34 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate ||
                                                        Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit35 <- tryCatch(lme4::glmer(Phenotype ~ rate + Group +
                                    (rate || Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit36 <- tryCatch(lme4::glmer(Phenotype ~ rate + (rate ||
                                                        Group) + (rate || Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    fit37 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (rate || Group), family = link, data = celldf),
                      error = function(e) NULL)
    fit38 <- tryCatch(lme4::glmer(Phenotype ~ rate + Subgroup +
                                    (rate || Group/Subgroup), family = link, data = celldf),
                      error = function(e) NULL)
    print(paste0("group1 ", selectedcelltype, " fitting models have been calculated successful",
                 ""))
    fit11_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (1 | Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit12_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (1 | Group), family = link, data = celldf), error = function(e) NULL)
    fit13_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit14_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (1 | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit15_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (1 | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit16_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (1 | Group) + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit17_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (1 | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit18_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (1 | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit21_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit22_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate | Group), family = link, data = celldf), error = function(e) NULL)
    fit23_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (rate | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit24_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit25_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit26_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate | Group) + (rate | Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit27_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (rate | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit28_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit31_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate || Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit32_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate || Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit33_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (rate || Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit34_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate || Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit35_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Group + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit36_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      (rate || Group) + (rate || Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit37_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (rate || Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit38_1 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Subgroup + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    print(paste0("group2 ", selectedcelltype, " fitting models have been calculated successful",
                 ""))
    fit11_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (1 | Subgroup), family = link, data = celldf), error = function(e) NULL)
    fit12_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (1 | Group), family = link, data = celldf), error = function(e) NULL)
    fit13_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit14_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (1 | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit15_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (1 | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit16_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (1 | Group) + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit17_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (1 | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit18_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (1 | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit21_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit22_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate | Group), family = link, data = celldf), error = function(e) NULL)
    fit23_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (rate | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit24_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit25_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit26_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate | Group) + (rate | Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit27_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (rate | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit28_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit31_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate || Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit32_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate || Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit33_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (rate || Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit34_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate || Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit35_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Group + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit36_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      (rate || Group) + (rate || Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit37_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (rate || Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit38_2 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov2 +
                                      Subgroup + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    print(paste0("group3 ", selectedcelltype, " fitting models have been calculated successful",
                 ""))
    fit11_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit12_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (1 | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit13_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (1 | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit14_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (1 | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit15_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (1 | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit16_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (1 | Group) + (1 | Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit17_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (1 | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit18_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (1 | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit21_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate | Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit22_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate | Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit23_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (rate | Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit24_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate | Group/Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit25_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit26_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate | Group) + (rate | Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit27_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (rate | Group), family = link,
                                    data = celldf), error = function(e) NULL)
    fit28_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (rate | Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit31_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate || Subgroup), family = link, data = celldf),
                        error = function(e) NULL)
    fit32_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate || Group), family = link, data = celldf),
                        error = function(e) NULL)
    fit33_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (rate || Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit34_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit35_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Group + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit36_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + (rate || Group) + (rate || Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    fit37_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (rate || Group), family = link,
                                    data = celldf), error = function(e) NULL)
    fit38_3 <- tryCatch(lme4::glmer(Phenotype ~ rate + Cov1 +
                                      Cov2 + Subgroup + (rate || Group/Subgroup), family = link,
                                    data = celldf), error = function(e) NULL)
    print(paste0("group4 ", selectedcelltype, " fitting models have been calculated successful",
                 ""))

    getfixfomula <- function(model1){
      formula <- summary(model1)
      formula <- paste0(link," fix effect model:",as.character(formula$terms)[3])
      return(formula)
    }
    getpvalue <- function(model1){
      pvalue <- summary(model1)
      pvalue <- pvalue$coefficients["rate","Pr(>|z|)"]
      return(pvalue)
    }

    fitall2 <- list(fit01,fit02,fit03,fit04,fit11, fit12, fit13, fit14, fit15, fit16,
                   fit17, fit18, fit21, fit22, fit23, fit24, fit25,
                   fit26, fit27, fit28, fit31, fit32, fit33, fit34,
                   fit35, fit36, fit37, fit38, fit11_1, fit12_1, fit13_1,
                   fit14_1, fit15_1, fit16_1, fit17_1, fit18_1, fit21_1,
                   fit22_1, fit23_1, fit24_1, fit25_1, fit26_1, fit27_1,
                   fit28_1, fit31_1, fit32_1, fit33_1, fit34_1, fit35_1,
                   fit36_1, fit37_1, fit38_1, fit11_2, fit12_2, fit13_2,
                   fit14_2, fit15_2, fit16_2, fit17_2, fit18_2, fit21_2,
                   fit22_2, fit23_2, fit24_2, fit25_2, fit26_2, fit27_2,
                   fit28_2, fit31_2, fit32_2, fit33_2, fit34_2, fit35_2,
                   fit36_2, fit37_2, fit38_2, fit11_3, fit12_3, fit13_3,
                   fit14_3, fit15_3, fit16_3, fit17_3, fit18_3, fit21_3,
                   fit22_3, fit23_3, fit24_3, fit25_3, fit26_3, fit27_3,
                   fit28_3, fit31_3, fit32_3, fit33_3, fit34_3, fit35_3,
                   fit36_3, fit37_3, fit38_3)
    fitall2 <- Filter(Negate(is.null), fitall2)


    fitfix <- list(fit01,fit02,fit03,fit04)
    fitfix <- Filter(Negate(is.null), fitfix)
    if(length(fitfix)>0){
      # fitresfix <- do.call(anova, fitfix)
      # fitresfix <- as.data.frame(fitresfix)
      formula.result <- unlist(lapply(fitfix, getfixfomula))
      AIC.result <- unlist(lapply(fitfix, AIC))
      pvalue.result <- unlist(lapply(fitfix, getpvalue))
      fitresfixall <- data.frame(id=paste0(0,seq(1:length(fitfix))),model=paste0("model 0",seq(1:length(fitfix))),celltype=selectedcelltype,
                                 formula=  formula.result , AIC=AIC.result,Pvalue=pvalue.result  )
    }


    fitall <- list(fit11, fit12, fit13, fit14, fit15, fit16,
                   fit17, fit18, fit21, fit22, fit23, fit24, fit25,
                   fit26, fit27, fit28, fit31, fit32, fit33, fit34,
                   fit35, fit36, fit37, fit38, fit11_1, fit12_1, fit13_1,
                   fit14_1, fit15_1, fit16_1, fit17_1, fit18_1, fit21_1,
                   fit22_1, fit23_1, fit24_1, fit25_1, fit26_1, fit27_1,
                   fit28_1, fit31_1, fit32_1, fit33_1, fit34_1, fit35_1,
                   fit36_1, fit37_1, fit38_1, fit11_2, fit12_2, fit13_2,
                   fit14_2, fit15_2, fit16_2, fit17_2, fit18_2, fit21_2,
                   fit22_2, fit23_2, fit24_2, fit25_2, fit26_2, fit27_2,
                   fit28_2, fit31_2, fit32_2, fit33_2, fit34_2, fit35_2,
                   fit36_2, fit37_2, fit38_2, fit11_3, fit12_3, fit13_3,
                   fit14_3, fit15_3, fit16_3, fit17_3, fit18_3, fit21_3,
                   fit22_3, fit23_3, fit24_3, fit25_3, fit26_3, fit27_3,
                   fit28_3, fit31_3, fit32_3, fit33_3, fit34_3, fit35_3,
                   fit36_3, fit37_3, fit38_3)
    fitall <- Filter(Negate(is.null), fitall)
    fitres <- do.call(anova, fitall)
    fitres <- as.data.frame(fitres)
    i = 1
    fitres$Pvalue = 0
    fitres$formula = 0
    fitres$id = seq(1:length(fitall))
    for (i in 1:length(fitall)) {
      sum <- summary(fitall[[i]])
      coff <- as.data.frame(sum$coefficients)
      tryCatch({
        fitres$Pvalue[i] <- coff["rate", grepl("Pr",
                                               names(coff))]
      }, error = function(e) {
        fitres$Pvalue[i] <<- coff["rate", grepl("value",
                                                names(coff))]
      })
      fitres$formula[i] <- paste(link, "mix effect model:",
                                 as.character(sum[["call"]][["formula"]])[3])
    }
    fitres$model <- paste("model", seq(1:length(fitres$formula)))
    fitres$celltype <- selectedcelltype
    fitres_sorted <- fitres[order(fitres$AIC), ]
    fitres_sorted <- fitres_sorted[, c("id", "model", "celltype",
                                       "formula", "AIC", "Pvalue")]

    fitres_sorted_all <- rbind(fitresfixall,fitres_sorted)
    fitres_sorted_all <- fitres_sorted_all[order(fitres_sorted_all$AIC), ]

    print(paste0(length(fitres_sorted_all$AIC)," ",selectedcelltype, " fitting models have been calculated successful",
                 ""))
    fitres_sorted_all$id <- seq(1:length(fitres_sorted_all$id))
    finalresult <- list(fitall2, fitres_sorted_all, cwas.data)
    names(finalresult) <- c("Models", "Models_info", "Rawdata")
    return(finalresult)
  }, error = function(e) {
    cat("Error in iteration", i, ":", conditionMessage(e),
        "\n")
    return("NotSuccess")
  })
}
