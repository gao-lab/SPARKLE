
OR_plot_compare <- function (data, show_variable = c("Celltype", "Model", "OR(95%CI)",
                                             "SPARKLE (p-value)", "SPARKLE (p-value) adjusted","Wilcoxon (p-value)"), plot.position = 4, Tablename = "Tablename")
{

  Wilcoxon_function <- function (cwas.data)
  {


    mydata <- cwas.data
    cellnameall <- unique(mydata$Celltype)
    plotlist <- list()
    result.all <- c()
    for (i in 1:length(cellnameall)) {
      cellname <- cellnameall[i]
      df <- mydata
      df <- subset(df, Celltype == cellname)
      group1 <- df[df$Phenotype =="Control", "rate"]
      group2 <- df[df$Phenotype =="Disease", "rate"]
      test_result <- wilcox.test(group1, group2)
      test_result$p.value
      result.df <- data.frame(cellname, test_result$p.value)
      result.all <- rbind(result.all, result.df)
    }
    result.all$p.value_Non_adjusted = result.all$test_result.p.value
    result.all$adjusted.p.vaule_Bonferroni <- p.adjust(result.all$test_result.p.value,
                                                       method = "bonferroni")
    result.all$adjusted.p.vaule_Benjamini_Hochberg <- p.adjust(result.all$test_result.p.value,
                                                               method = "BH")
    result.all$adjusted.p.vaule_Holm <- p.adjust(result.all$test_result.p.value,
                                                 method = "holm")
    result.all$adjusted.p.vaule_Benjamini_Yekutieli <- p.adjust(result.all$test_result.p.value,
                                                                method = "BY")
    result.all$adjusted.p.vaule_Hommel <- p.adjust(result.all$test_result.p.value,
                                                   method = "hommel")
    result.all$adjusted.p.vaule_Hochberg <- p.adjust(result.all$test_result.p.value,
                                                     method = "hochberg")
    result.all$test_result.p.value = NULL
    return(list(result.all, plotlist))
  }


  raw.sparkle.data  <- data[["All_models"]][[1]][["inputdata"]]

  cwas.data <- data[["All_models"]][[1]][["Rawdata"]]




  atri_info <- attributes(cwas.data)
  data <- data[[2]]


  wilicox2 <- Wilcoxon_function (raw.sparkle.data)


  wilico <- wilicox2[[1]]
  wilico$celltype <- wilico$cellname
  wilico$Wilcoxon <- wilico$p.value_Non_adjusted
  wilico$sig <- ifelse(wilico$Wilcoxon<0.001,"***",
                       ifelse(wilico$Wilcoxon<0.01,"**",
                              ifelse(wilico$Wilcoxon<0.05,"*","")))

  wilico$Wilcoxon <- ifelse(wilico$Wilcoxon<0.01,"<0.01",round(wilico$Wilcoxon,2))
  wilico$Wilcoxon <- paste0(wilico$Wilcoxon,wilico$sig)
  wilico <- wilico[,c("celltype","Wilcoxon")]
  wilico$"Wilcoxon (p-value)" <- wilico$Wilcoxon


  data <- merge(data,wilico,by="celltype")

  data$formula <- gsub("binomial link function:","",data$formula)

  show_variable2 <- show_variable[1:plot.position - 1]
  show_variable2[plot.position] <- c("                        ")
  num <- length(show_variable)
  num1 <- plot.position + 1
  num2 <- length(show_variable) + 1
  show_variable2[num1:num2] <- show_variable[plot.position:num]
  numeric_cols <- sapply(data, is.numeric)
  data1 <- data
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) round(x,
                                                                     2))
  data1[numeric_cols] <- lapply(data[numeric_cols], function(x) format(x,
                                                                       scientific = TRUE, digits = 2))
  data$"                        " <- paste(rep(" ", length(data[, 1])), collapse = " ")
  data$"  " <- paste(rep(" ", length(data[, 1])), collapse = " ")
  data$"OR(95%CI)" <- paste0(data1$OR, " ", "[", data1$CI_lower,
                             ",", data1$CI_upper, "]")
  data$OR <- ifelse(data$OR > 20, 10, data$OR)
  data$"OR(95%CI)"[grepl("Inf", data1$OR)] <- "Not reliable"
  data$CI_upper <- ifelse(data$CI_upper > 10, 10, data$CI_upper)
  data$"SPARKLE (p-value)" <- ifelse(data$Pvalue < 0.01, "<0.01", data$Pvalue)
  data$"SPARKLE (p-value)" <- paste(data$"SPARKLE (p-value)", data$significance)
  data$"SPARKLE (p-value)" <- paste0("  ", data$"SPARKLE (p-value)")
  data$"SPARKLE (p-value)"[grepl("Inf", data1$OR)] <- "Not reliable"
  tryCatch({
    data$"SPARKLE (p-value) adjusted" <- ifelse(data$adjustedPvalue <
                                                  0.01, "<0.01", data$adjustedPvalue)
    data$"SPARKLE (p-value) adjusted" <- paste(data$"SPARKLE (p-value) adjusted",
                                               data$significance.adjusted)
    data$"SPARKLE (p-value) adjusted" <- paste0("  ", data$"SPARKLE (p-value) adjusted")
    data$"SPARKLE (p-value) adjusted"[grepl("Inf", data1$OR)] <- "Not reliable"
  }, error = function(e) {
    print(paste("Warning Message:", "No Adjusted p value in this data"))
    show_variable2 <<- show_variable2[show_variable2 !=
                                        "SPARKLE (p-value) adjusted"]
  })
  data$Pvalue[grepl("Inf", data1$OR)] <- "Not reliable"
  if (is.null(data$formula)) {
    data$formula <- c("Integrated model")
  }
  data$Model <- data$formula
  data$Celltype <- data$celltype
  data$Significance <- data$significance
  tm <- forestploter::forest_theme(base_size = 10, ci_pch = 15,
                                   ci_col = "blue4", ci_fill = "blue4", ci_alpha = 0.8,
                                   ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2, refline_lwd = 2,
                                   refline_lty = "dashed", refline_col = "red4", vertline_lwd = 1,
                                   vertline_lty = "dashed", vertline_col = "grey20", footnote_cex = 0.6,
                                   footnote_fontface = "italic", footnote_col = "red4")
  footnote1 <- c("* p<0.05  ** p<0.01  *** p<0.001 ")
  info_parts <- list(paste("Phenotype:", ifelse(is.null(atri_info$Phenotype),
                                                "Not Provided", atri_info$Phenotype), "  ", "Control_label:",
                           ifelse(is.null(atri_info$Control_label), "Not Provided",
                                  atri_info$Control_label), "  ", "Disease_label:",
                           ifelse(is.null(atri_info$Disease_label), "Not Provided",
                                  atri_info$Disease_label)), paste("Group:", ifelse(is.null(atri_info$Group),
                                                                                    "Not Provided", atri_info$Group), "  ", "Subgroup:",
                                                                   ifelse(is.null(atri_info$Subgroup), "Not Provided",
                                                                          atri_info$Subgroup), "  ", "Covariate 1:", ifelse(is.null(atri_info$Covariate1),
                                                                                                                            "Not Provided", atri_info$Covariate1), "  ", "Covariate 2:",
                                                                   ifelse(is.null(atri_info$Covariate2), "Not Provided",
                                                                          atri_info$Covariate2)))
  info <- paste(info_parts, collapse = "\n")
  footnote2 <- paste0(footnote1, info)
  pp <- forestploter::forest(data[, show_variable2], est = data$OR,
                             lower = data$CI_lower, upper = data$CI_upper, sizes = 0.8,
                             ci_column = 4, ref_line = 1, arrow_lab = c("Low risk",
                                                                        "High Risk"), xlim = c(-1, 10), ticks_at = c(0,
                                                                                                                     1, 3, 5), theme = tm, footnote = footnote2)
  pp <- forestploter::edit_plot(pp, row = which(data$Pvalue <
                                                  0.0500001), gp = grid::gpar(col = "red4", fontface = "italic"))
  pp <- forestploter::insert_text(pp, text = Tablename, col = 1:5,
                                  part = "header", gp = grid::gpar(fontface = "bold"))
  pp <- forestploter::add_border(pp, part = "header", where = "bottom")
  print(pp)
  return(pp)
}


