or_model <- function (data, show_variable = c("Celltype", "Covariant", "OR(95%CI)",
                                              "P-value", "P-value adjusted"), plot.position = 4, Tablename = "Tablename")
{


  cwas.data <- data
  atri_info <- attributes(cwas.data)

  show_variable2 <- show_variable[1:plot.position - 1]
  show_variable2[plot.position] <- c(" ")
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
  data$" " <- paste(rep(" ", length(data[, 1])), collapse = " ")
  data$"  " <- paste(rep(" ", length(data[, 1])), collapse = " ")
  data$"OR(95%CI)" <- paste0(data1$OR, " ", "[", data1$CI_lower,
                             ",", data1$CI_upper, "]")
  data$OR <- ifelse(data$OR > 20, 10, data$OR)
  data$"OR(95%CI)"[grepl("Inf", data1$OR)] <- "Not reliable"
  data$CI_upper <- ifelse(data$CI_upper > 10, 10, data$CI_upper)
  data$"P-value" <- ifelse(data$Pvalue < 0.01, "<0.01", data$Pvalue)
  data$"P-value" <- paste(data$"P-value", data$significance)
  data$"P-value" <- paste0("  ", data$"P-value")
  data$"P-value"[grepl("Inf", data1$OR)] <- "Not reliable"
  tryCatch({
    data$"P-value adjusted" <- ifelse(data$adjustedPvalue <
                                        0.01, "<0.01", data$adjustedPvalue)
    data$"P-value adjusted" <- paste(data$"P-value adjusted",
                                     data$significance.adjusted)
    data$"P-value adjusted" <- paste0("  ", data$"P-value adjusted")
    data$"P-value adjusted"[grepl("Inf", data1$OR)] <- "Not reliable"
  }, error = function(e) {
    print(paste("Warning Message:", "No Adjusted p value in this data"))
    show_variable2 <<- show_variable2[show_variable2 !=
                                        "P-value adjusted"]
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
  #footnote2 <- paste0(footnote1, info)
  footnote2 <- footnote1
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
