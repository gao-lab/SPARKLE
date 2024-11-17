#' Fill Missing Celltypes with Zero and Retain Sample Metadata
#'
#' This function fills missing cell type rates in the provided dataset by adding any missing
#' combinations of samples and cell types, setting their rate to 0. Additionally, it retains
#' and matches all other metadata related to the samples.
#'
#' @param data A dataframe that contains at least sample, celltype, and rate columns.
#' @param sample_col A string specifying the column name for sample identifiers. Default is "Sample".
#' @param celltype_col A string specifying the column name for cell types. Default is "Celltype".
#' @param rate_col A string specifying the column name for the rate values. Default is "rate".
#'
#' @return A dataframe where missing sample-celltype combinations have been filled in with
#'         a rate of 0, and all metadata related to the samples has been retained and matched.
#' @export
cwas_fill_missing_celltype <- function(data, sample_col = "Sample", celltype_col = "Celltype") {
  # 提取与Sample相关的元数据列（非数值列）

  # 识别所有数值型列
  numeric_cols <- names(data)[sapply(data, is.numeric)]

  meta_cols <- setdiff(names(data), c(celltype_col,numeric_cols))


  # 提取数值列和元数据列，确保元数据唯一
  meta_data <- dplyr::distinct(data[, meta_cols])


  # 找到所有的 sample 和 celltype 组合
  complete_samples <- expand.grid(Sample = unique(data[[sample_col]]),
                                  Celltype = unique(data[[celltype_col]]))

  # 合并实际数据，对所有数值型列进行补全，缺失值补0
  df_filled <- merge(complete_samples, data[, c(sample_col, celltype_col, numeric_cols)],
                     by = c(sample_col, celltype_col), all.x = TRUE)

  # 对所有数值型列补0
  df_filled[numeric_cols] <- lapply(df_filled[numeric_cols], function(x) {
    x[is.na(x)] <- 0
    return(x)
  })

  # 确保没有重复的行
  df_filled <- unique(df_filled)
  meta_data_unique <- unique(meta_data)

  # 根据Sample匹配元数据信息，保留所有列
  df_filled <- merge(meta_data_unique, df_filled, by = "Sample")


  # 返回补全后的数据
  myattr <-attributes( df_filled)
  myattrinfo <- attributes(data)

  myattr$Sample <- df_filled$Sample
  myattr$Phenotype <-  df_filled$Phenotype
  myattr$Group <- df_filled$Group
  myattr$Subgroup <- df_filled$Subgroup
  myattr$Cov1 <- df_filled$Covariate1
  myattr$Cov2 <- df_filled$Covariate2
  myattr$Control_label <- myattrinfo$Control_label
  myattr$Disease_label <- myattrinfo$Disease_label
  myattr$Cellrate <- df_filled$Cellrate

#  myattr <<- myattr
  attributes(df_filled) <- myattr

  return(df_filled)
}

