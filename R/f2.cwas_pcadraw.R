#' Generate PCA Plot for CWAS Data
#'
#' This function generates a PCA (Principal Component Analysis) plot for the given CWAS data based on cell type rates and a specified grouping variable.
#'
#' @param cwas.data The CWAS data containing information on samples, cell types, rates, and the grouping variable.
#' @param div_group The name of the grouping variable used for color-coding in the PCA plot (default is "Phenotype").
#'
#' @return A PCA plot visualizing the CWAS data.
#'
#' @export
#'
#' @examples cwas_pcadraw(cwas.test.data)
cwas_pcadraw <- function(cwas.data,div_group="Phenotype"){

  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  myattr <- attributes(cwas.data)
  genelist <-  myattr$genelist
  cwas.data <- cwas.data[,!colnames(cwas.data) %in%genelist]

  df <- cwas.data[,c("Sample","Celltype","rate",div_group)]
  names(df) <- c("Sample","Celltype","rate","div_group")

  df$div_group <- ifelse(df$div_group=="Disease",myattr$Disease_label,myattr$Control_label)
  df <- dplyr::distinct(df)
  wide_df <- tidyr::spread(df, key = Celltype, value = rate)
  wide_df[is.na(wide_df)] <- 0
  data <- t(wide_df[, !(names(wide_df) %in% c("Sample", "div_group"))])
  p1 <- draw_pca(data,factor(wide_df$div_group))
  p1
  return(p1)
}
