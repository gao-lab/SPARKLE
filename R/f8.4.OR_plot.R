#' OR Plot Function
#'
#' This function generates a forest plot to visualize odds ratios (OR) and confidence intervals (CI) from provided data.
#'
#' @param data A list containing required data for generating the forest plot.
#' @param show_variable A character vector specifying the variables to display in the forest plot.
#'   Default is c("Celltype", "Model", "OR(95%CI)", "P-value", "P-value adjusted").
#' @param plot.position An integer specifying the position to insert an empty space in the `show_variable` vector.
#'   Default is 4, which corresponds to inserting an empty space after "P-value".
#' @param Tablename A character string specifying the table name to display in the forest plot header.
#'   Default is "Tablename".
#'
#' @return A forest plot visualizing odds ratios and confidence intervals based on the provided data.
#'
#' @examples
#' # Example usage:
#' OR_plot(data = cwas.test.data.single.model.best)
#'
#' @keywords plot
#' @export

OR_plot <- function(data,show_variable=c("Celltype","Model","OR(95%CI)","P-value","P-value adjusted"),plot.position=4,Tablename="Tablename"){

  raw.sparkle.data  <- data[["All_models"]][[1]][["inputdata"]]

  if(is.null(raw.sparkle.data)){

    pp <- OR_plot_origin(data = data,show_variable =show_variable,plot.position =plot.position,Tablename =Tablename  )

}else{

    pp <- OR_plot_compare(data = data,show_variable =c("Celltype", "Model", "OR(95%CI)",
                                                "SPARKLE (p-value)", "SPARKLE (p-value) adjusted","Wilcoxon (p-value)"),
                   plot.position =plot.position,Tablename =Tablename  )

  }
  return(pp)

}
