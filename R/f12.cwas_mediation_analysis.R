#' Conduct CWAS Mediation Analysis
#'
#' This function conducts CWAS (Cell-type Weighted Analysis of Single-cell data) mediation analysis based on the specified method.
#'
#' @param cwas.data A data frame containing the CWAS results.
#' @param Best_model The best model selected for mediation analysis.
#' @param X.cell Vector specifying the cell types used as predictors in the mediation analysis. Default is NULL.
#' @param X.gene Vector specifying the genes used as predictors in the mediation analysis. Default is NULL.
#' @param mediate.cell Vector specifying the cell types to be used as mediators in the mediation analysis.
#' @param mediate.gene Vector specifying the genes to be used as mediators in the mediation analysis.
#' @param method The method for conducting mediation analysis. It should be one of "cell_mediate_cell", "gene_mediate_cell", or "cell_mediate_gene".
#' @param num_cores Number of CPU cores to be used for parallel computation. Default is NULL.
#'
#' @return The results of CWAS mediation analysis.
#'
#' @export
#'
cwas_mediation_analysis <- function(cwas.data,Best_model, X.cell=NULL,X.gene=NULL,mediate.cell=NULL,mediate.gene=NULL,method=NULL,num_cores =NULL, fix_effect = NULL,
                                    random_effect = NULL ){


  if (is.null(method)) {
    stop("method cannot be NULL,please input [cell_mediate_cell] or [gene_mediate_cell] or[cell_mediate_gene]", call. = FALSE)
  }



  if(method=="cell_mediate_cell"){ ##### F3.1 细胞-细胞 调节效应计算  #####

      final <- cwas_cell_mediate_cell_analysis(cwas.data,Best_model,X.cell,selected_celltype=mediate.cell,fix_effect = fix_effect,random_effect = random_effect)

  }else if(method=="gene_mediate_cell"){
    if(is.null(num_cores)){#### F3.2.1 基因-细胞 调节效应（顺序计算） #####
      final <-cwas_gene_mediate_cell_analysis_singlecore(cwas.data,Best_model,X.cell,Mediation = mediate.gene, fix_effect = fix_effect,random_effect = random_effect)
    }else{#### F3.2.2 基因-细胞 调节效应（并行计算） #####
      final <-cwas_gene_mediate_cell_analysis_para(cwas.data,Best_model,X.cell,  Mediation=mediate.gene, fix_effect = fix_effect,random_effect = random_effect, num_cores = num_cores)
    }

  }else if(method=="cell_mediate_gene"){

    if(is.null(num_cores)){#### F3.3.1 细胞-基因 调节效应（顺序计算） #####
      final <- cwas_cell_mediate_gene_analysis_singlecore(cwas.data,Best_model,X.gene,selected_celltype=mediate.cell, fix_effect = fix_effect,random_effect = random_effect)

    }else{#### F3.3.2 细胞-基因 调节效应（并行计算） #####

      final <- cwas_cell_mediate_gene_analysis_para(cwas.data, Best_model, X.gene, selected_celltype = mediate.cell, num_cores = num_cores, fix_effect = fix_effect,random_effect = random_effect)

    }


  }else{

    cat("method cannot be NULL,please input [cell_mediate_cell] or [gene_mediate_cell] or[cell_mediate_gene]")

  }

  sf.rds(final,file = paste0("012_Mediation_analysis_",method))

  return(final)


}
