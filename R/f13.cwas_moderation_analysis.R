#' cwas_moderation_analysis
#'
#' This function performs moderation analysis between cells or genes. It calculates the moderation effects between cells, genes, or a combination of both based on the specified method.
#'
#' @param cwas.data A data frame containing the input data for the analysis.
#' @param Best_model A model object containing information about the best model.
#' @param X.cell A string specifying the cell type for the moderation effect analysis. Default is NULL.
#' @param X.gene A string specifying the gene for the moderation effect analysis. Default is NULL.
#' @param Moderation.cell A string specifying the moderating cell type. Default is NULL.
#' @param Moderation.gene A string specifying the moderating gene. Default is NULL.
#' @param method A string specifying the method to use for the analysis. Acceptable values are "cell_moderation_cell", "gene_moderation_cell", or "cell_moderation_gene".
#' @param num_cores An integer specifying the number of cores to use for parallel computation. If NULL, sequential computation is performed. Default is NULL.
#'
#' @return An object containing the results of the moderation analysis.
#'
#' @details
#' - If the `method` parameter is "cell_moderation_cell", the function performs cell-to-cell moderation effect analysis.
#' - If the `method` parameter is "gene_moderation_cell", the function performs gene-to-cell moderation effect analysis. Depending on the `num_cores` parameter, it decides whether to perform sequential or parallel computation.
#' - If the `method` parameter is "cell_moderation_gene", the function performs cell-to-gene moderation effect analysis. Depending on the `num_cores` parameter, it decides whether to perform sequential or parallel computation.
#'
#' @export
cwas_moderation_analysis <- function(cwas.data,Best_model, X.cell=NULL,X.gene=NULL,Moderation.cell=NULL,Moderation.gene=NULL,method=NULL,num_cores =NULL ){


  if (is.null(method)) {
    stop("method cannot be NULL,please input [cell_moderation_cell] or [gene_moderation_cell] or[cell_moderation_gene]", call. = FALSE)
  }

  if(method=="cell_moderation_cell"){ ##### F3.1 细胞-细胞 调节效应计算  #####

    final <- cwas_cell_moderation_cell_analysis(cwas.data,Best_model,X.cell,selected_celltype=Moderation.cell)

  }else if(method=="gene_moderation_cell"){
    if(is.null(num_cores)){#### F3.2.1 基因-细胞 调节效应（顺序计算） #####
      final <-cwas_gene_moderation_cell_analysis(cwas.data,Best_model,X.cell,Moderation.cell,selected_gene=Moderation.gene)
    }else{#### F3.2.2 基因-细胞 调节效应（并行计算） #####

      final <-cwas_gene_moderation_cell_analysis_para(cwas.data,Best_model,X.cell,Moderation.cell,selected_gene=Moderation.gene, num_cores = num_cores)
    }

  }else if(method=="cell_moderation_gene"){

    if(is.null(num_cores)){#### F3.3.1 细胞-基因 调节效应（顺序计算） #####
      final <- cwas_cell_moderation_gene_analysis(cwas.data,Best_model,X.gene,selected_celltype=Moderation.cell)

    }else{#### F3.3.2 细胞-基因 调节效应（并行计算） #####
      final <- cwas_cell_moderation_gene_analysis_para(cwas.data, Best_model, X.gene, selected_celltype = Moderation.cell, num_cores = num_cores)
    }


  }else{

    cat("method cannot be NULL,please input [cell_moderation_cell] or [gene_moderation_cell] or[cell_moderation_gene]")

  }
  sf.rds(final,file = paste0("013_Moderation_analysis_",method))

  return(final)


}
