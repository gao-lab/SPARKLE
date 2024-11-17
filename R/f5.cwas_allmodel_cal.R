#' Calculate Results for Models of a Single Cell Type
#'
#' This function calculates results for all models containing a specific cell type in the CWAS data.
#'
#' @param cwas.data The CWAS data containing model information.
#' @param link_function The link function to use for model fitting (default is "binomial").
#'
#' @return A list of model results for each cell type.
#'
#' @export
#'

cwas_allmodel_cal <- function(cwas.data,link_function="binomial",method="mix_effect",num_cores=NULL){

  if(is.null(num_cores)){
    mdlist <- cwas_allmodel_cal_singlecore(cwas.data=cwas.data,link_function=link_function,method=method)
  }else{
    mdlist <- cwas_allmodel_cal_multicore(cwas.data=cwas.data, link_function = link_function, method =method, num_cores=num_cores)

  }


  return(mdlist)
}
