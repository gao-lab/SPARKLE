#' Calculate All Models with Two Cell Types as Variables
#'
#' This function computes all possible models involving combinations of two specified cell types from the CWAS data.
#'
#' @param cwas.data The CWAS dataset containing the rates and phenotypes of cell types.
#' @param selected_celltype A character vector specifying the cell types to consider for modeling (default is NULL, which uses all unique cell types in the dataset).
#' @param linkfunction The link function to use for modeling (default is "binomial").
#' @param interaction Logical indicating whether to consider interaction effects between cell types (default is FALSE).
#'
#' @return A list of calculated models for each combination of two cell types.
#'
#' @export
#'

cwas_2celltype_allmodel_cal <-  function(cwas.data, selected_celltype= NULL,linkfunction="binomial",interaction=F){

  #' @examples cwas_2celltype_allmodel_cal(cwas.test.data)

  attr.info <- attributes(cwas.data)
  genelist <-  attr.info$genelist
  cwas.data <- cwas.data[,!colnames(cwas.data) %in%genelist]
  ### 生成所有的两种细胞的组合情况
  #library(gtools)
  if(is.null(selected_celltype)){ selected_celltype <- unique(cwas.data$Celltype)}
  selected_celltype <- as.character(selected_celltype)
  cellcom <- gtools::permutations(n = length(selected_celltype), r = 2, v = selected_celltype)
  #cellcom <- t(combn(selected_celltype,2))
  cellcom <- as.matrix(cellcom)
  colnames(cellcom) <- c("celltype1","celltype2")
  mdlist <- list()
  comparelist <- list()

  ### 数据长转宽表格
  cwas.data$num=NULL
  cwas.data$allnum=NULL
  wide_data <- tidyr::pivot_wider(data = cwas.data,
                           id_cols =setdiff(names(cwas.data), c("Celltype", "rate")),
                           names_from = Celltype,
                           values_from = rate,
                           values_fill = 0)


  dfnum <- length(wide_data[1,])
  celltypenum <- dfnum-length(selected_celltype)+1

  ### 计算每一种细胞类型组合的时候的结果

  for(i in 1:length(cellcom[,1])){
    #cell2df <- wide_data %>% dplyr::select(colnames(wide_data)[-c(celltypenum:dfnum)],cellcom[i,1],cellcom[i,2])
    # 获取需要保留的列名
    selected_cols <- colnames(wide_data)[-c(celltypenum:dfnum)]

    # 从 wide_data 中选择指定的列以及 cellcom[i,1] 和 cellcom[i,2] 列
    cell2df <- wide_data[, c(selected_cols, cellcom[i, 1], cellcom[i, 2])]

    cell2df <- as.data.frame(cell2df)
    num_cols <- ncol(cell2df)
    colnames(cell2df)[num_cols] <- "rate2"
    colnames(cell2df)[num_cols-1] <- "rate"



    if(cell2df$Phenotype[1]%in%c("Disease",	"Control")){cell2df$Phenotype <- ifelse(cell2df$Phenotype=="Control",0,1)}


    if(interaction==T){
      print("Cell type interaction in all models calculation")

      tryCatch({
        mdlist[[i]] <- interaction_2model_cal_v5(cell2df,selectedcelltype = paste(cellcom[i,1],cellcom[i,2]),link = linkfunction)
      }, error = function(e) {
        # 如果出现错误，提示无法完成混合模型计算,执行固定效应模型
        cat("Error: Mixed effect model was not applicable for this data.\n")
      })


    }else{
      print("2 Cell type all models calculation")


      tryCatch({
        comparelist[[i]] <- compare_2model_cal_v5(cell2df,selectedcelltype =   paste(cellcom[i,1],cellcom[i,2]),link = linkfunction)

      }, error = function(e) {
        # 如果出现错误，提示无法完成混合模型计算,执行固定效应模型
        cat("Error: Mixed effect model was not applicable for this data.\n")
      })




    }


    print(i)

  }



  if(interaction==T){
    names(mdlist) <-  paste(cellcom[,1],cellcom[,2])
    mdlist <- mdlist[!sapply(mdlist, is.null)]
    sf.rds(mdlist,file = "04CWAS_2_cell_interaction_model")
    return(mdlist)
  }else{
    names(comparelist) <-  paste(cellcom[,1],cellcom[,2])
    comparelist <- comparelist[!sapply(comparelist, is.null)]
    sf.rds(comparelist,file = "03CWAS_2_cell_all_model")
    return(comparelist)
  }



}

