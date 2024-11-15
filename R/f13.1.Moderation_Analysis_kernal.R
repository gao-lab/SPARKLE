#' Moderation Analysis Function
#'
#' This function performs moderation analysis on CWAS data. It supports different types of moderation including cell-cell, gene-cell, and gene-gene mediation effects.
#'
#' @param cwas.data A data frame containing CWAS data.
#' @param Best_model A list containing the best model chosen for the analysis.
#' @param X.cell A character string specifying the cell type to be used as the independent variable (IV). Default is NULL.
#' @param X.gene A character string specifying the gene to be used as the independent variable (IV). Default is NULL.
#' @param Moderation.cell A character string specifying the cell type to be used as the moderator. Default is NULL.
#' @param Moderation.gene A character string specifying the gene to be used as the moderator. Default is NULL.
#' @param nsim An integer specifying the number of simulations for bootstrapping. Default is 1000.
#' @param seed An integer specifying the seed for random number generation. Default is 1.
#' @param CI A character string specifying the confidence interval method. Default is "boot".
#'
#' @return A list containing the results of the moderation analysis.
#'
#' @export
#'
moderation_analysis_auto <- function(cwas.data,Best_model, X.cell=NULL,X.gene=NULL,Moderation.cell=NULL,Moderation.gene=NULL,nsim=1000, seed=1,CI = "boot"){


  formula_process_string_vector <- function(input_vector,ratename,clean=T) {

    # 初始化结果向量
    Fixeffect <- character(length(input_vector))
    Random_effect <- character(length(input_vector))
    T <- logical(length(input_vector))

    for (i in seq_along(input_vector)) {
      element <- input_vector[i]

      if (element == "rate") {
        # 如果是 "rate"，不做任何操作
        Fixeffect[i] <- element
        Random_effect[i] <- ""

      } else if (grepl("\\(", element) && grepl("\\)", element)) {
        # 如果包含 "()"，将该字符赋值给 Random_effect，并设置 T 为 TRUE
        Fixeffect[i] <- ""
        Random_effect[i] <- element

      } else {
        # 否则，赋值给 Fixeffect
        Fixeffect[i] <- element
        Random_effect[i] <- ""

      }
    }
    Fixeffect <- Fixeffect[Fixeffect != ""]
    Random_effect <- Random_effect[Random_effect != ""]

    Group_info <- gsub("1|rate|\\(|\\)|\\||\\|\\|", "", Random_effect)
    Group_info <- gsub(" ", "", Group_info)
    if (grepl("/", Group_info[1])) {
      Group_info <- strsplit(Group_info, "/")[[1]]
    } else {
      Group_info
    }

    Covariate <- gsub("rate|\\(|\\)|\\||\\|\\|", "", Fixeffect)
    Covariate <- Covariate[Covariate != ""]
    # 返回一个列表，包含 Fixeffect, Random_effect 和 T
    ratename<- gsub("-", "_",ratename)

    if(length(Random_effect)>1){
      Random_effect <- Random_effect[2]
    }
    Cleaned_random_effect <- gsub("rate", ratename, Random_effect)
        # if(clean){Cleaned_random_effect <- gsub("rate", ratename, Random_effect)
    # }else{Cleaned_random_effect <- Random_effect }
    return(list(Fixeffect = Fixeffect, Random_effect = Cleaned_random_effect,Group_info=Group_info,Covariate=Covariate))
  }

  replace_double_pipe <- function(vector) {
    # 使用gsub函数替换所有"||"为"|"
    vector <- gsub("\\|\\|", "|", vector)
    return(vector)
  }

  getfomula <- function(X.cell,Cleaninfo=T){
    names(Best_model[["Chosen_model"]])  <- names(Best_model[["All_models"]])
    names(Best_model[["Chosen_model"]]) <- gsub("-", "_", names(Best_model[["Chosen_model"]]) )
    X.cell.fomula <- as.character(Best_model[["Chosen_model"]][[X.cell]]@call[["formula"]])[3]
    X.cell.fomula <- strsplit(X.cell.fomula, " \\+ ")[[1]]
    X.cell.fomula <-replace_double_pipe(X.cell.fomula )
    X.cell.fomula.result <- formula_process_string_vector(X.cell.fomula,ratename = X.cell,clean = Cleaninfo)

    return(X.cell.fomula.result)

  }



  cwas.data$Celltype <- gsub("-", "_", cwas.data$Celltype)
  Moderation.cell<- gsub("-", "_", Moderation.cell)

  celltypename <- unique(cwas.data$Celltype)
  attr.info <- attributes(cwas.data)
  genelist <-  attr.info$genelist
  X.cell.fomula.result=NULL
  Moderation.fomula.result=NULL

  if(!is.null(X.cell)){

    X.cell<- gsub("-", "_", X.cell)

    if(is.null(Moderation.gene)){# 细胞自变量，细胞调节变量
      Moderation <- Moderation.cell
      cwas.data <- cwas.data[,!colnames(cwas.data) %in%genelist]
      cwas.data <- cwas.data[cwas.data$Celltype %in%c(X.cell,Moderation), ]
      X.cell.fomula.result <- getfomula(X.cell)
      Moderation.fomula.result <- getfomula(Moderation)
      wide_df <- tidyr::spread(cwas.data, key = Celltype, value = rate)
      wide_df[is.na(wide_df)] <- 0
      X.info <-X.cell

    }else{# 细胞自变量，同种或者不同种细胞调节变量

      if(is.null(Moderation.cell)){Moderation.cell=X.cell}

      cwas.data.X <- cwas.data[cwas.data$Celltype %in%c(X.cell), ]
      cwas.data.X <- cwas.data.X[,!colnames(cwas.data.X) %in%genelist]
      cwas.data.Moderation <- cwas.data[cwas.data$Celltype %in%c(Moderation.cell), ]
      cwas.data.Moderation <- cwas.data.Moderation[,colnames(cwas.data.Moderation) %in%c("Sample","Celltype","rate",Moderation.gene)]
      # cwas.data.Moderation$Moderation=cwas.data.Moderation$rate*cwas.data.Moderation[,Moderation.gene]#基因总量

      cwas.data.Moderation$Moderation= cwas.data.Moderation[,Moderation.gene]#基因均值

      cwas.data.Moderation <- cwas.data.Moderation[,c("Sample","Moderation")]
      X.cell.fomula.result <- getfomula(X.cell)
      Moderation.fomula.result$Random_effect <- ""
      wide_df <-merge(cwas.data.X,cwas.data.Moderation,by="Sample")
      wide_df[is.na(wide_df)] <- 0
      colnames(wide_df) <- gsub("rate", X.cell, colnames(wide_df))
      Moderation.name <- paste0(Moderation.gene,"_",Moderation.cell)
      colnames(wide_df) <- gsub("Moderation", Moderation.name, colnames(wide_df))
      X.info <- X.cell
      Moderation <- Moderation.name

    }

  }else if(is.null(X.cell)&!is.null(X.gene)){

    if(is.null(Moderation.cell)){stop("Please Input Moderation Cell label !",
                                      call. = FALSE) }else{
                                        # 基因自变量，细胞调节变量
                                        diff_set <- genelist[!genelist %in% X.gene]
                                        cwas.data <- cwas.data[,!colnames(cwas.data) %in%diff_set]
                                        cwas.data <- cwas.data[cwas.data$Celltype %in%Moderation.cell, ]
                                        Moderation.fomula.result <- getfomula(Moderation.cell)
                                        X.cell.fomula.result$Random_effect <- ""
                                        wide_df <- cwas.data
                                        wide_df[is.na(wide_df)] <- 0
                                        colnames(wide_df) <- gsub("rate", Moderation.cell, colnames(wide_df))
                                        X.info <-X.gene
                                        Moderation <- Moderation.cell

                                      }

  }else{
    stop("Please Input X.cell or X.gene label !",     call. = FALSE)

  }

  BM <-bruceR::PROCESS(wide_df, y="Phenotype", x=X.info, mods =Moderation,ci = CI,clusters=X.cell.fomula.result$Group_info,hlm.re.y=X.cell.fomula.result$Random_effect,hlm.re.m =Moderation.fomula.result$Random_effect ,covs = X.cell.fomula.result$Covariate, nsim=nsim, seed=seed)

  print( BM$results[[1]]$jn[[1]]$plot)


  return(BM)
}
