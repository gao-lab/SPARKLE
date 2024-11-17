#' Mediation Analysis Auto
#'
#' This function performs mediation analysis using cell types or genes as independent and/or mediating variables. It processes the input data and leverages the `bruceR::PROCESS` function to conduct the analysis.
#'
#' @param cwas.data Data frame containing cell types, gene expression, and phenotype data.
#' @param Best_model List containing model selection information.
#' @param X.cell Character. Name of the independent cell type variable. Default is NULL.
#' @param X.gene Character. Name of the independent gene variable. Default is NULL.
#' @param Moderation.cell Character. Name of the mediating cell type variable. Default is NULL.
#' @param Moderation.gene Character. Name of the mediating gene variable. Default is NULL.
#' @param nsim Integer. Number of simulations. Default is 100.
#' @param seed Integer. Random seed for reproducibility. Default is 1.
#' @param CI Character. Method for confidence interval estimation. Default is "boot".
#' @param fix_effect Vector. Fixed effect covariates. Default is NULL.
#' @param random_effect Vector. Random effect covariates. Default is NULL.
#'
#' @return A list containing the results of the mediation analysis from the `bruceR::PROCESS` function.
#'
#' @examples
#' result <- mediation_analysis_auto(
#'   cwas.data = cwas.test.data,
#'   Best_model = cwas.test.data.single.model.best,
#'   X.cell = "MTF",
#'   Moderation.cell = "PIF",
#'   Moderation.gene = "TGFB1",
#'   nsim = 500,
#'   seed = 42,
#'   CI = "boot"
#' )
#'
#' @details
#' The function follows these main steps:
#' 1. `formula_process_string_vector`: Helper function to process formula strings.
#' 2. `getfomula`: Helper function to extract formula information from `Best_model`.
#' 3. Data preprocessing: Cleans and formats the data based on the specified independent and mediating variables.
#' 4. Mediation analysis: Uses `bruceR::PROCESS` to perform mediation analysis and returns the results.
#'
#' If `X.cell` is specified, the function will treat the cell type as the independent variable. If `X.gene` is specified, it will treat the gene as the independent variable. The mediating variable can be either a cell type or a gene, specified by `Moderation.cell` and `Moderation.gene` respectively.
#'
#' @note
#' Either `X.cell` or `X.gene` must be provided. If `X.gene` is used, `Moderation.cell` must be specified. The data frame `cwas.data` must contain a `Sample` column for merging operations.
#'
#' @export
#'
mediation_analysis_auto <- function(cwas.data,Best_model, X.cell=NULL,X.gene=NULL,Moderation.cell=NULL,Moderation.gene=NULL,nsim=100, seed=1,CI = "boot",fix_effect=NULL,random_effect=NULL){

  ###### F0 预备函数 ########
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
    if (grepl("/", Group_info)) {
      Group_info <- strsplit(Group_info, "/")[[1]]
    } else {
      Group_info
    }

    Covariate <- gsub("rate|\\(|\\)|\\||\\|\\|", "", Fixeffect)
    Covariate <- Covariate[Covariate != ""]
    # 返回一个列表，包含 Fixeffect, Random_effect 和 T
    ratename<- gsub("-", "_",ratename)
    if(clean){Cleaned_random_effect <- gsub("rate", ratename, Random_effect)
    }else{Cleaned_random_effect <- Random_effect }
    return(list(Fixeffect = Fixeffect, Random_effect = Cleaned_random_effect,Group_info=Group_info,Covariate=Covariate))
  }


  getfomula <- function(X.cell,Cleaninfo=T){
    names(Best_model[["Chosen_model"]])  <- names(Best_model[["All_models"]])
    names(Best_model[["Chosen_model"]]) <- gsub("-", "_", names(Best_model[["Chosen_model"]]) )
    X.cell.fomula <- as.character(Best_model[["Chosen_model"]][[X.cell]]@call[["formula"]])[3]
    X.cell.fomula <- strsplit(X.cell.fomula, " \\+ ")[[1]]
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

    if(is.null(Moderation.gene)){# 细胞自变量，细胞中介变量
      Moderation <- Moderation.cell
      cwas.data <- cwas.data[,!colnames(cwas.data) %in%genelist]
      cwas.data <- cwas.data[cwas.data$Celltype %in%c(X.cell,Moderation), ]
      #  X.cell.fomula.result <- getfomula(X.cell)
      # Moderation.fomula.result <- getfomula(Moderation)
      wide_df <- tidyr::spread(cwas.data, key = Celltype, value = rate)
      wide_df[is.na(wide_df)] <- 0
      X.info <-X.cell

    }else{# 细胞自变量，同种或者不同种细胞的基因作为中介变量

      if(is.null(Moderation.cell)){Moderation.cell=X.cell}

      cwas.data.X <- cwas.data[cwas.data$Celltype %in%c(X.cell), ]
      cwas.data.X <- cwas.data.X[,!colnames(cwas.data.X) %in%genelist]
      cwas.data.Moderation <- cwas.data[cwas.data$Celltype %in%c(Moderation.cell), ]
      cwas.data.Moderation <- cwas.data.Moderation[,colnames(cwas.data.Moderation) %in%c("Sample","Celltype","rate",Moderation.gene)]
      # cwas.data.Moderation$Moderation=cwas.data.Moderation$rate*cwas.data.Moderation[,Moderation.gene]#基因总量
      cwas.data.Moderation$Moderation= cwas.data.Moderation[,Moderation.gene]#基因均值

      cwas.data.Moderation <- cwas.data.Moderation[,c("Sample","Moderation")]
      #X.cell.fomula.result <- getfomula(X.cell)
      Moderation.fomula.result$Random_effect <- ""
      wide_df <-merge(cwas.data.X,cwas.data.Moderation,by="Sample")
      wide_df[is.na(wide_df)] <- 0
      colnames(wide_df) <- gsub("rate", X.cell, colnames(wide_df))
      Moderation.name <- paste0(Moderation.gene,"_",Moderation.cell)
      colnames(wide_df) <- gsub("Moderation", Moderation.name, colnames(wide_df))
      X.info <- X.cell
      Moderation <- Moderation.name

    }

  }else if(is.null(X.cell)&!is.null(X.gene)){# 基因自变量，细胞为中介变量

    if(is.null(Moderation.cell)){stop("Please Input Moderation Cell label !",
                                      call. = FALSE) }else{
                                        # 基因自变量，细胞调节变量
                                        diff_set <- genelist[!genelist %in% X.gene]
                                        cwas.data <- cwas.data[,!colnames(cwas.data) %in%diff_set]
                                        cwas.data <- cwas.data[cwas.data$Celltype %in%Moderation.cell, ]
                                        #Moderation.fomula.result <- getfomula(Moderation.cell)
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

  BM <-bruceR::PROCESS(wide_df, y="Phenotype", x=X.info, meds = Moderation,ci = CI, nsim=nsim, seed=seed,clusters =random_effect,covs = fix_effect  )

  print( BM$results[[1]]$jn[[1]]$plot)


  return(BM)
}
