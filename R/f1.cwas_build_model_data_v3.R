#' Prepare CWAS Data for Modeling
#'
#' This function prepares input data containing sample metadata, phenotype, cell type, and optional additional variables into a standard format suitable for subsequent analysis. The input can be a Seurat object or a metadata dataframe. It extracts necessary columns and computes cell rates if not provided.
#'
#' @param inputdata A Seurat object or a metadata dataframe containing sample metadata, phenotype, cell type, and optional additional variables.
#' @param Sample The column name specifying sample identifiers.
#' @param Phenotype The column name specifying the phenotype variable.
#' @param Celltype The column name specifying the cell type variable.
#' @param Group (Optional) The column name specifying group information.
#' @param Subgroup (Optional) The column name specifying subgroup information.
#' @param Covariate1 (Optional) The column name specifying the first covariate.
#' @param Covariate2 (Optional) The column name specifying the second covariate.
#' @param selected_celltype (Optional) A vector of specific cell types to include in the analysis.
#' @param Control_label The label representing control samples (default is "Control").
#' @param Disease_label The label representing disease samples (default is "Disease").
#' @param Cellrate (Optional) The column name specifying precomputed cell rates (if available).
#'
#' @return A standardized dataframe suitable for CWAS analysis, containing sample identifiers, phenotype, group information, subgroup information, covariates, and cell rates.
#'
#' @export
#'
#' @examples  cwas_build_model_data(inputdata = cwas.test,Sample = "GSM",Phenotype = "sample.condition",Celltype = "Subtype",Group = "GSE",Control_label="lymph.node",Disease_label="tumor",Cellrate = "Fraction")



#inputdata$celltype
#inputdata <- readRDS("I:\\BaiduSyncdisk\\DT20230523已备份\\DT20230801\\20231104CWAS\\20231212axisFrame\\重新注释sc20231217.rds")

# genelist=c("COL1A1","TGFBR1")
#
# Sample="orig.ident"
# Phenotype="Pathology"
# Control_label="Health"
# Disease_label="Fibrosis"
# Group="Tissue"
# Subgroup="GSE"
# selected_celltype=NULL
# Covariate1=NULL
# Covariate2=NULL
# Cellrate=NULL
# Celltype="celltype"

cwas_build_model_data <- function(inputdata,Sample="orig.ident",Phenotype,Celltype,Group=NULL,Subgroup=NULL,Covariate1=NULL,Covariate2=NULL,selected_celltype=NULL,Control_label="Control",Disease_label="Disease",Cellrate=NULL,genelist=NULL,geneset_score=NULL,add_therapeutic_targets=F){

  #usethis::use_package(package="dplyr",type="Import")

  data_path <- system.file("extdata", "TTDdatabase240523.rdata", package = "SPARKLE")
  system.file(package = "SPRAKLE")
  # 检查文件是否存在
  if (!file.exists(data_path)) {
    stop("Data file not found in the package, please reinstall the SPARKLE in the right path.")
  }

  # 加载数据集
  load(data_path)

  Target_names <- unique(TTDdatabase$Target_info$Target_gene)


  ##### 判断数据类型是否是dataframe或者Seurat对象 #####

  ##### Seurat对象 #####
  if (inherits(inputdata, "Seurat")) {
    celldfall <- inputdata@meta.data


    celldf <- celldfall[,c(Sample,Phenotype,Celltype)]
    colnames(celldf) <- c("Sample","Phenotype","Celltype")
    if(length(unique(celldf$Phenotype))<2){
      stop("Phenotype less than 1 CANNOT calculate, please provide !",
           call. = FALSE)
    }


    if(add_therapeutic_targets){genelist <- c(genelist,Target_names)}

    celldf.celltype <- dplyr::group_by(celldf, Celltype, Sample)
    celldf.celltype <- dplyr::summarize(celldf.celltype, num = dplyr::n())
    celldf.allnum <- dplyr::group_by(celldf, Sample)
    celldf.allnum <- dplyr::summarize(celldf.allnum, allnum = dplyr::n())
    celldf.celltype.result <- merge(celldf.celltype,celldf.allnum,by="Sample")
    celldf.celltype.result$rate <- celldf.celltype.result$num/ celldf.celltype.result$allnum

    celldf.celltype.result$num=NULL
    celldf.celltype.result$allnum=NULL

    rownames(celldf)=NULL

    if(is.null(Group)){
      print("Warning: No group infomation")
    }else{  celldf$Group <- celldfall[,Group]}

    if(is.null(Subgroup)){
      print("Warning: No Subgroup infomation")
    }else{ celldf$Subgroup <- celldfall[,Subgroup]}

    if(is.null(Covariate1)){
      print("Warning: No Covariate1 infomation")
    }else{ celldf$Cov1 <- celldfall[,Covariate1]}

    if(is.null(Covariate2)){
      print("Warning: No Covariate2 infomation")
    }else{ celldf$Cov2 <- celldfall[,Covariate2]}

    celldf$Celltype=NULL
    celldf <- unique(celldf)

    celldf.celltype.result.meta <- merge(celldf.celltype.result,celldf,by="Sample", all.x = TRUE)


    celldf <- celldf.celltype.result.meta
    celldf$Phenotype <- ifelse(celldf$Phenotype==Control_label,"Control",ifelse(celldf$Phenotype==Disease_label,"Disease","Others"))
    #table(celldf$Phenotype,celldf$Sample)
    celldf <- celldf[celldf$Phenotype != "Others", ]

  } else if (is.data.frame(inputdata)) {

    ##### data.frame对象 #####
    celldfall <- inputdata

    if(is.null(Cellrate)){


      celldf <- celldfall[,c(Sample,Phenotype,Celltype)]
      colnames(celldf) <- c("Sample","Phenotype","Celltype")
      if(length(unique(celldf$Phenotype))<2){
        stop("Phenotype less than 1 CANNOT calculate, please provide !",
             call. = FALSE)
      }


      if(add_therapeutic_targets){genelist <- c(genelist,Target_names)}

      celldf.celltype <- dplyr::group_by(celldf, Celltype, Sample)
      celldf.celltype <- dplyr::summarize(celldf.celltype, num = dplyr::n())
      celldf.allnum <- dplyr::group_by(celldf, Sample)
      celldf.allnum <- dplyr::summarize(celldf.allnum, allnum = dplyr::n())
      celldf.celltype.result <- merge(celldf.celltype,celldf.allnum,by="Sample")
      celldf.celltype.result$rate <- celldf.celltype.result$num/ celldf.celltype.result$allnum

      celldf.celltype.result$num=NULL
      celldf.celltype.result$allnum=NULL

      rownames(celldf)=NULL

      if(is.null(Group)){
        print("Warning: No group infomation")
      }else{  celldf$Group <- celldfall[,Group]}

      if(is.null(Subgroup)){
        print("Warning: No Subgroup infomation")
      }else{ celldf$Subgroup <- celldfall[,Subgroup]}

      if(is.null(Covariate1)){
        print("Warning: No Covariate1 infomation")
      }else{ celldf$Cov1 <- celldfall[,Covariate1]}

      if(is.null(Covariate2)){
        print("Warning: No Covariate2 infomation")
      }else{ celldf$Cov2 <- celldfall[,Covariate2]}

      celldf$Celltype=NULL
      celldf <- unique(celldf)

      celldf.celltype.result.meta <- merge(celldf.celltype.result,celldf,by="Sample", all.x = TRUE)


      celldf <- celldf.celltype.result.meta
      celldf$Phenotype <- ifelse(celldf$Phenotype==Control_label,"Control",ifelse(celldf$Phenotype==Disease_label,"Disease","Others"))
      #table(celldf$Phenotype,celldf$Sample)
      celldf <- celldf[celldf$Phenotype != "Others", ]


    }else{

      celldf <- celldfall[,c(Sample,Phenotype,Celltype,Cellrate)]
      colnames(celldf) <- c("Sample","Phenotype","Celltype","rate")

      if(is.null(Group)){
        print("Warning: No group infomation")
      }else{  celldf$Group <- celldfall[,Group]}

      if(is.null(Subgroup)){
        print("Warning: No Subgroup infomation")
      }else{ celldf$Subgroup <- celldfall[,Subgroup]}

      if(is.null(Covariate1)){
        print("Warning: No Covariate1 infomation")
      }else{ celldf$Cov1 <- celldfall[,Covariate1]}

      if(is.null(Covariate2)){
        print("Warning: No Covariate2 infomation")
      }else{ celldf$Cov2 <- celldfall[,Covariate2]}


      ## celldfPhenotype <- factor(celldf$Phenotype, levels=c(Control_label, Disease_label))

      celldf$Phenotype <- ifelse(celldf$Phenotype==Control_label,"Control",ifelse(celldf$Phenotype==Disease_label,"Disease","Others"))
      celldf <- celldf[celldf$Phenotype != "Others", ]


    }


  }else {
    return("Wrong input format, please input Seruat object or metadata dataframe file !")
  }


  ####


  #### 增加基因信息 ####



  if(is.null(genelist)){print("No gene infomation added")}else{

    if(inherits(inputdata, "Seurat")){

      allgene <- rownames(inputdata)
      diff_vec1 <- setdiff(genelist, allgene)
      if(length(diff_vec1)>0){print(paste(diff_vec1,"are not in the Seurat Object, Please check!")) }
      genelist <- intersect(genelist, allgene)
      #library(Seurat)
      geneinfo <- Seurat::FetchData(object =inputdata, c(Sample,Celltype,genelist))
      colnames(geneinfo) <- c("Sample","Celltype",genelist)

      group_vars <- geneinfo[, c("Sample", "Celltype")]
      # 获取基因表达数据
      gene_data <- geneinfo[, genelist]
      # 计算均值
      result <- stats::aggregate(gene_data, by = list(Sample = group_vars$Sample, Celltype = group_vars$Celltype), FUN = mean)
      colnames(result) <- c("Sample","Celltype",genelist)

      merged_df <- merge(celldf, result, by = c("Sample", "Celltype"))

      celldf <- merged_df

    }else{
      print("No Seruat objected, gene infomation add failed")
    }


  }

  #### 增加基因集打分信息 ####



  if(is.null(geneset_score)){print("No geneset score infomation added")}else{

    if(inherits(inputdata, "Seurat")){

      Seruat.name <- colnames(inputdata@meta.data)
      diff_vec1 <- setdiff(geneset_score, Seruat.name)
      if(length(diff_vec1)>0){print(paste(diff_vec1,"are not in the Seurat Object, Please check!")) }
      geneset_score <- intersect(geneset_score, Seruat.name)
      #library(Seurat)
      geneinfo <- Seurat::FetchData(object =inputdata, c(Sample,Celltype,geneset_score))
      colnames(geneinfo) <- c("Sample","Celltype",geneset_score)

      group_vars <- geneinfo[, c("Sample", "Celltype")]
      # 获取基因表达数据
      gene_data <- geneinfo[, geneset_score]
      # 计算均值
      result <- stats::aggregate(gene_data, by = list(Sample = group_vars$Sample, Celltype = group_vars$Celltype), FUN = mean)
      colnames(result) <- c("Sample","Celltype",geneset_score)

      merged_df <- merge(celldf, result, by = c("Sample", "Celltype"))

      celldf <- merged_df

    }else{
      print("No Seruat objected, gene infomation add failed")
    }


  }

  ##### 增加用户计算的Label信息 #####

  if(is.null(selected_celltype)){}else{celldf <- celldf[which(celldf$Celltype %in% selected_celltype),]}

  myattr <- attributes(celldf)
  myattr$Sample <- Sample
  myattr$Phenotype <- Phenotype
  myattr$Group <- Group
  myattr$Subgroup <- Subgroup
  myattr$Cov1 <- Covariate1
  myattr$Cov2 <- Covariate2
  myattr$Control_label <- Control_label
  myattr$Disease_label <- Disease_label
  myattr$Cellrate <- Cellrate
  myattr$genelist <- intersect(genelist, allgene)
  myattr$geneset_score <- intersect(geneset_score, colnames(inputdata@meta.data))
  attributes(celldf) <- myattr
  ####判断是否有空值####


  if (any(is.na(celldf))) {

    print("Warning:There are missing values in fibro.cwas.data2, add [No label] in the data")

    celldf[is.na(celldf)] <- "No label"

  }


  if (length(unique(celldf$Phenotype))==1) {

    print("Warning:There is only one Phenotype in this data,please check your input !")

  }


  celldf$Celltype <- gsub("-","_",celldf$Celltype)

  celldf <- filter_cell(celldf)

  sf.rds(celldf, "01CWAS_data")

  return(celldf)



}

