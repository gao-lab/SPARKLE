

#### F3.2.1 基因-细胞 中介 #####

cwas_gene_mediate_cell_analysis_singlecore <- function(cwas.data,Best_model,X.cell,Mediation=NULL,fix_effect,random_effect){


  attri_info <- attributes(cwas.data)

  if(length(X.cell)>1){   stop("Please Input only one X.cell label or put into the selected_celltype!",     call. = FALSE)}

  if(is.null(X.cell)){   stop("Please Input only one X.cell label !",     call. = FALSE)   }

  if (is.null(Mediation)) {
    if (!is.null(attri_info$geneset_score)) {
      Mediation <- attri_info$geneset_score
    } else if (!is.null(attri_info$genelist)) {
      Mediation <- attri_info$genelist
    } else {
      stop("Please Input Mediation label!", call. = FALSE)
    }
  }



  cellcom <-  data.frame(X.cell,Mediation)
  colnames(cellcom) <- c("X.cell","Mediation")
  cellcom$Name <- paste0("X:",cellcom[,1]," Mediation:",cellcom[,2])
  modlist <- list()
  res.all <- c()

  for (i in 1:length(cellcom$Mediation)) {
    print(cellcom[i,])

    error_occurred <- FALSE  # 初始化错误标志变量

    modlist[[i]] <- tryCatch({
      mediation_analysis_auto(
        cwas.data = cwas.data,
        Best_model = Best_model,
        X.cell = as.character( cellcom[i, 1]),
        Moderation.cell = as.character(cellcom[i, 1]),
        Moderation.gene = as.character(cellcom[i, 2]),
        fix_effect = fix_effect,
        random_effect = random_effect
      )
    }, error = function(e) {
      # 打印错误信息并设置错误标志变量
      message(paste0("Error in iteration ", i, ": ", e$message))
      error_occurred <<- TRUE
      return(NULL)
    })

    # 如果发生错误，跳过本次循环
    if (error_occurred) {
      next
    }

    print(i)
    print(paste0(cellcom$Name[i], " Model calculation have been done"))

    res <- modlist[[i]]$results[[1]]$mediation

    res.reshaped <- data.frame(
      cellcom$Name[i],
      res$Effect[1],
      res[1, 8],
      res$pval[1],
      res$pval[2],
      res$pval[3]
    )

    colnames(res.reshaped) <- c("Model", "Effect", "CI", "Indirect", "Direct", "Total")

    res.all <- rbind(res.all, res.reshaped)
  }



  names(modlist) <- cellcom$Name
  res.all <- res.all[order(res.all$Indirect), ]
  pp <- Mediation_forest_plot(res.all)
  print(pp)
  final.mod <- list(modlist,res.all,pp)
  names(final.mod) <- c("Allresults","Summary","Forest_plot")


 ###sf.rds(final.mod, "09cwas_gene_mediation_cell_analysis")

  return(final.mod)

}

#### F3.2.2  基因-细胞 中介（并行计算） #####


library(doParallel)
library(foreach)
library(plyr) # Load plyr package for rbind.fill function

cwas_gene_mediate_cell_analysis_para <- function(cwas.data, Best_model, X.cell, Mediation = NULL, num_cores = 4,fix_effect,random_effect) {

  library(doParallel)
  library(foreach)
  library(plyr)
  attri_info <- attributes(cwas.data)

  if (length(X.cell) > 1) {
    stop("Please Input only one X.cell label or put into the selected_celltype!", call. = FALSE)
  }

  if (is.null(X.cell)) {
    stop("Please Input only one X.cell label !", call. = FALSE)
  }

  if (is.null(Mediation)) {
    if (!is.null(attri_info$geneset_score)) {
      Mediation <- attri_info$geneset_score
    } else if (!is.null(attri_info$genelist)) {
      Mediation <- attri_info$genelist
    } else {
      stop("Please Input Mediation label!", call. = FALSE)
    }
  }

  cellcom <- data.frame(X.cell, Mediation)
  colnames(cellcom) <- c("X.cell", "Mediation")
  cellcom$Name <- paste0("X:", cellcom[, 1], " Mediation:", cellcom[, 2])
  modlist <- list()
  res.all <- data.frame()

  # Define the function for parallel execution
  analysis_func <- function(i) {
    result <- list()
    error_occurred <- FALSE

    mod <- tryCatch({
      mediation_analysis_auto(
        cwas.data = cwas.data,
        Best_model = Best_model,
        X.cell = cellcom[i, 1],
        Moderation.cell = cellcom[i, 1],
        Moderation.gene = cellcom[i, 2],
        fix_effect = fix_effect,
        random_effect = random_effect
      )
    }, error = function(e) {
      message(paste0("Error in iteration ", i, ": ", e$message))
      error_occurred <<- TRUE
      return(NULL)
    })

    if (error_occurred) {
      return(NULL)
    }

    result$mod <- mod
    result$res <- mod$results[[1]]$mediation

    # Ensure consistent structure by padding missing values if necessary
    reshaped <- data.frame(
      Model = cellcom$Name[i],
      Effect = ifelse(length(result$res$Effect) >= 1, result$res$Effect[1], NA),
      CI = ifelse(ncol(result$res) >= 8, result$res[1, 8], NA),
      Indirect = ifelse(length(result$res$pval) >= 1, result$res$pval[1], NA),
      Direct = ifelse(length(result$res$pval) >= 2, result$res$pval[2], NA),
      Total = ifelse(length(result$res$pval) >= 3, result$res$pval[3], NA)
    )

    return(as.data.frame(reshaped))
  }

  # Setup parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Export necessary data and functions to the cluster
 # clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "mediation_analysis_auto","fix_effect","random_effect"))
  parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "mediation_analysis_auto", "fix_effect", "random_effect"), envir = environment())
  # Run the analysis in parallel
  results <- foreach::foreach(i = 1:nrow(cellcom), .packages = c("doParallel", "foreach")) %dopar% {
    analysis_func(i)
  }

  # Stop the cluster
  parallel::stopCluster(cl)

  # Filter out NULL results and combine the rest
  results <- Filter(Negate(is.null), results)
  if (length(results) > 0) {
    res.all <- do.call(rbind.fill, results)
  }

  # Sorting the results
  res.all <- res.all[order(res.all$Indirect), ]

  # Generate Mediation Effect Plot
  pp <- Mediation_forest_plot(res.all)
  print(pp)

  # Prepare final result
  final.mod <- list(modlist, res.all, pp)
  names(final.mod) <- c("Allresults", "Summary", "Forest_plot")

 ###sf.rds(final.mod, "10cwas_gene_mediation_cell_analysis_Parallel")
  return(final.mod)
}

#cwas_gene_mediate_cell_analysis_para(cwas.data = cwas.test.data,Best_model = cwas.test.data.single.model.best,X.cell = "MTF",Mediation = "TGFB1",num_cores = 2,fix_effect = "Group",random_effect = "Subgroup")

#### 3.3.1  细胞-基因 中介（顺序计算） #####


cwas_cell_mediate_gene_analysis_singlecore <- function(cwas.data,Best_model,X.gene,selected_celltype=NULL,fix_effect,random_effect){


  if(length(X.gene)>1){   stop("Please Input only one X.gene label  !",     call. = FALSE)}

  if(is.null(X.gene)){    stop("Please Input only one X.gene label !",     call. = FALSE)  }

  if(is.null(selected_celltype)){ selected_celltype <- unique(cwas.data$Celltype)}

  cellcom <-  data.frame(X.gene,selected_celltype)
  colnames(cellcom) <- c("X.gene","Moderation")
  cellcom$Name <- paste0("X:",cellcom[,1]," Mediation:",cellcom[,2])
  modlist <- list()
  res.all <- c()
  for(i in 1:length(cellcom$X.gene)){
    print(cellcom[i,])
    modlist[[i]] <-  mediation_analysis_auto(cwas.data=cwas.data,Best_model=Best_model,X.gene=cellcom[i,1],Moderation.cell = cellcom[i,2],fix_effect=fix_effect,random_effect=random_effect)
    print(i)
    print(paste0(cellcom$Name[i] ," Model calculation have been done"))
    res <- modlist[[i]]$results[[1]]$mediation
    res.reshaped <- data.frame(cellcom$Name[i],res$Effect[1],res[1,8],res$pval[1],res$pval[2],res$pval[3])
    colnames(res.reshaped) <- c("Model","Effect","CI","Indirect","Direct","Total")
    res.all <- rbind(res.all,res.reshaped)
  }

  names(modlist) <- cellcom$Name
  res.all <- res.all[order(res.all$Indirect), ]
  pp <- Mediation_forest_plot(res.all)
  print(pp)
  final.mod <- list(modlist,res.all,pp)
  names(final.mod) <- c("Allresults","Summary","Forest_plot")




 ###sf.rds(final.mod, "09cwas_cell_mediation_cell_analysis")

  return(final.mod)

}


#### 3.3.2  细胞-基因 中介（并行计算） #####

library(doParallel)
library(foreach)
library(plyr)  # Load plyr package for rbind.fill function

cwas_cell_mediate_gene_analysis_para <- function(cwas.data, Best_model, X.gene, selected_celltype = NULL, num_cores = 4,fix_effect,random_effect) {

  library(doParallel)
  library(foreach)
  library(plyr)

  if (is.null(X.gene)) {
    stop("Please Input only one X.gene label!", call. = FALSE)
  }

  if (is.null(selected_celltype)) {
    selected_celltype <- unique(cwas.data$Celltype)
  }

  #cellcom <- data.frame(X.gene, selected_celltype)
  cellcom <- expand.grid(X.gene = X.gene, selected_celltype = selected_celltype)
  colnames(cellcom) <- c("X.gene", "Mediation")
  cellcom$Name <- paste0("X:", cellcom[, 1], " Mediation:", cellcom[, 2])
  modlist <- list()
  res.all <- data.frame()


  # Define the function for parallel execution
  analysis_func <- function(i) {
    result <- list()
    error_occurred <- FALSE

    mod <- tryCatch({
      mediation_analysis_auto(
        cwas.data = cwas.data,
        Best_model = Best_model,
        X.gene = as.character(cellcom[i, 1]),
        Moderation.cell = cellcom[i, 2],
        fix_effect = fix_effect,
        random_effect = random_effect
      )
    }, error = function(e) {
      message(paste0("Error in iteration ", i, ": ", e$message))
      error_occurred <<- TRUE
      return(NULL)
    })

    if (error_occurred) {
      return(NULL)
    }

    result$mod <- mod
    result$res <- mod$results[[1]]$mediation

    # Ensure consistent structure by padding missing values if necessary
    reshaped <- data.frame(
      Model = cellcom$Name[i],
      Effect = ifelse(length(result$res$Effect) >= 1, result$res$Effect[1], NA),
      CI = ifelse(ncol(result$res) >= 8, result$res[1, 8], NA),
      Indirect = ifelse(length(result$res$pval) >= 1, result$res$pval[1], NA),
      Direct = ifelse(length(result$res$pval) >= 2, result$res$pval[2], NA),
      Total = ifelse(length(result$res$pval) >= 3, result$res$pval[3], NA)
    )

    return(as.data.frame(reshaped))
  }

  # Setup parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Export necessary data and functions to the cluster
  #clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "mediation_analysis_auto","fix_effect","random_effect"))
  parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "mediation_analysis_auto", "fix_effect", "random_effect"), envir = environment())
  # Run the analysis in parallel
  results <- foreach::foreach(i = 1:nrow(cellcom), .packages = c("plyr", "foreach")) %dopar% {
    analysis_func(i)
  }

  # Stop the cluster
  parallel::stopCluster(cl)

  # Filter out NULL results and combine the rest
  results <- Filter(Negate(is.null), results)
  if (length(results) > 0) {
    res.all <- do.call(rbind.fill, results)
  }

  # Sorting the results
  res.all <- res.all[order(res.all$Indirect), ]

  # Generate Mediation Effect Plot
  pp <- Mediation_forest_plot(res.all)
  print(pp)

  # Prepare final result
  final.mod <- list(modlist, res.all, pp)
  names(final.mod) <- c("Allresults", "Summary", "Forest_plot")

  # Save the result
 ###sf.rds(final.mod, "09cwas_cell_mediation_cell_analysis")

  return(final.mod)
}
