
##### F3.1 细胞-细胞 调节效应计算  #####


cwas_cell_moderation_cell_analysis <- function(cwas.data,Best_model,X.cell=NULL,selected_celltype=NULL){

  if(length(X.cell)>1){   stop("Please Input only one X.cell label or put into the selected_celltype!",     call. = FALSE)}

  if(is.null(X.cell)){

    final.mod <-  cell_moderation_cell_analysis_all(cwas.data,Best_model,selected_celltype =selected_celltype)

  }else{

    if(is.null(selected_celltype)){ selected_celltype <- unique(cwas.data$Celltype)}
    selected_celltype <- selected_celltype[!selected_celltype==X.cell]
    cellcom <-  data.frame(X.cell,selected_celltype)
    colnames(cellcom) <- c("X.cell","Moderation")
    cellcom$Name <- paste0("X:",cellcom[,1]," Moderation:",cellcom[,2])
    modlist <- list()
    res.all <- c()
    for(i in 1:length(cellcom$X.cell)){
      print(cellcom[i,])
      modlist[[i]] <-  moderation_analysis_auto(cwas.data=cwas.data,Best_model=Best_model,X.cell=cellcom[i,1],Moderation.cell = cellcom[i,2])
      print(i)
      print(paste0(cellcom$Name[i] ," Model calculation have been done"))
      res <- modlist[[i]]$results[[1]]$simple.slopes
      res.reshaped <- data.frame(cellcom$Name[i],res$Effect[2],res$`[95% CI]`[2],res$pval[1],res$pval[2],res$pval[3])
      colnames(res.reshaped) <- c("Model","Effect","CI","Pvalue(-SD)","Pvalue(Mean)","Pvalue(+SD)")
      res.all <- rbind(res.all,res.reshaped)
    }

    names(modlist) <- cellcom$Name
    res.all <- res.all[order(res.all$`Pvalue(Mean)`), ]
    pp <- Moderation_forest_plot(res.all)
    print(pp)
    final.mod <- list(modlist,res.all,pp)
    names(final.mod) <- c("Allresults","Summary","Forest_plot")


  }


  ###sf.rds(final.mod, "08cwas_cell_moderation_cell_analysis")

  return(final.mod)

}


#### F3.2.1 基因-细胞 调节效应（顺序计算） #####


cwas_gene_moderation_cell_analysis <-  function(cwas.data,Best_model,X.cell,Moderation.cell,selected_gene=NULL){

  res.all <- c()
  cellcom <-  data.frame(X.cell,Moderation.cell,selected_gene)
  colnames(cellcom) <- c("X.cell","Moderation.cell","Moderation.gene")
  cellcom$Name <- paste0("X:",cellcom$X.cell," Moderation gene:",cellcom$Moderation.gene,"(",cellcom$Moderation.cell,")")
  modlist <- list()
  for(i in 1:length(selected_gene)){
    print(cellcom[i,])
    modlist[[i]] <-  moderation_analysis_auto(cwas.data=cwas.data,Best_model=Best_model,X.cell=cellcom[i,1],Moderation.cell = cellcom[i,2],Moderation.gene =cellcom[i,3] )
    print(i)
    print(paste0(cellcom$Name[i] ," Model calculation have been done"))
    res <- modlist[[i]]$results[[1]]$simple.slopes
    res.reshaped <- data.frame(cellcom$Name[i],res$Effect[2],res$`[95% CI]`[2],res$pval[1],res$pval[2],res$pval[3])
    colnames(res.reshaped) <- c("Model","Effect","CI","Pvalue(-SD)","Pvalue(Mean)","Pvalue(+SD)")
    res.all <- rbind(res.all,res.reshaped)
  }

  names(modlist) <- cellcom$Name
  res.all <- res.all[order(res.all$`Pvalue(Mean)`), ]
  pp <- Moderation_forest_plot(res.all)
  print(pp)
  final.mod <- list(modlist,res.all,pp)
  names(final.mod) <- c("Allresults","Summary","Forest_plot")

  ###sf.rds(final.mod, "09cwas_gene_moderation_cell_analysis")

  return(final.mod)

}


#### F3.2.2 基因-细胞 调节效应（并行计算） #####


cwas_gene_moderation_cell_analysis_para <- function(cwas.data, Best_model, X.cell, Moderation.cell, selected_gene = NULL, num_cores = 4) {

  library(doParallel)
  library(foreach)
  library(plyr)


  # Validate input
  if (is.null(selected_gene)) {
    stop("Please provide selected_gene list!", call. = FALSE)
  }

  # Prepare combinations
  cellcom <- data.frame(X.cell, Moderation.cell, selected_gene)
  colnames(cellcom) <- c("X.cell", "Moderation.cell", "Moderation.gene")
  cellcom$Name <- paste0("X:", cellcom$X.cell, " Moderation gene:", cellcom$Moderation.gene, "(", cellcom$Moderation.cell, ")")

  # Placeholder for results
  modlist <- list()
  res.all <- data.frame()

  # Define the function to be executed in parallel
  analysis_func <- function(i) {
    result <- list()
    error_occurred <- FALSE

    # Perform moderation analysis
    mod <- tryCatch({
      moderation_analysis_auto(
        cwas.data = cwas.data,
        Best_model = Best_model,
        X.cell = cellcom[i, 1],
        Moderation.cell = cellcom[i, 2],
        Moderation.gene = cellcom[i, 3]
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
    result$res <- mod$results[[1]]$simple.slopes

    # Ensure consistent structure by padding missing values if necessary
    reshaped <- data.frame(
      Model = cellcom$Name[i],
      Effect = ifelse(length(result$res$Effect) >= 2, result$res$Effect[2], NA),
      CI = ifelse(ncol(result$res) >= 3, result$res$`[95% CI]`[2], NA),
      `Pvalue(-SD)` = ifelse(length(result$res$pval) >= 1, result$res$pval[1], NA),
      `Pvalue(Mean)` = ifelse(length(result$res$pval) >= 2, result$res$pval[2], NA),
      `Pvalue(+SD)` = ifelse(length(result$res$pval) >= 3, result$res$pval[3], NA)
    )

    return(as.data.frame(reshaped))
  }

  # Setup parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Export necessary data and functions to the cluster
  #parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "moderation_analysis_auto"))
  parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "moderation_analysis_auto"), envir = environment())


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
  res.all$`Pvalue(Mean)` <- res.all$Pvalue.Mean.
  res.all <- res.all[order(res.all$`Pvalue(Mean)`), ]

  # Generate Effect Plot
  pp <- Moderation_forest_plot(res.all)
  print(pp)

  # Prepare final result
  final.mod <- list(modlist, res.all, pp)
  names(final.mod) <- c("Allresults", "Summary", "Forest_plot")

  # Save the result
  ###sf.rds(final.mod, "09cwas_gene_moderation_cell_analysis")

  return(final.mod)
}

#### F3.3.1 细胞-基因 调节效应（顺序计算） #####

cwas_cell_moderation_gene_analysis <-  function(cwas.data,Best_model,X.gene,selected_celltype=NULL){
  res.all <- c()
  if(is.null(selected_celltype)){ selected_celltype <- unique(cwas.data$Celltype)}
  cellcom <-  data.frame(X.gene, selected_celltype)
  colnames(cellcom) <- c("X.gene","Moderation.cell")
  cellcom$Name <- paste0("X:",cellcom$X.gene," Moderation cell:",cellcom$Moderation.cell)

  modlist <- list()
  for(i in 1:length(selected_celltype)){
    print(cellcom[i,])
    modlist[[i]] <-  moderation_analysis_auto(cwas.data=cwas.data,Best_model=Best_model,X.gene = cellcom[i,1],Moderation.cell = cellcom[i,2])
    print(i)
    print(paste0(cellcom$Name[i] ," Model calculation have been done"))
    res <- modlist[[i]]$results[[1]]$simple.slopes
    res.reshaped <- data.frame(cellcom$Name[i],res$Effect[2],res$`[95% CI]`[2],res$pval[1],res$pval[2],res$pval[3])
    colnames(res.reshaped) <- c("Model","Effect","CI","Pvalue(-SD)","Pvalue(Mean)","Pvalue(+SD)")
    res.all <- rbind(res.all,res.reshaped)
  }

  names(modlist) <- cellcom$Name
  res.all <- res.all[order(res.all$`Pvalue(Mean)`), ]
  pp <- Moderation_forest_plot(res.all)
  print(pp)
  final.mod <- list(modlist,res.all,pp)
  names(final.mod) <- c("Allresults","Summary","Forest_plot")

  ###sf.rds(final.mod, "09cwas_gene_moderation_cell_analysis")

  return(final.mod)

}


#### F3.3.2 细胞-基因 调节效应（并行计算） #####


# library(parallel)
# library(doParallel)
# library(plyr)

cwas_cell_moderation_gene_analysis_para <- function(cwas.data, Best_model, X.gene, selected_celltype = NULL, num_cores = 4) {

  library(doParallel)
  library(foreach)
  library(plyr)

  if (is.null(X.gene)) {
    stop("Please Input only one X.gene label!", call. = FALSE)
  }

  if (is.null(selected_celltype)) {
    selected_celltype <- unique(cwas.data$Celltype)
  }

  # Generate the full combination of X.gene and selected_celltype
  cellcom  <- expand.grid(X.gene = X.gene, Moderation.cell = selected_celltype)
  cellcom$Name <- paste0("X:", cellcom$X.gene, " Moderation cell:", cellcom$Moderation.cell)
  modlist <- list()
  res.all <- data.frame()



  # Define the function for parallel execution
  analysis_func <- function(i) {
    result <- list()
    error_occurred <- FALSE

    mod <- tryCatch({
      moderation_analysis_auto(
        cwas.data = cwas.data,
        Best_model = Best_model,
        X.gene = as.character(cellcom[i, 1]),
        Moderation.cell = cellcom[i, 2]
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
    result$res <- mod$results[[1]]$simple.slopes

    # Ensure consistent structure by padding missing values if necessary
    reshaped <- data.frame(
      Model = cellcom$Name[i],
      Effect = ifelse(length(result$res$Effect) >= 2, result$res$Effect[2], NA),
      CI = ifelse(ncol(result$res) >= 3, result$res$`[95% CI]`[2], NA),
      `Pvalue(-SD)` = ifelse(length(result$res$pval) >= 1, result$res$pval[1], NA),
      `Pvalue(Mean)` = ifelse(length(result$res$pval) >= 2, result$res$pval[2], NA),
      `Pvalue(+SD)` = ifelse(length(result$res$pval) >= 3, result$res$pval[3], NA)
    )

    return(as.data.frame(reshaped))
  }

  # Setup parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Export necessary data and functions to the cluster
 # parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "moderation_analysis_auto"))
  parallel::clusterExport(cl, varlist = c("cwas.data", "Best_model", "cellcom", "moderation_analysis_auto"), envir = environment())
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
  res.all$`Pvalue(Mean)` <- res.all$Pvalue.Mean.
  res.all <- res.all[order(res.all$`Pvalue(Mean)`), ]

  # Generate Effect Plot
  pp <- Moderation_forest_plot(res.all)
  print(pp)

  # Prepare final result
  final.mod <- list(modlist, res.all, pp)
  names(final.mod) <- c("Allresults", "Summary", "Forest_plot")

  # Save the result
  ###sf.rds(final.mod, "09cwas_gene_moderation_cell_analysis")

  return(final.mod)
}

################################################################ F4 调节效应（并行计算） #############################################################################################


