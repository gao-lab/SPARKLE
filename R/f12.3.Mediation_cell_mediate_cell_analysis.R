
cwas_cell_mediate_cell_analysis <- function(cwas.data,Best_model,X.cell,selected_celltype=NULL,fix_effect,random_effect){

  if(length(X.cell)>1){   stop("Please Input only one X.cell label or put into the selected_celltype!",     call. = FALSE)}

  if(is.null(X.cell)){

    stop("Please Input only one X.cell label !",     call. = FALSE)

  }else{

    if(is.null(selected_celltype)){ selected_celltype <- unique(cwas.data$Celltype)}
    selected_celltype <- selected_celltype[!selected_celltype==X.cell]
    cellcom <-  data.frame(X.cell,selected_celltype)
    colnames(cellcom) <- c("X.cell","Moderation")
    cellcom$Name <- paste0("X:",cellcom[,1]," Mediation:",cellcom[,2])
    modlist <- list()
    res.all <- c()
    for(i in 1:length(cellcom$X.cell)){
      print(cellcom[i,])
      modlist[[i]] <-  mediation_analysis_auto(cwas.data=cwas.data,Best_model=Best_model,X.cell=cellcom[i,1],Moderation.cell = cellcom[i,2], fix_effect = fix_effect,
                                               random_effect = random_effect)
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


  }


  sf.rds(final.mod, "09cwas_cell_mediation_cell_analysis")

  return(final.mod)

}
