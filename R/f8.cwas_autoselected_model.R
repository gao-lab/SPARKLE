#' Automatically Select Best Model and Generate Forest Plot with Odds Ratios
#'
#' This function automatically selects the best model from a list of models and generates a forest plot with odds ratios (ORs).
#'
#' @param mdlist A list of models obtained from previous analyses.
#' @param method Method for selecting the best model ("Best", "Equally Best", or "Integration") (default is "Best").
#'   - "Best": Selects the model with the lowest AIC (Akaike Information Criterion) value as the best model.
#'   - "Equally Best": Computes multiple models with equivalent AIC values based on KL (Kullback-Leibler) information, then selects the model with corrected positive results.
#'   - "Integration": Integrates all credible models to derive a combined inference.
#' @param interaction Logical indicating whether to consider interaction effects (default is FALSE).
#'
#' @return A list containing the selected best model and additional information.
#'
#' @export
#'
#' @examples  cwas_autoselected_model(cwas.test.data.single.model,method="Best")

cwas_autoselected_model <- function(mdlist,method=c("Best"),interaction=F){



  #usethis::use_package(package="tidyr",type="Import")
  #usethis::use_package(package="FactoMineR",type="Import")
  #usethis::use_package(package="leaps",type="Import")
  #usethis::use_package(package="lme4",type="Import")
  #usethis::use_package(package="stats",type="Import")
  #usethis::use_package(package="ggrepel",type="Import")
  #usethis::use_package(package="gtsummary",type="Import")
  #usethis::use_package(package="grid",type="Import")
  #usethis::use_package(package="forestploter",type="Import")
  #usethis::use_package(package="magrittr",type="Import")
  #usethis::use_package(package="MuMIn",type="Import")



  if(interaction==F){### 非交互项

    print("All models Selection")

    if(method=="Best"){

      mymodel <- autoselectedmodel_Best(mdlist)
      #OR_plot(mymodel[[2]],Tablename =method )

      tryCatch( mymodel$OR_plot_result <-OR_plot(data = mymodel,Tablename = paste(method,"model for all cell types") ), error = function(e) NULL)
      return(mymodel)

    }else if(method=="Equally Best"){
      mymodel <- autoselectedmodel_Equally_Best(mdlist)
      tryCatch( mymodel$OR_plot_result <-  OR_plot(mymodel,Tablename = paste(method,"model for all cell types") ), error = function(e) NULL)
      return(mymodel)

    }else if(method=="Integration"){
      mymodel <- autoselectedmodel_Integration(mdlist)
      mymodel <-list(NULL,mymodel)
      tryCatch( mymodel$OR_plot_result <-  OR_plot(mymodel,Tablename = paste(method,"model for all cell types") ), error = function(e) NULL)
      return(mymodel)

    }else{
      print("please selected method from Best or Equally Best or Integration")

    }

  }else{### 交互项

    print("All interaction models Selection")

    if(method=="Best"){

      mymodel <- autoselectedmodel_Best(mdlist,effect_cell_rate="rate:rate2")

      tryCatch(mymodel$OR_plot_result <- OR_plot(mymodel,Tablename =method ), error = function(e) NULL)
      return(mymodel)

    }else if(method=="Equally Best"){
      mymodel <- autoselectedmodel_Equally_Best(mdlist,effect_cell_rate="rate:rate2")
      tryCatch( mymodel$OR_plot_result <-  OR_plot(mymodel,Tablename =method ), error = function(e) NULL)
      return(mymodel)

    }else if(method=="Integration"){
      mymodel <- autoselectedmodel_Integration(mdlist,effect_cell_rate="rate:rate2")
      mymodel <-list(NULL,mymodel)
      tryCatch(mymodel$OR_plot_result <-  OR_plot(mymodel[[2]],Tablename =method ), error = function(e) NULL)
      return(mymodel)

    }else{
      print("please selected method from Best or Equally Best or Integration")

    }




  }




}




