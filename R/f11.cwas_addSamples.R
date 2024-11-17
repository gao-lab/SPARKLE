#' Simulate Power Curve with Increasing Sample Size
#'
#'This function simulates a power curve by varying the sample size (X) and the power (Y) based on specified parameters.
#'
#' @param cwas.data.selected A list containing the selected models of CWAS analysis.
#' @param selected_celltype  A character vector specifying the cell types of interest.
#' @param yournsim   Number of simulations to perform for each sample size increment (default is 10).
#' @param OR_cutoff_up Upper cutoff value for odds ratio (OR) (optional).
#' @param OR_cutoff_down Lower cutoff value for odds ratio (OR) (optional).
#' @param CCI Logical indicating whether to consider cell change index (CCI) (default is TRUE).
#' @param interaction Logical indicating whether to consider interaction effects (default is FALSE).
#' @param increased_sample_info The method to increase sample size, such as "rate" or other specified factors.
#'
#' @return A simulated power curve plot for each specified cell type with varying sample sizes.
#' @export
#'

cwas_addSamples <- function(cwas.data.selected,selected_celltype,yournsim=10,OR_cutoff_up=NULL,OR_cutoff_down=NULL,CCI=T,interaction=F,increased_sample_info="rate"){


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
  #usethis::use_package(package="simr",type="Import")

  cell_types <- cwas.data.selected[[2]]$celltype
  pclist <- list()


  for(i in 1:length(selected_celltype)){
    cell_index <- which(cell_types == selected_celltype[i])
    fit <- cwas.data.selected[[1]][[cell_index]]
    environment(attributes(fit@frame)$formula) <- globalenv()
    celldf <<-fit@frame

    # 计算均值差值
   # result_diff <<- fit@frame %>%  group_by(Phenotype) %>%    summarize(mean_rate = mean(rate))

    grouped_data <- split(fit@frame, fit@frame$Phenotype)

    # 计算每个组的平均 rate
    result_diff <- lapply(grouped_data, function(df) {
      mean_rate <- mean(df$rate)
      data.frame(Phenotype = unique(df$Phenotype), mean_rate = mean_rate)
    })

    # 将结果合并为一个数据框
    result_diff <- do.call(rbind, result_diff)

    # 将结果赋给 result_diff
    result_diff <- as.data.frame(result_diff)


    # 是否考虑细胞变化量
    ifelse(CCI, mean_difference <<- diff(result_diff$mean_rate),mean_difference<<-1)

    if(interaction){
      object <- "rate:rate2"
      if(!is.null(OR_cutoff_up)){
        if(fixef(fit)[object]>0){
          fixef(fit)[object] <- log(OR_cutoff_up)/mean_difference #指定固定效应x为-0.1
          print(paste(object,"risk were set as",OR_cutoff_up))
        }else{
          fixef(fit)[object] <- log(OR_cutoff_down)/mean_difference #指定固定效应x为-0.1
          print(paste(object,"risk were set as",OR_cutoff_down))
        }
      }

    }else{
      object <- "rate"
      if(!is.null(OR_cutoff_up)){
        if(fixef(fit)[object]>0){
          fixef(fit)[object] <- log(OR_cutoff_up)/mean_difference #指定固定效应x为-0.1
          print(paste(object,"risk were set as",OR_cutoff_up))
        }else{
          fixef(fit)[object] <- log(OR_cutoff_down)/mean_difference #指定固定效应x为-0.1
          print(paste(object,"risk were set as",OR_cutoff_down))
        }
      }
    }

    model_ext_subj <- simr::extend(fit, within=increased_sample_info, n=10)
    pclist[[i]] <- simr::powerCurve(model_ext_subj, within=increased_sample_info, breaks=1:10)
    plot(pclist[[i]])
  }

  sf.rds(pclist,"05CWAS_add_sample")

}


