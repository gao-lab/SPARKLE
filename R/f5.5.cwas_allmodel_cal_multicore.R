cwas_allmodel_cal_multicore <- function(cwas.data, link_function = "binomial", method = "mix_effect",num_cores=2) {
  if (any(is.na(cwas.data))) {
    stop("Stop: There are missing values in your input data", call. = FALSE)
  }

  attr.info <- attributes(cwas.data)
  genelist <- attr.info$genelist
  cwas.data <- cwas.data[, !colnames(cwas.data) %in% genelist]

  allcelltype <- unique(cwas.data$Celltype)


  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # 导出必要的数据和函数到集群
  parallel::clusterExport(cl, varlist = c("cwas.data", "link_function", "method"),
                          envir = environment())

  mdlist <- foreach::foreach(i = 1:length(allcelltype), .packages = 'SPARKLE') %dopar% {
    celldf <- cwas.data[cwas.data$Celltype == allcelltype[i], ]
    result <- "null"

    if (method == "mix_effect") {
      tryCatch({
        result <- cwas_model_cal_mix(cwas.data = celldf,
                                     selectedcelltype = allcelltype[i], link = link_function)
      }, error = function(e) {
        cat("Error: Mixed effect model was not applicable for this data. Switching to fixed model.\n")
      })

      # 确保结果检查正确
      if (is.character(result) && result == "NotSuccess") {
        result <- cwas_model_cal_fix(cwas.data = celldf,
                                     selectedcelltype = allcelltype[i], link = link_function)
      }
    } else if (method == "fix_effect") {
      result <- cwas_model_cal_fix(cwas.data = celldf,
                                   selectedcelltype = allcelltype[i], link = link_function)
    } else if (method == "all") {
      result <- cwas_model_cal_mix2(cwas.data = celldf,
                                    selectedcelltype = allcelltype[i], link = link_function)
    } else {
      stop("Please select method for [mix_effect],[fix_effect],or[all]")
    }

    result$inputdata <- cwas.data
    return(result)
  }

  parallel::stopCluster(cl)  # 关闭并行计算
  names(mdlist) <- allcelltype
  #mdlist2<<-mdlist
  # 保存结果
  sf.rds(mdlist, file =  "02CWAS_all_single_model")
  return(mdlist)
}
