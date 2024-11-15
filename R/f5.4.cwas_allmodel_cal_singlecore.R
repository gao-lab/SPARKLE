
cwas_allmodel_cal_singlecore <- function(cwas.data,link_function="binomial",method="mix_effect"){



  if (any(is.na(cwas.data))) {

    stop("Stop:There are missing values in your input data", call. = FALSE)

  }

  attr.info <- attributes(cwas.data)
  genelist <-  attr.info$genelist
  cwas.data <- cwas.data[,!colnames(cwas.data) %in%genelist]
  mdlist <- list()

  allcelltype <- unique(cwas.data$Celltype)


  for (i in 1:length(allcelltype)) {

    #celldf <- cwas.data %>% dplyr::select_all() %>% dplyr::filter(Celltype==allcelltype[i])

    celldf <- data.frame()
    celldf <- cwas.data[cwas.data$Celltype == allcelltype[[i]], ]

    mdlist[[i]]="null"

    #print(length(celldf$Sample))
    if(method=="mix_effect"){

      # 尝试执行混合模型
      tryCatch({
        mdlist[[i]] <-  cwas_model_cal_mix(cwas.data = celldf, selectedcelltype = allcelltype[i], link = link_function)
      }, error = function(e) {
        # 如果出现错误，提示无法完成混合模型计算,执行固定效应模型
        cat("Error: Mixed effect model was not applicable for this data. Switching to fixed model.\n")

      })


      judegement <- mdlist[[i]]=="NotSuccess"

      if(judegement[1]){
        mdlist[[i]]<- cwas_model_cal_fix(cwas.data = celldf, selectedcelltype = allcelltype[i], link = link_function)
      }

      print(i)
    }else if(method=="fix_effect"){

      mdlist[[i]]<- cwas_model_cal_fix(cwas.data = celldf, selectedcelltype = allcelltype[i], link = link_function)

    }else if(method=="all"){

      mdlist[[i]] <-  cwas_model_cal_mix2(cwas.data = celldf, selectedcelltype = allcelltype[i], link = link_function)

    }else{

      print("Please select method for [mix_effect],[fix_effect],or[all]")

    }
    mdlist[[i]]$inputdata <- cwas.data

  }

  names( mdlist) <- allcelltype


  sf.rds(mdlist,file = "02CWAS_all_single_model")

  return(mdlist)
}
