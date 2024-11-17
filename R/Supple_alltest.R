

consistent_change2 <- function(cwas.test.data,dataname="Fibroblast"){

  cwas.test.data1 <- cwas.test.data
  cwas.test.data1$Group <- c(1)
  cwas.test.data1$Subgroup <- NULL
  info1 <- inconsitancy_calculation1(cwas.test.data1,1)

  cwas.test.data2 <- cwas.test.data
  cwas.test.data2$Subgroup <- NULL
  info2 <- inconsitancy_calculation1(cwas.test.data2,2)

  cwas.test.data3 <- cwas.test.data
  cwas.test.data3$Group <- cwas.test.data3$Subgroup
  cwas.test.data2$Subgroup <- NULL
  info3 <- inconsitancy_calculation1(cwas.test.data3,3)

  cwas.test.data4 <- cwas.test.data
  info4 <- inconsitancy_calculation1(cwas.test.data4,4)

  info <- rbind(info1,info2)

  info <- rbind(info,info3)

  info <- rbind(info,info4)

  info$dataname=dataname

  return(info)

}

inconsitancy_calculation1 <- function(cwas.data,groupname){


  signed_rank <-  function(x){sign(x)*rank(abs(x))}

  cwas.data$Phenotype <- as.factor(cwas.data$Phenotype)

  celltypeall <- unique(cwas.data$Celltype)
  dfall <- c()
  # 线性模型：拟合rate
  for(i in 1:length(celltypeall)){

    cwas.data1 <- cwas.data%>% dplyr::select_all() %>% dplyr::filter(Celltype==celltypeall[i])

    if(groupname==1){
      linear_model <- lm( signed_rank(rate) ~ Phenotype  , data = cwas.data1)
    }else if(groupname==2){
      linear_model <- lm( signed_rank(rate) ~ Phenotype+Group  , data = cwas.data1)
    }else if(groupname==3){
      linear_model <- lm( signed_rank(rate) ~ Phenotype+ Subgroup , data = cwas.data1)
    }else if(groupname==4){
      linear_model <- lm( signed_rank(rate) ~ Phenotype+Group+Subgroup , data = cwas.data1)
    }else{ print("Wrong groupname")}

    linear_model.res <- summary(linear_model)
    linear.P <-  linear_model.res$coefficients["PhenotypeDisease","Pr(>|t|)"]

    logistic_model <- lm( rate  ~ Phenotype, data = cwas.data1)
    logistic_model.res <- summary(logistic_model)
    logistic.P <- logistic_model.res$coefficients["PhenotypeDisease","Pr(>|t|)"]

   df1 <-  data.frame(celltypeall[i],linear.P,logistic.P)
   dfall <- rbind(dfall,df1)


  }



  class <- cwas_classic_comparision(cwas.data)
  wl <- class[[1]][c(1,2,4)]

  all <- cbind(dfall,wl)

  jugep <- function(drop.pvalue2,orig.pvalue){
    if(drop.pvalue2>0.05&orig.pvalue>0.05){
      return(1)
    }else if(drop.pvalue2<0.05&orig.pvalue<0.05){
      return(1)
    }else if(drop.pvalue2>0.05&orig.pvalue<0.05){
      return(0)
    }else if(drop.pvalue2<0.05&orig.pvalue>0.05){
      return(0)
    }
  }
  all$Consistent1 <- 1

  for(i in 1:length(celltypeall)){

    all$Consistent1[i] <- jugep(all$linear.P[i],all$p.value_Non_adjusted[i])
    all$Consistent2[i] <- jugep(all$logistic.P[i],all$adjusted.p.vaule_Benjamini_Hochberg[i])


  }



  rownames(all) <-  all$celltype

  all$celltype=NULL
  all$cellname=NULL
  all$celltypeall.i.=NULL
  all$logistic.P=NULL

  names(all) <- c("SPARKLE.pvalue","Wilcoxon.pvalue","Wilcoxon.adj.pvalue","Consistance1","Consistance2")

  all$inconsistency1 <- (1-mean(all$Consistance1))*100
  all$inconsistency2 <- (1-mean(all$Consistance2))*100
  all$groupname <- groupname
  all$cellname <- rownames(all)
  return(all)

}


